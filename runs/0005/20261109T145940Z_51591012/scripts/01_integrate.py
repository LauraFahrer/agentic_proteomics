"""Cross-run synthesis: build a per-line 'axis scorecard' from the derived signatures of
runs 0001-0004, quantify how strongly the analyses converge, and show the convergence is
NOT independent validation (all axes are computed from the same protein-abundance matrix).
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, ttest_ind

from common import (CAT, FIG, INK2, ORDER, RUN0002, RUN0003, RUN0004, SUBTYPE, SUBTYPE_COLORS,
                    TAB, diverging_cmap, read_gmt, savefig, style)

style()
out = {}
DATA = __import__("os").environ.get("PXD048538_DATA", "data_cache")
GS = __import__("pathlib").Path(DATA) / "genesets"

# shared protein matrix (committed by run 0002)
L = pd.read_csv(RUN0002 / "results/tables/protein_matrix_cellline_mean_log2.csv.gz", index_col=0)
gene = L.pop("gene").fillna("")
sym = gene.str.split(";").str[0]
lines = [c for c in ORDER if c in L.columns]
sub = pd.Series({l: SUBTYPE[l] for l in lines})
Lc = L[L.notna().sum(axis=1) >= 12]
Z = Lc.sub(Lc.mean(axis=1), axis=0).div(Lc.std(axis=1) + 1e-9, axis=0)
Z.index = sym.reindex(Lc.index)
Z = Z[~Z.index.duplicated()]
hall = read_gmt(GS / "h.all.v2023.1.Hs.symbols.gmt", set(Z.index), 12, 400)
def hs(k): return Z.loc[list(hall[k])].mean()
def gset(genes):
    g = [x for x in genes if x in Z.index]
    return Z.loc[g].mean()

Lg = L.copy(); Lg.index = sym.values; Lg = Lg.groupby(level=0).mean()
floor = np.nanmin(Lg.values) - 1
def gset_floor(genes):
    g = [x for x in genes if x in Lg.index]
    M = Lg.loc[g].fillna(floor)
    return M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1) + 1e-9, axis=0).mean()

# --- per-line axis scores, each tied to the run(s) that derived it ---
axes_scores = pd.DataFrame({
    "GI-lineage TF (HNF1A/4A/FOXA2) [r1]": gset_floor(["HNF1A", "HNF4A", "FOXA2"]),
    "Oxidative metabolism (OXPHOS+FAO) [r2,r4]": (hs("HALLMARK_OXIDATIVE_PHOSPHORYLATION") + hs("HALLMARK_FATTY_ACID_METABOLISM")) / 2,
    "Glutathione/antioxidant detox [r2]": gset(["GCLC", "GCLM", "GPX1", "GPX4", "GSTP1", "GSTK1", "AKR1A1", "TXNRD1", "G6PD", "ABCC1"]),
    "Cytoplasmic translation/ribosome [r2,r3]": gset([g for g in Z.index if g.startswith(("RPL", "RPS", "EIF"))]),
    "RTK-adaptor pY module (EPHB2/SHC1/RASA1) [r3]": gset(["EPHB2", "SHC1", "GAREM1", "SHB", "RASA1", "ITGB4"]),
    "OXPHOS complexes (resp chain/cytochrome) [r4]": gset([g for g in Z.index if g.startswith(("NDUF", "COX", "UQCR", "SDH", "ATP5"))]),
    "MYC targets / proliferation [r2,r3]": (hs("HALLMARK_MYC_TARGETS_V1") + hs("HALLMARK_E2F_TARGETS")) / 2,
    "Proteostasis / UPR [r2,r3]": hs("HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
})
A = axes_scores.T[lines]
A.to_csv(TAB / "per_line_axis_scores.csv", float_format="%.3f")

# subtype association of each axis (the metabolic axes should split EMT vs Metab)
rows = []
for name in A.index:
    a, b = A.loc[name, sub[sub == "EMT"].index], A.loc[name, sub[sub == "Metabolism"].index]
    t, p = ttest_ind(a, b, equal_var=False)
    rows.append({"axis": name, "EMT_minus_Metab": float(a.mean() - b.mean()), "welch_p": float(p)})
axis_sub = pd.DataFrame(rows)
axis_sub.to_csv(TAB / "axis_subtype_association.csv", index=False)
out["axis_subtype"] = axis_sub.round(3).to_dict("records")

# convergence: correlation among the axes across the 14 lines
corr = np.corrcoef(A.values)
corr = pd.DataFrame(corr, index=A.index, columns=A.index)
# the metabolic-axis cluster (lineage, oxidative, detox, OXPHOS complexes) vs translation/MYC
metab_axes = [a for a in A.index if any(k in a for k in ["lineage", "Oxidative", "detox", "OXPHOS complexes"])]
prolif_axes = [a for a in A.index if any(k in a for k in ["translation", "MYC", "Proteostasis"])]
mm = corr.loc[metab_axes, metab_axes].values
out["mean_r_within_metabolic_axes"] = float(mm[np.triu_indices(len(metab_axes), 1)].mean())
out["mean_r_metabolic_vs_proliferation"] = float(corr.loc[metab_axes, prolif_axes].values.mean())
json.dump(out, open(TAB / "integration_summary.json", "w"), indent=2, default=float)
print(json.dumps(out, indent=2, default=float))

# ---- figure 1: per-line scorecard + axis convergence ----
fig = plt.figure(figsize=(13.5, 6.0))
gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1], wspace=0.55)
ax = fig.add_subplot(gs[0, 0])
col = list(sub.sort_values().index)
im = ax.imshow(A[col].values, cmap=diverging_cmap(), vmin=-1.3, vmax=1.3, aspect="auto")
ax.set_xticks(range(len(col)), col, rotation=90)
for t in ax.get_xticklabels():
    t.set_color(SUBTYPE_COLORS[sub[t.get_text()]])
ax.set_yticks(range(len(A)), A.index, fontsize=7.5)
fig.colorbar(im, ax=ax, shrink=0.6, label="per-line z-score")
ax.set_title("A  Cross-run axis scorecard (14 lines)\n[r1-r4 = which run derived the axis]", fontsize=9.5)
ax = fig.add_subplot(gs[0, 1])
order_ax = metab_axes + prolif_axes + [a for a in A.index if a not in metab_axes + prolif_axes]
Ccc = corr.loc[order_ax, order_ax]
im = ax.imshow(Ccc.values, cmap=diverging_cmap(), vmin=-1, vmax=1)
short = [a.split(" [")[0][:22] for a in order_ax]
ax.set_xticks(range(len(short)), short, rotation=90, fontsize=6.5)
ax.set_yticks(range(len(short)), short, fontsize=6.5)
fig.colorbar(im, ax=ax, shrink=0.6, label="Pearson r across 14 lines")
ax.set_title("B  Axes co-vary (same data!)\nmetabolic axes mean r = %.2f" % out["mean_r_within_metabolic_axes"], fontsize=9.5)
savefig(fig, "01_cross_run_scorecard")
plt.close(fig)
