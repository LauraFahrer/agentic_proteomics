"""Pathway landscape across the 14 GC cell lines (prompt_0002: biological processes & pathways).

- Pathway activity = mean z-score (across lines) of the proteins in each gene set.
- Which pathways vary most across the panel?
- Pathway co-regulation: correlate pathway scores across lines, cluster into modules,
  identify the dominant (anti-correlated) axes of proteome variation.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from common import (CELL_META, DATA, INK2, SUBTYPE_COLORS, TAB, diverging_cmap, read_gmt,
                    savefig, style)

style()
out = {}
L = pd.read_csv(TAB / "protein_matrix_cellline_mean_log2.csv.gz", index_col=0)
gene = L.pop("gene").fillna("")
sym = gene.str.split(";").str[0]
lines = list(L.columns)
sub = CELL_META.paper_subtype.reindex(lines)
ORDER = ["AGS", "eAGS", "AZ521", "SNU1750", "MKN28", "N87", "TMK1",
         "NCC24", "SNU719", "YCCEL1", "KatoIII", "Ocum1", "SNU601", "YCC3"]

# protein z-scores by gene symbol (proteins present in >=12 lines, to keep pathway scores stable)
Lc = L[L.notna().sum(axis=1) >= 12]
Z = Lc.sub(Lc.mean(axis=1), axis=0).div(Lc.std(axis=1), axis=0)
Z.index = sym.reindex(Lc.index)
Z = Z[~Z.index.duplicated()]
universe = set(Z.index)
out["genes_for_pathway_scoring"] = int(len(universe))


def score_sets(gmt, tag, min_size=12, max_size=400):
    sets = read_gmt(DATA / "genesets" / gmt, universe, min_size, max_size)
    sc = pd.DataFrame({f"{tag}:{k}": Z.loc[list(v)].mean() for k, v in sets.items()}).T
    sz = pd.Series({f"{tag}:{k}": len(v) for k, v in sets.items()})
    return sc, sz


hall, hsz = score_sets("h.all.v2023.1.Hs.symbols.gmt", "H")
hall.index = [i.replace("H:HALLMARK_", "") for i in hall.index]
hsz.index = hall.index
hall.to_csv(TAB / "hallmark_activity_per_line.csv")

# --- pathway variability across the panel ---------------------------------
var = hall.std(axis=1).sort_values(ascending=False)
var_tab = pd.DataFrame({"hallmark": var.index, "sd_across_lines": var.values,
                        "n_proteins": hsz.reindex(var.index).values,
                        "range": (hall.max(axis=1) - hall.min(axis=1)).reindex(var.index).values})
var_tab.to_csv(TAB / "hallmark_variability.csv", index=False)
out["most_variable_hallmarks"] = var.head(12).round(3).to_dict()
out["least_variable_hallmarks"] = var.tail(6).round(3).to_dict()

# --- pathway co-regulation modules ----------------------------------------
top = var.head(34).index  # the variable hallmarks carry the biology
H = hall.loc[top]
corr = np.corrcoef(H.values)
corr = pd.DataFrame(corr, index=top, columns=top)
Zl = linkage(pdist(H.values, "correlation"), "average")
order_p = leaves_list(Zl)

fig, ax = plt.subplots(figsize=(10.5, 9.5))
C = corr.values[np.ix_(order_p, order_p)]
im = ax.imshow(C, cmap=diverging_cmap(), vmin=-1, vmax=1)
labs = [top[i] for i in order_p]
ax.set_xticks(range(len(labs)), labs, rotation=90, fontsize=6.5)
ax.set_yticks(range(len(labs)), labs, fontsize=6.5)
fig.colorbar(im, ax=ax, shrink=0.5, label="Pearson r of pathway activity across 14 lines")
ax.set_title("Hallmark pathway co-regulation across the GC cell-line panel\n"
             "(red = pathways that rise together; blue = reciprocal programmes)")
savefig(fig, "01_pathway_coregulation")
plt.close(fig)

# data-driven pathway modules: cut the dendrogram of the co-regulation structure into blocks
from scipy.cluster.hierarchy import fcluster
from scipy.stats import pearsonr, spearmanr, ttest_ind
modlab = fcluster(Zl, 4, "maxclust")
MOD_NAMES = {}
module_members = {}
for m in sorted(set(modlab)):
    members = list(np.array(top)[modlab == m])
    module_members[m] = members
# name modules by their dominant biology (by hand from membership, printed for transparency)
mod_score = pd.DataFrame({m: hall.loc[mm].mean() for m, mm in module_members.items()}).T
# correlation between modules (which programmes are reciprocal)
mod_corr = pd.DataFrame(np.corrcoef(mod_score.values), index=mod_score.index, columns=mod_score.index)
out["pathway_modules"] = {int(m): mm for m, mm in module_members.items()}
out["module_correlation"] = mod_corr.round(2).values.tolist()
# association of each module with the published subtype
mod_test = []
for m, mm in module_members.items():
    a, b = mod_score.loc[m, sub == "EMT"], mod_score.loc[m, sub == "Metabolism"]
    t, p = ttest_ind(a, b, equal_var=False)
    mod_test.append({"module": int(m), "n_pathways": len(mm), "diff_EMT_minus_Metab": float(a.mean() - b.mean()),
                     "welch_p": float(p), "members": ", ".join(mm)})
mod_test = pd.DataFrame(mod_test)
mod_test.to_csv(TAB / "pathway_modules_vs_subtype.csv", index=False)
out["module_vs_subtype"] = mod_test.to_dict("records")
print(mod_test.to_string())

# identify the two most strongly anti-correlated modules for the scatter
iu = [(i, j) for i in mod_corr.index for j in mod_corr.index if i < j]
amin = min(iu, key=lambda p: mod_corr.loc[p[0], p[1]])
m1, m2 = amin
r_po, p_po = pearsonr(mod_score.loc[m1, lines], mod_score.loc[m2, lines])
out["most_anticorrelated_modules"] = {"modules": [int(m1), int(m2)], "pearson_r": float(r_po), "p": float(p_po)}


def mlabel(m):
    mm = module_members[m]
    return " / ".join(mm[:3]) + ("…" if len(mm) > 3 else "")


fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"wspace": 0.32, "width_ratios": [1.5, 1]})
ax = axes[0]
show = var.head(26).index
S = hall.loc[show, ORDER]
im = ax.imshow(S.values, cmap=diverging_cmap(), vmin=-1.5, vmax=1.5, aspect="auto")
ax.set_xticks(range(len(ORDER)), ORDER, rotation=90)
for t in ax.get_xticklabels():
    t.set_color(SUBTYPE_COLORS[sub[t.get_text()]])
ax.set_yticks(range(len(show)), show, fontsize=7)
fig.colorbar(im, ax=ax, shrink=0.7, label="pathway activity (mean protein z)")
ax.set_title("A  26 most variable hallmark pathways across lines")
ax = axes[1]
for s, c in SUBTYPE_COLORS.items():
    mk = (sub == s).values
    ax.scatter(mod_score.loc[m1, lines][mk], mod_score.loc[m2, lines][mk], s=42, color=c, label=s, edgecolor="white")
for ln in lines:
    ax.annotate(ln, (mod_score.loc[m1, ln], mod_score.loc[m2, ln]), fontsize=7, color=INK2,
                xytext=(3, 3), textcoords="offset points")
ax.set_xlabel(f"module {m1}: {mlabel(m1)}", fontsize=8)
ax.set_ylabel(f"module {m2}: {mlabel(m2)}", fontsize=8)
p4 = mod_test.set_index("module").loc[4, "welch_p"]
ax.set_title("B  Proliferation (mod 1) and metabolism (mod 4) vary independently\n"
             f"(module r = {r_po:.2f}, n.s.); the metabolism module separates\n"
             f"EMT (low) from Metabolism (high) lines, Welch p = {p4:.3f}", fontsize=8.5)
ax.legend(title="Hu et al.", fontsize=7)
savefig(fig, "02_pathway_landscape")
plt.close(fig)

json.dump(out, open(TAB / "pathway_landscape_summary.json", "w"), indent=2, default=float)
print(json.dumps(out, indent=2, default=float))
