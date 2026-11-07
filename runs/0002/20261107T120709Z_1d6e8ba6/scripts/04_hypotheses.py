"""Pathway-level hypotheses for prompt_0002, each with the supporting test and a figure.

H1  Metabolism-subtype lines run a coordinated glutathione/antioxidant + xenobiotic
    detox programme (co-regulated with OXPHOS) -> predicted redox/ferroptosis vulnerability.
H2  MYC-target protein load is tightly coupled to the proteostasis (UPR) and DNA-repair
    machinery across the panel, independent of subtype -> proteostasis vulnerability.
H3  The panel is not organised along a Warburg glycolysis<->OXPHOS trade-off; a subset is
    OXPHOS-high -> candidate OXPHOS-inhibitor responders.
STRING functional enrichment is used as an orthogonal check on the detox gene set.
"""
import json
import urllib.parse
import urllib.request

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, ttest_ind

from common import (CAT, CELL_META, DATA, INK2, SUBTYPE_COLORS, TAB, diverging_cmap, read_gmt,
                    savefig, style)

style()
out = {}
L = pd.read_csv(TAB / "protein_matrix_cellline_mean_log2.csv.gz", index_col=0)
gene = L.pop("gene").fillna(""); sym = gene.str.split(";").str[0]
lines = list(L.columns); sub = CELL_META.paper_subtype.reindex(lines)
Lc = L[L.notna().sum(axis=1) >= 12]
Z = Lc.sub(Lc.mean(axis=1), axis=0).div(Lc.std(axis=1), axis=0)
Z.index = sym.reindex(Lc.index); Z = Z[~Z.index.duplicated()]
hall = read_gmt(DATA / "genesets/h.all.v2023.1.Hs.symbols.gmt", set(Z.index), 12, 400)
def hs(k): return Z.loc[list(hall[k])].mean()
ORDER = ["AGS", "eAGS", "AZ521", "SNU1750", "MKN28", "N87", "TMK1",
         "NCC24", "SNU719", "YCCEL1", "KatoIII", "Ocum1", "SNU601", "YCC3"]

# ======================================================================= H1
DETOX = {"glutathione synthesis": ["GCLC", "GCLM", "GSS"],
         "GSH peroxidases/transferases": ["GPX1", "GPX4", "GSTP1", "GSTK1", "GSTM3", "GSTO1", "MGST1", "MGST2"],
         "aldehyde/ketone reduction": ["ALDH3A1", "ALDH3A2", "ALDH4A1", "ALDH1A1", "AKR1A1", "AKR1B1", "AKR1C1"],
         "NADPH/antioxidant": ["NQO1", "TXN", "TXNRD1", "PRDX1", "PRDX6", "G6PD"],
         "efflux": ["ABCC1", "ABCC3"]}
detox_genes = [g for v in DETOX.values() for g in v if g in Z.index]
ds = Z.loc[detox_genes].mean()
a, b = ds[sub == "EMT"], ds[sub == "Metabolism"]
r_do, p_do = pearsonr(ds[lines], hs("HALLMARK_OXIDATIVE_PHOSPHORYLATION")[lines])
out["H1_detox"] = {"n_proteins": len(detox_genes), "diff_EMT_minus_Metab": float(a.mean() - b.mean()),
                   "welch_p": float(ttest_ind(a, b, equal_var=False)[1]),
                   "detox_vs_OXPHOS_r": float(r_do), "detox_vs_OXPHOS_p": float(p_do),
                   "proteins": detox_genes}

# STRING functional enrichment of the detox gene set (orthogonal check)
try:
    data = urllib.parse.urlencode({"identifiers": "\r".join(detox_genes), "species": 9606,
                                   "caller_identity": "agentic_proteomics_research"}).encode()
    req = urllib.request.Request("https://string-db.org/api/json/enrichment", data=data)
    enr = json.load(urllib.request.urlopen(req, timeout=60))
    enr = pd.DataFrame(enr)
    enr = enr[enr.category.isin(["KEGG", "RCTM", "Process", "Hallmark", "WikiPathways"])]
    enr = enr.sort_values("fdr")[["category", "term", "description", "number_of_genes", "fdr"]]
    enr.to_csv(TAB / "H1_string_enrichment_detox.csv", index=False)
    out["H1_string_top"] = enr.head(10)[["category", "description", "fdr"]].values.tolist()
except Exception as e:
    out["H1_string_error"] = str(e)

# ======================================================================= H2
myc = (hs("HALLMARK_MYC_TARGETS_V1") + hs("HALLMARK_MYC_TARGETS_V2")) / 2
upr = hs("HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
dnar = hs("HALLMARK_DNA_REPAIR")
out["H2_MYC_proteostasis"] = {"MYC_vs_UPR_r": float(pearsonr(myc[lines], upr[lines])[0]),
                              "MYC_vs_UPR_p": float(pearsonr(myc[lines], upr[lines])[1]),
                              "MYC_vs_DNArepair_r": float(pearsonr(myc[lines], dnar[lines])[0]),
                              "MYC_vs_DNArepair_p": float(pearsonr(myc[lines], dnar[lines])[1]),
                              "MYC_high_lines": myc.sort_values(ascending=False).head(4).round(2).to_dict()}

# ======================================================================= H3
ox = hs("HALLMARK_OXIDATIVE_PHOSPHORYLATION"); gly = hs("HALLMARK_GLYCOLYSIS")
out["H3_oxphos_glycolysis"] = {"OXPHOS_vs_GLYCOLYSIS_r": float(pearsonr(ox[lines], gly[lines])[0]),
                               "OXPHOS_vs_GLYCOLYSIS_p": float(pearsonr(ox[lines], gly[lines])[1]),
                               "OXPHOS_high_lines": ox.sort_values(ascending=False).head(4).round(2).to_dict()}
print(json.dumps(out, indent=2, default=float))
json.dump(out, open(TAB / "hypotheses_summary.json", "w"), indent=2, default=float)

# ---------------------------------------------------------------- figure H1
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={"wspace": 0.3, "width_ratios": [1.7, 1]})
ax = axes[0]
rows, ylab, grp = [], [], []
for cat, gl in DETOX.items():
    for g in gl:
        if g in Z.index:
            rows.append(Z.loc[g, ORDER]); ylab.append(g); grp.append(cat)
M = pd.DataFrame(rows, index=ylab)
im = ax.imshow(M.values, cmap=diverging_cmap(), vmin=-2, vmax=2, aspect="auto")
ax.set_xticks(range(len(ORDER)), ORDER, rotation=90)
for t in ax.get_xticklabels():
    t.set_color(SUBTYPE_COLORS[sub[t.get_text()]])
ax.set_yticks(range(len(ylab)), ylab, fontsize=7)
prev = None
for i, gname in enumerate(grp):
    if gname != prev:
        ax.axhline(i - 0.5, color="white", lw=2); prev = gname
fig.colorbar(im, ax=ax, shrink=0.6, label="protein z-score across lines")
ax.set_title("H1  Glutathione / antioxidant / xenobiotic-detox proteins\n(label colour = Hu et al. subtype)", fontsize=10)
ax = axes[1]
for s, c in SUBTYPE_COLORS.items():
    m = (sub == s).values
    ax.scatter(hs("HALLMARK_OXIDATIVE_PHOSPHORYLATION")[lines][m], ds[lines][m], s=42, color=c, label=s, edgecolor="white")
for ln in lines:
    ax.annotate(ln, (hs("HALLMARK_OXIDATIVE_PHOSPHORYLATION")[ln], ds[ln]), fontsize=7, color=INK2,
                xytext=(3, 3), textcoords="offset points")
ax.set_xlabel("OXPHOS activity (protein z)")
ax.set_ylabel("detox/antioxidant programme (protein z)")
ax.set_title(f"detox vs OXPHOS: r = {r_do:.2f} (p = {p_do:.1g})\nMetab > EMT: Welch p = {out['H1_detox']['welch_p']:.1g}", fontsize=10)
ax.legend(title="Hu et al.", fontsize=7)
savefig(fig, "05_H1_detox_programme")
plt.close(fig)

# ---------------------------------------------------------------- figure H2/H3
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), gridspec_kw={"wspace": 0.34})
def scat(ax, x, y, xl, yl, title):
    for s, c in SUBTYPE_COLORS.items():
        m = (sub == s).values
        ax.scatter(x[lines][m], y[lines][m], s=40, color=c, label=s, edgecolor="white")
    for ln in lines:
        ax.annotate(ln, (x[ln], y[ln]), fontsize=6.5, color=INK2, xytext=(3, 3), textcoords="offset points")
    r, p = pearsonr(x[lines], y[lines])
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(f"{title}\nr = {r:.2f}, p = {p:.1g}", fontsize=9.5)
    return r, p
scat(axes[0], myc, upr, "MYC-target activity", "UPR / proteostasis activity", "H2  MYC ↔ proteostasis")
scat(axes[1], myc, dnar, "MYC-target activity", "DNA-repair activity", "H2  MYC ↔ DNA repair")
scat(axes[2], gly, ox, "glycolysis activity", "OXPHOS activity", "H3  glycolysis vs OXPHOS\n(no Warburg trade-off)")
axes[0].legend(title="Hu et al.", fontsize=7)
savefig(fig, "06_H2_H3_metabolic_coupling")
plt.close(fig)
