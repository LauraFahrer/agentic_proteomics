"""Signaling hypotheses for prompt_0003: assemble the evidence into a signaling-module
figure and a machine-readable summary.

H1  Metabolism/GI-lineage subtype = convergent RTK network (MET, FGFR2/3, RON, Eph, HER2)
    -> RAS-ERK via SHC1/GAREM1/RASA1 adaptor pY -> MEK dependency + RTK redundancy.
H2  EMT subtype = PI3K-AKT-mTOR / SGK (AGC) Ser/Thr signaling -> mTOR/SGK dependency.
H3  The robust subtype signaling difference is tyrosine-level; the bulk Ser/Thr phosphoproteome
    difference is abundance-coupled (0 protein-normalised sites at q<0.05).
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from common import CAT, CELL_META, INK2, SUBTYPE_COLORS, TAB, diverging_cmap, savefig, style

style()
out = {}
order = ["AGS", "eAGS", "AZ521", "SNU1750", "MKN28", "N87", "TMK1",
         "NCC24", "SNU719", "YCCEL1", "KatoIII", "Ocum1", "SNU601", "YCC3"]
sub = CELL_META.paper_subtype
EMT = [l for l in order if sub[l] == "EMT"]; MET = [l for l in order if sub[l] == "Metabolism"]

yp = pd.read_csv(TAB / "yp_linemean.csv.gz").set_index("site")
dyp = pd.read_csv(TAB / "diff_YP_EMT_vs_Metabolism.csv.gz")

# --- H1: RTK -> adaptor -> RAS-ERK tyrosine module (significant pY sites) ---
# map the significant pY sites to readable gene_site labels using the diff table's gene column
sig = dyp[dyp.q < 0.05].copy()
ACC2GENE = {"P29323": "EPHB2", "P29353": "SHC1", "P20936": "RASA1", "P16144": "ITGB4",
            "Q9H706": "GAREM1", "O00560": "SDCBP", "Q8WWI1": "LMO7"}
def nice(site, gene):
    if isinstance(gene, str) and gene != "nan":
        return f"{gene}_{site.split('_')[-1]}"
    acc = site.split("_")[0].split("-")[0]
    return f"{ACC2GENE.get(acc, acc)}_{site.split('_')[-1]}"
sig["label"] = [nice(s, g) for s, g in zip(sig.site, sig.GeneSymbol)]
# flag RTK/adaptor module members
MODULE = ["EPHB2", "SHC1", "RASA1", "ITGB4", "GAREM1", "SHB", "PEBP1", "SASH1", "SDCBP"]
sig["module"] = sig.label.str.split("_").str[0].isin(MODULE)
out["H1_sig_pY_sites"] = int(len(sig))
out["H1_sig_pY_all_higher_in_Metabolism"] = bool((sig.logFC < 0).mean() > 0.9)
out["H1_adaptor_module_sites"] = sig[sig.module].label.tolist()

# Tyr kinase activity (already computed by phosX)
yk = pd.read_csv(TAB / "kinase/EMT_vs_Metabolism_Y_kinases.tsv", sep="\t", index_col=0)
rtk = ["MET", "FGFR2", "FGFR3", "MST1R", "EPHB1", "EPHA3", "EPHB3", "EPHA8", "HER2", "TRKB", "TRKA", "JAK3", "DDR2", "CSFR"]
out["H1_RTK_activity_higher_in_Metab_FDR"] = {k: float(yk.loc[k, "FDR q value"]) for k in rtk if k in yk.index}

# --- H2: AGC (AKT/SGK/mTOR) Ser/Thr activity higher in EMT ---
stk = pd.read_csv(TAB / "kinase/EMT_vs_Metabolism_ST_kinases.tsv", sep="\t", index_col=0)
agc = ["AKT1", "AKT2", "AKT3", "SGK1", "SGK3", "P70S6K", "RSK3", "P90RSK"]
mapk = ["MEK5", "MEK2", "MEK1", "ZAK", "MEKK2", "KHS1", "KHS2", "GCK"]
out["H2_AGC_higher_in_EMT_KS"] = {k: float(stk.loc[k, "KS"]) for k in agc if k in stk.index}
out["H2_MAPK_higher_in_Metab_KS"] = {k: float(stk.loc[k, "KS"]) for k in mapk if k in stk.index}

# --- H3: protein-normalised Ser/Thr has no subtype difference ---
diffsum = json.load(open(TAB / "diff_phospho_summary.json"))
out["H3_PP_q<0.05_raw"] = diffsum["EMT_vs_Metabolism"]["PP_q<0.05"]
out["H3_PP_q<0.05_protnorm"] = diffsum["EMT_vs_Metabolism"]["PP_protnorm_q<0.05"]
out["H3_YP_q<0.05"] = diffsum["EMT_vs_Metabolism"]["YP_q<0.05"]
json.dump(out, open(TAB / "hypotheses_summary.json", "w"), indent=2, default=float)
print(json.dumps(out, indent=2, default=float))

# ---- figure: signaling module heatmaps + kinase-axis summary ----
fig = plt.figure(figsize=(13, 6.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1.2, 1.1], wspace=0.5)

# (A) key differential pY sites across lines (z-scored)
ax = fig.add_subplot(gs[0, 0])
sites = sig.sort_values("logFC").site.tolist()
M = yp.loc[[s for s in sites if s in yp.index], order]
Mz = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1), axis=0)
lab = [nice(s, dyp.set_index("site").GeneSymbol.get(s)) for s in M.index]
im = ax.imshow(Mz.values, cmap=diverging_cmap(), vmin=-1.8, vmax=1.8, aspect="auto")
ax.set_xticks(range(len(order)), order, rotation=90, fontsize=7)
for t in ax.get_xticklabels():
    t.set_color(SUBTYPE_COLORS[sub[t.get_text()]])
ax.set_yticks(range(len(lab)), lab, fontsize=6)
ax.set_title("A  Differential pY sites (q<0.05)\nacross lines (RTK/adaptor module in bold)", fontsize=9)
for i, l in enumerate(lab):
    if l.split("_")[0] in MODULE:
        ax.get_yticklabels()[i].set_fontweight("bold")
fig.colorbar(im, ax=ax, shrink=0.5, label="pY z-score")

# (B) Tyr kinase activity (RTKs higher in Metab)
ax = fig.add_subplot(gs[0, 1])
yk2 = yk.sort_values("KS")
selk = pd.concat([yk2.head(14), yk2.tail(3)])
ax.barh(range(len(selk)), selk.KS, color=[CAT[1] if v > 0 else CAT[0] for v in selk.KS], height=0.72)
ax.set_yticks(range(len(selk)), selk.index, fontsize=7)
ax.axvline(0, color="#c3c2b7", lw=0.6)
for i, (v, q) in enumerate(zip(selk.KS, selk["FDR q value"])):
    if q < 0.1:
        ax.text(v - 0.01, i, "*", va="center", ha="right", color=INK2, fontsize=11)
ax.set_xlabel("Tyr-kinase enrichment (KS)", fontsize=8)
ax.set_title("B  Tyrosine kinome\n(− = higher in Metabolism)", fontsize=9)

# (C) AKT/SGK vs MEK/MAP3K S/T axis
ax = fig.add_subplot(gs[0, 2])
rows = agc + ["|"] + mapk
vals = [stk.loc[k, "KS"] if k in stk.index else np.nan for k in agc] + [0] + \
       [stk.loc[k, "KS"] if k in stk.index else np.nan for k in mapk]
cols = [CAT[1]] * len(agc) + ["white"] + [CAT[0]] * len(mapk)
ax.barh(range(len(rows)), vals, color=cols, height=0.72)
ax.set_yticks(range(len(rows)), rows, fontsize=7)
ax.axvline(0, color="#c3c2b7", lw=0.6)
ax.set_xlabel("Ser/Thr-kinase enrichment (KS)", fontsize=8)
ax.set_title("C  AGC (EMT, +) vs\nMAP2K/MAP3K (Metab, −)", fontsize=9)
fig.suptitle("Subtype-specific signalling: EMT = PI3K-AKT-mTOR/SGK; Metabolism = RTK → RAS-ERK (MEK)",
             fontsize=11, y=1.02)
savefig(fig, "04_signaling_modules")
plt.close(fig)
