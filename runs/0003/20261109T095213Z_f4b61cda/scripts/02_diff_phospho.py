"""Differential phosphorylation between the proteome subtypes (and EBV+ vs EBV-), on the
phosphoproteome (PP, S/T/Y) and tyrosine (YP) layers. Phosphosites are analysed both raw and
protein-normalised (phospho - parent protein), to separate signalling from abundance change.
Writes seqrnk files (flanking-sequence ranked lists) for phosX kinase-activity inference.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CAT, CELL_META, INK2, TAB, savefig, style
from phospho_common import load_PP, load_WP_linemean, load_YP, write_seqrnk
from stats import moderated_t

style()
SEQ = TAB.parent / "seqrnk"; SEQ.mkdir(exist_ok=True)
order = list(CELL_META.index)
ppm, ppx = load_PP(order)
ypm, ypx = load_YP(order)
wp_acc, _ = load_WP_linemean(order)

EMT = [l for l in order if CELL_META.paper_subtype[l] == "EMT"]
MET = [l for l in order if CELL_META.paper_subtype[l] == "Metabolism"]
EBV = ["NCC24", "SNU719", "YCCEL1"]
NONEBV = [l for l in order if l not in EBV + ["eAGS"]]
CONTRASTS = {"EMT_vs_Metabolism": (EMT, MET, 4), "EBV_vs_nonEBV": (EBV, NONEBV, 3)}

# protein level per PP site (for normalisation)
prot = wp_acc.reindex(ppm.UniprotID.values); prot.index = ppx.index

summary = {}
for name, (g1, g2, mn) in CONTRASTS.items():
    # ---- PP (S/T/Y global phosphoproteome) ----
    dpp = moderated_t(ppx, g1, g2, min_per_group=mn)
    dpp = dpp.join(ppm[["site", "GeneSymbol", "restype", "win10", "UniprotID"]])
    # protein-normalised logFC where parent protein present
    nn = (ppx - prot)
    dnn = moderated_t(nn, g1, g2, min_per_group=mn)[["logFC", "q"]].rename(
        columns={"logFC": "logFC_protnorm", "q": "q_protnorm"})
    dpp = dpp.join(dnn)
    dpp.sort_values("p").to_csv(TAB / f"diff_PP_{name}.csv.gz", float_format="%.5g", index=False)
    # ---- YP (tyrosine IP) ----
    dyp = moderated_t(ypx, g1, g2, min_per_group=(2 if name == "EBV_vs_nonEBV" else 2))
    dyp = dyp.join(ypm[["site", "GeneSymbol", "restype", "win10", "UniprotID"]])
    dyp.sort_values("p").to_csv(TAB / f"diff_YP_{name}.csv.gz", float_format="%.5g", index=False)
    # ---- seqrnk files for phosX ----
    st = dpp[dpp.restype.isin(["S", "T"]) & dpp.win10.notna()]
    n_st = write_seqrnk(zip(st.win10, st.logFC), SEQ / f"{name}_ST.seqrnk")
    yrows = pd.concat([dpp[dpp.restype == "Y"][["win10", "logFC"]], dyp[["win10", "logFC"]]]).dropna()
    n_y = write_seqrnk(zip(yrows.win10, yrows.logFC), SEQ / f"{name}_Y.seqrnk")
    summary[name] = {"PP_tested": int(len(dpp)), "PP_q<0.05": int((dpp.q < 0.05).sum()),
                     "PP_protnorm_q<0.05": int((dpp.q_protnorm < 0.05).sum()),
                     "YP_tested": int(len(dyp)), "YP_q<0.05": int((dyp.q < 0.05).sum()),
                     "seqrnk_ST": n_st, "seqrnk_Y": n_y,
                     "PP_top_up_g1": dpp[dpp.q < 0.05].sort_values("logFC", ascending=False).site.head(12).tolist(),
                     "PP_top_up_g2": dpp[dpp.q < 0.05].sort_values("logFC").site.head(12).tolist(),
                     "YP_top_up_g1": dyp[dyp.q < 0.1].sort_values("logFC", ascending=False).site.head(10).tolist(),
                     "YP_top_up_g2": dyp[dyp.q < 0.1].sort_values("logFC").site.head(10).tolist()}
    print("=" * 20, name, summary[name]["PP_q<0.05"], "PP /", summary[name]["YP_q<0.05"], "YP sig")

json.dump(summary, open(TAB / "diff_phospho_summary.json", "w"), indent=2, default=float)
print(json.dumps(summary, indent=2, default=float))


# volcano figure for the primary contrast (PP protein-normalised + YP)
def repel(ax, items, fontsize=6):
    ax.figure.canvas.draw()
    inv = ax.transData.inverted(); placed = []
    for x, y, txt in sorted(items, key=lambda t: -t[1]):
        xd, yd = ax.transData.transform((x, y))
        for k in range(44):
            cand = yd + 12.0 * ((k + 1) // 2) * (1 if k % 2 else -1)
            if all(abs(cand - py) > 11 or abs(xd - px) > 60 for px, py in placed):
                yd = cand; break
        placed.append((xd, yd)); xa, ya = inv.transform((xd, yd))
        ax.annotate(txt, (x, y), xytext=(xa, ya), textcoords="data", fontsize=fontsize, color=INK2,
                    arrowprops=dict(arrowstyle="-", color="#c3c2b7", lw=0.4, shrinkA=0, shrinkB=2))


def volcano(ax, d, lfc, q, title, extra=()):
    sig = d[q] < 0.05
    ax.scatter(d[lfc][~sig], -np.log10(d.p[~sig]), s=4, color="#c3c2b7", linewidths=0)
    ax.scatter(d[lfc][sig & (d[lfc] > 0)], -np.log10(d.p[sig & (d[lfc] > 0)]), s=7, color=CAT[1], linewidths=0)
    ax.scatter(d[lfc][sig & (d[lfc] < 0)], -np.log10(d.p[sig & (d[lfc] < 0)]), s=7, color=CAT[0], linewidths=0)
    lab = d[sig].assign(sc=lambda x: -np.log10(x.p) * x[lfc].abs()).sort_values("sc", ascending=False).head(10)
    lab = pd.concat([lab, d[d.site.isin(extra)]]).drop_duplicates("site")
    repel(ax, [(r[lfc], -np.log10(r.p), r.site) for _, r in lab.iterrows()])
    ax.axvline(0, color="#c3c2b7", lw=0.6); ax.set_xlabel(lfc); ax.set_ylabel("-log10 p")
    ax.set_title(title, fontsize=9)


dpp = pd.read_csv(TAB / "diff_PP_EMT_vs_Metabolism.csv.gz")
dyp = pd.read_csv(TAB / "diff_YP_EMT_vs_Metabolism.csv.gz")
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"wspace": 0.28})
volcano(axes[0], dpp.dropna(subset=["logFC_protnorm", "q_protnorm"]), "logFC_protnorm", "q_protnorm",
        "A  Phosphoproteome, protein-normalised\n(EMT - Metabolism)")
volcano(axes[1], dyp, "logFC", "q", "B  Tyrosine phosphoproteome\n(EMT - Metabolism)")
fig.text(0.5, -0.02, "positive = higher in EMT subtype; negative = higher in Metabolism subtype",
         ha="center", fontsize=8, color=INK2)
savefig(fig, "02_diff_phospho_volcano")
plt.close(fig)
