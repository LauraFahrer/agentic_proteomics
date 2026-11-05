"""Differential abundance and rank-based pathway analysis for the main contrasts.

Line-level contrasts use the cell-line mean matrix (one value per line = the
independent biological unit; replicates are cultures of the same line).
The isogenic AGS vs eAGS contrast uses the 3 biological replicates per line.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (CAT, CELL_META, DATA, INK2, MUTED, TAB, read_gmt, savefig, style)
from stats import moderated_t, rank_set_test

style()
X = pd.read_csv(TAB / "protein_matrix_log2.csv.gz", index_col=0)
gene = X.pop("gene").fillna("")
L = pd.read_csv(TAB / "protein_matrix_cellline_mean_log2.csv.gz", index_col=0)
L.pop("gene")
sym = gene.str.split(";").str[0]
lines = list(L.columns)

EMT = [l for l in lines if CELL_META.paper_subtype[l] == "EMT"]
MET = [l for l in lines if CELL_META.paper_subtype[l] == "Metabolism"]
EBV_NAT = ["NCC24", "SNU719", "YCCEL1"]
EBV_NEG = [l for l in lines if l not in EBV_NAT + ["eAGS"]]
CL_EBVLIKE = ["TMK1", "NCC24", "SNU719", "YCCEL1"]
CL_DIFFUSE = ["KatoIII", "Ocum1", "SNU601", "YCC3"]
SRCC = ["KatoIII", "NCC24", "Ocum1", "SNU601"]

contrasts = {
    "paperEMT_vs_paperMetabolism": (L, EMT, MET, 3),
    "EBVnatural_vs_EBVneg": (L, EBV_NAT, EBV_NEG, 3),
    "clusterEBVlike_vs_rest": (L, CL_EBVLIKE, [l for l in lines if l not in CL_EBVLIKE], 3),
    "clusterDiffuse_vs_rest": (L, CL_DIFFUSE, [l for l in lines if l not in CL_DIFFUSE], 3),
    "SRCC_vs_nonSRCC": (L, SRCC, [l for l in lines if l not in SRCC], 3),
    "TMK1_vs_EBVneg_other": (None, None, None, None),  # handled separately below
    "eAGS_vs_AGS_replicates": (X, [f"eAGS_R{i}" for i in (1, 2, 3)],
                               [f"AGS_R{i}" for i in (1, 2, 3)], 2),
}

universe = set(sym[sym != ""])
gsets = {}
for f, tag in [("h.all.v2023.1.Hs.symbols.gmt", "H"), ("c2.cp.kegg.v7.5.1.symbols.gmt", "KEGG"),
               ("c2.cp.reactome.v7.5.1.symbols.gmt", "REACTOME"), ("c5.go.bp.v7.5.1.symbols.gmt", "GOBP"),
               ("c5.go.cc.v7.5.1.symbols.gmt", "GOCC")]:
    for k, v in read_gmt(DATA / "genesets" / f, universe, min_size=10, max_size=300).items():
        gsets[k] = v

summary = {}
results = {}
for name, (M, g1, g2, mn) in contrasts.items():
    if M is None:
        continue
    r = moderated_t(M, g1, g2, min_per_group=mn)
    r.insert(0, "gene", gene.reindex(r.index))
    r = r.sort_values("p")
    r.to_csv(TAB / f"diff_{name}.csv.gz", float_format="%.5g")
    results[name] = r
    # gene-level statistic: mean t over protein groups of the same leading symbol
    gt = r.assign(sym=sym.reindex(r.index)).query("sym != ''").groupby("sym").t.mean()
    gs = rank_set_test(gt, gsets)
    gs.to_csv(TAB / f"geneset_{name}.csv", index=False, float_format="%.4g")
    summary[name] = {
        "group1": g1, "group2": g2, "proteins_tested": int(len(r)), "prior_df": float(r.attrs["d0"]),
        "n_q<0.05": int((r.q < 0.05).sum()), "n_q<0.05_|logFC|>1": int(((r.q < 0.05) & (r.logFC.abs() > 1)).sum()),
        "top_up": r[(r.q < 0.05) & (r.logFC > 0)].head(15).gene.tolist(),
        "top_down": r[(r.q < 0.05) & (r.logFC < 0)].head(15).gene.tolist(),
        "genesets_up_q<0.05": gs[(gs.q < 0.05) & (gs.z > 0)].head(12).set.tolist(),
        "genesets_down_q<0.05": gs[(gs.q < 0.05) & (gs.z < 0)].head(12).set.tolist(),
    }
    print("=" * 30, name, summary[name]["n_q<0.05"], summary[name]["n_q<0.05_|logFC|>1"])
    print("UP", summary[name]["top_up"])
    print("DN", summary[name]["top_down"])
    print("GS+", summary[name]["genesets_up_q<0.05"])
    print("GS-", summary[name]["genesets_down_q<0.05"])

# ---- presence / absence (on-off) proteins for the EBV-like cluster and natural EBV lines
pres = L.notna()
onoff = []
for name, grp in [("EBVnatural", EBV_NAT), ("clusterEBVlike", CL_EBVLIKE),
                  ("clusterDiffuse", CL_DIFFUSE), ("paperEMT", EMT), ("paperMetabolism", MET)]:
    rest = [l for l in lines if l not in grp and not (name == "EBVnatural" and l == "eAGS")]
    only = pres[grp].all(axis=1) & ~pres[rest].any(axis=1)
    absent = ~pres[grp].any(axis=1) & pres[rest].all(axis=1)
    for p in L.index[only]:
        onoff.append((name, "present_only_in_group", p, gene[p], float(L.loc[p, grp].mean())))
    for p in L.index[absent]:
        onoff.append((name, "absent_only_in_group", p, gene[p], float(L.loc[p, rest].mean())))
onoff = pd.DataFrame(onoff, columns=["group", "pattern", "protein", "gene", "mean_log2_where_present"])
onoff.to_csv(TAB / "presence_absence_proteins.csv", index=False)
print(onoff.groupby(["group", "pattern"]).gene.apply(lambda s: ", ".join(s[:40])).to_string())

# ---- concordance: in-vitro EBV infection (eAGS vs AGS) vs natural EBV+ lines -----
a = results["eAGS_vs_AGS_replicates"]
b = results["EBVnatural_vs_EBVneg"]
j = a[["logFC", "q", "gene"]].join(b[["logFC", "q"]], lsuffix="_eAGS", rsuffix="_EBVnat", how="inner")
from scipy.stats import fisher_exact, spearmanr
rho, prho = spearmanr(j.logFC_eAGS, j.logFC_EBVnat)
sigA = j.q_eAGS < 0.05
sigB = j.q_EBVnat < 0.05
both_up = (sigA & sigB & (j.logFC_eAGS > 0) & (j.logFC_EBVnat > 0))
both_dn = (sigA & sigB & (j.logFC_eAGS < 0) & (j.logFC_EBVnat < 0))
opp = sigA & sigB & (np.sign(j.logFC_eAGS) != np.sign(j.logFC_EBVnat))
summary["eAGS_vs_natural_EBV_concordance"] = {
    "n_common": int(len(j)), "spearman_rho_logFC": float(rho), "p": float(prho),
    "n_sig_both_up": int(both_up.sum()), "n_sig_both_down": int(both_dn.sum()), "n_sig_opposite": int(opp.sum()),
    "genes_both_up": j[both_up].sort_values("logFC_EBVnat", ascending=False).gene.tolist()[:40],
    "genes_both_down": j[both_dn].sort_values("logFC_EBVnat").gene.tolist()[:40],
}
j.to_csv(TAB / "eAGS_vs_naturalEBV_logFC.csv.gz", float_format="%.4g")
print(json.dumps(summary["eAGS_vs_natural_EBV_concordance"], indent=1))

json.dump(summary, open(TAB / "differential_summary.json", "w"), indent=2, default=float)


# ---- figures ------------------------------------------------------------------
def repel_labels(ax, items, fontsize=6.5):
    """Greedy vertical de-collision of (x, y, text) labels in data coords."""
    import numpy as np
    fig = ax.figure
    fig.canvas.draw()
    items = sorted(items, key=lambda t: (-t[1], t[0]))
    placed = []  # (x_disp, y_disp) of placed label anchors
    inv = ax.transData.inverted()
    for x, y, txt in items:
        xd, yd = ax.transData.transform((x, y))
        step = 12.0
        for k in range(40):
            cand = yd + (step * ((k + 1) // 2) * (1 if k % 2 else -1))
            if all(abs(cand - py) > 11 or abs(xd - px) > 46 for px, py in placed):
                yd = cand
                break
        placed.append((xd, yd))
        xa, ya = inv.transform((xd, yd))
        ax.annotate(txt, (x, y), xytext=(xa, ya), textcoords="data", fontsize=fontsize,
                    color=INK2, ha="left",
                    arrowprops=dict(arrowstyle="-", color="#c3c2b7", lw=0.4,
                                    shrinkA=0, shrinkB=2))


def volcano(ax, r, title, label_n=12, extra=()):
    sig = r.q < 0.05
    ax.scatter(r.logFC[~sig], -np.log10(r.p[~sig]), s=4, color="#c3c2b7", linewidths=0)
    up = sig & (r.logFC > 0)
    dn = sig & (r.logFC < 0)
    ax.scatter(r.logFC[up], -np.log10(r.p[up]), s=6, color=CAT[1], linewidths=0)
    ax.scatter(r.logFC[dn], -np.log10(r.p[dn]), s=6, color=CAT[0], linewidths=0)
    lab = r[sig].assign(score=lambda d: -np.log10(d.p) * d.logFC.abs()).sort_values("score",
                                                                                    ascending=False)
    lab = pd.concat([lab.head(label_n), r[r.gene.isin(extra)]]).drop_duplicates()
    repel_labels(ax, [(row.logFC, -np.log10(row.p), row.gene.split(";")[0]) for _, row in lab.iterrows()])
    ax.axvline(0, color="#c3c2b7", lw=0.6)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("−log10 p (moderated t)")
    ax.set_title(title + f"\n{int(up.sum())} up / {int(dn.sum())} down at q<0.05", fontsize=9)


fig, axes = plt.subplots(2, 2, figsize=(11, 9), gridspec_kw={"hspace": 0.45, "wspace": 0.28})
volcano(axes[0, 0], results["paperEMT_vs_paperMetabolism"], "A  Hu et al. EMT vs Metabolism (7 vs 7 lines)",
        extra=("VIM", "CDH1", "GLUL", "ALDH3A1"))
volcano(axes[0, 1], results["EBVnatural_vs_EBVneg"], "B  Natural EBV+ (3) vs EBV− (10) lines",
        extra=("CDKN2A", "CD274"))
volcano(axes[1, 0], results["clusterEBVlike_vs_rest"], "C  Cluster {TMK1, NCC24, SNU719, YCCEL1} vs rest")
volcano(axes[1, 1], results["eAGS_vs_AGS_replicates"], "D  eAGS vs AGS (isogenic, 3 vs 3 replicates)",
        extra=("CDKN2A",))
savefig(fig, "09_volcanoes")
plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 4.6))
ax.scatter(j.logFC_EBVnat, j.logFC_eAGS, s=3, color="#c3c2b7", linewidths=0)
ax.scatter(j.logFC_EBVnat[both_up | both_dn], j.logFC_eAGS[both_up | both_dn], s=10, color=CAT[1],
           linewidths=0, label="q<0.05 in both, same direction")
ax.scatter(j.logFC_EBVnat[opp], j.logFC_eAGS[opp], s=10, color=CAT[0], linewidths=0,
           label="q<0.05 in both, opposite")
jl = j[both_up | both_dn].assign(s=lambda d: d.logFC_EBVnat.abs() + d.logFC_eAGS.abs()) \
    .sort_values("s", ascending=False).head(14)
repel_labels(ax, [(row.logFC_EBVnat, row.logFC_eAGS, row.gene.split(";")[0]) for _, row in jl.iterrows()])
ax.axhline(0, color="#c3c2b7", lw=0.6); ax.axvline(0, color="#c3c2b7", lw=0.6)
ax.set_xlabel("log2FC natural EBV+ vs EBV− lines")
ax.set_ylabel("log2FC eAGS vs AGS")
ax.set_title(f"In-vitro vs natural EBV (Spearman ρ = {rho:.2f})")
ax.legend(fontsize=7, loc="lower right")
savefig(fig, "10_eAGS_vs_naturalEBV")
plt.close(fig)
