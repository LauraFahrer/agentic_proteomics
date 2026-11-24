#!/usr/bin/env python3
"""Step 04 - final subtypes, grid ranking of the paper's parameters, integrative consensus clustering
(self-implemented Monti/ConsensusClusterPlus procedure), Kim-Park metagene features (Fig. 2A, 4A/C, 5A, S2A/B, S3A, S4)."""
import glob
import json
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from common import *
import methods as M

ND = os.path.join(DC, "nmf")
grid = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(TAB, "03_nmf_grid_*_mean.csv")))])
grid.to_csv(os.path.join(TAB, "04_nmf_grid_all.csv"), index=False)

# ---- Fig. S2A / S4A / S4E
fig, axs = plt.subplots(2, 2, figsize=(10, 9))
for ax, d in zip(axs.flat, ["WP", "PP", "YP", "WPq"]):
    g = grid[grid.dataset == d]
    for fr, s in g.groupby("frac"):
        ax.plot(s.k, s.cophenetic, "o-", label=f"{int(fr*100)}% (n={s.n_features.iloc[0]})")
    ax.set_title(d); ax.set_xlabel("k"); ax.set_ylabel("cophenetic"); ax.set_xticks([2, 3, 4, 5]); ax.legend(fontsize=7)
savefig(fig, "04_FigS2A_S4AE_cophenetic_grid")

paper_set = pd.DataFrame(dict(dataset=["WP", "PP", "YP", "WPq"], frac=[0.1, 0.2, 1.0, 0.1], k=[3, 2, 3, 3],
                              paper_n_features=[568, 1691, 4780, np.nan],
                              paper_mean_silhouette=[0.87, 1.0, 0.92, np.nan]))
sel = paper_set.merge(grid, on=["dataset", "frac", "k"])
sel["rank_of_cophenetic_in_grid"] = [int((grid[grid.dataset == r.dataset].cophenetic > r.cophenetic).sum() + 1) for r in sel.itertuples()]
sel.to_csv(os.path.join(TAB, "04_paper_parameters_in_grid.csv"), index=False)
best = grid.loc[grid.groupby("dataset").cophenetic.idxmax()]
best.to_csv(os.path.join(TAB, "04_grid_best_cophenetic.csv"), index=False)
print(sel.to_string()); print(best.to_string())

assign = []
for r in sel.itertuples():
    z = np.load(os.path.join(ND, f"{r.dataset}_mean_{int(r.frac*100):03d}_k{r.k}.npz"), allow_pickle=True)
    smp = list(z["samples"]); pl = paper_label(r.dataset, smp)
    assign.append(pd.DataFrame(dict(dataset=r.dataset, sample=smp, line=[line_of(s) for s in smp],
                                    cluster=z["cl"], silhouette=z["sil"], paper_cluster=pl)))
    C = z["C"]; Z, _ = M.consensus_hclust(C, r.k); o = hierarchy.leaves_list(Z)
    fig, ax = plt.subplots(figsize=(8, 7.5))
    im = ax.imshow(C[np.ix_(o, o)], cmap="RdYlBu_r", vmin=0, vmax=1)
    lab = [f"{smp[i]} [p{pl[i]}|c{z['cl'][i]}]" for i in o]
    ax.set_yticks(range(len(o))); ax.set_yticklabels(lab, fontsize=5); ax.set_xticks([])
    ax.set_title(f"{r.dataset}: top {int(r.frac*100)}% (n={r.n_features}), k={r.k}, coph={r.cophenetic:.3f}\n"
                 f"[p = paper cluster | c = this run]", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.6, label="consensus")
    savefig(fig, f"04_consensus_{r.dataset}_{int(r.frac*100)}pct_k{r.k}")
A = pd.concat(assign); A.to_csv(os.path.join(TAB, "04_subtype_assignments.csv"), index=False)
comp = []
for d, g in A.groupby("dataset"):
    comp.append(dict(dataset=d, ARI_vs_paper=ari(g.cluster, g.paper_cluster), mean_silhouette=g.silhouette.mean(),
                     contingency=pd.crosstab(g.cluster, g.paper_cluster).to_dict()))
comp = pd.DataFrame(comp); comp.to_csv(os.path.join(TAB, "04_subtypes_vs_paper.csv"), index=False); print(comp.to_string())
wq = A[A.dataset == "WPq"].groupby("line").cluster.apply(lambda x: "".join(map(str, x)))
print("WPq per line:", wq.to_dict())

# ---- integrative consensus clustering (Methods 2.10; Fig. 5A)
yp_s = list(A[A.dataset == "YP"]["sample"])
rows, names = [], []
for d in ["WP", "PP", "YP"]:
    g = A[A.dataset == d].set_index("sample").loc[yp_s]
    for c in sorted(g.cluster.unique()):
        rows.append((g.cluster.values == c).astype(float)); names.append(f"{d}subtype{c}")
B = np.array(rows)
pd.DataFrame(B, index=names, columns=yp_s).to_csv(os.path.join(TAB, "04_integrative_binary_table.csv"))
cc = M.consensus_clustering(B, kmax=4, reps=500, p_item=0.8, seed=123456)
kst = []
for k, v in cc.items():
    C, cl = v["consensus"], v["classes"]
    kst.append(dict(k=k, sizes="/".join(map(str, np.bincount(cl)[1:])),
                    mean_within_consensus=np.mean([C[np.ix_(cl == c, cl == c)].mean() for c in np.unique(cl)])))
pd.DataFrame(kst).to_csv(os.path.join(TAB, "04_integrative_k_stats.csv"), index=False)
ic = cc[2]["classes"]
met = pd.Series([PAPER_INT[line_of(s)] == "Metabolism" for s in yp_s]).groupby(ic).mean().idxmax()
IS = pd.DataFrame(dict(sample=yp_s, line=[line_of(s) for s in yp_s], integrative_cluster=ic,
                       integrative_subtype=np.where(ic == met, "Metabolism", "EMT"),
                       paper_integrative=[PAPER_INT[line_of(s)] for s in yp_s]))
IS.to_csv(os.path.join(TAB, "04_integrative_assignments_samples.csv"), index=False)
L = IS.groupby("line").agg(integrative_subtype=("integrative_subtype", lambda x: x.iloc[0] if x.nunique() == 1 else "discordant"),
                           paper_integrative=("paper_integrative", "first")).reset_index()
L.to_csv(os.path.join(TAB, "03_integrative_clusters.csv"), index=False)
print(L.to_string()); print("integrative ARI vs paper:", ari(IS.integrative_subtype, IS.paper_integrative))
o = np.lexsort((IS.line, IS.integrative_subtype))
fig, ax = plt.subplots(figsize=(11, 3.5))
ax.imshow(B[:, o], cmap="gray_r", aspect="auto"); ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=7)
ax.set_xticks(range(len(o))); ax.set_xticklabels([f"{yp_s[i]} ({IS.integrative_subtype[i][:3]})" for i in o], rotation=90, fontsize=6)
ax.set_title("Integrative consensus clustering (500 resamplings, k = 2); black = member")
savefig(fig, "04_Fig5A_integrative_binary_heatmap")

# ---- Kim-Park metagene features (Fig. S3A)
paper42 = ["UBE2E2", "CD55", "LGALS1", "MOV10", "SLC2A3", "PTK7", "CRABP2", "ITGA2", "ANXA3", "EPS8L2", "SH3KBP1", "TACSTD2",
           "APOBEC3B", "ISG15", "BCAM", "AKR1C3", "STARD10", "AGR2", "AGR3", "ALDH3A1", "ALDH3A2", "VIL1", "S100A6", "ITGB4",
           "CTSE", "LGALS3BP", "LYZ", "LGALS4", "FXYD3", "ANXA5", "ALDH1A1", "LZTS1", "ANKRD45", "PHGDH", "ANXA6", "PFN2",
           "MAGEA4", "TGM2", "VIM", "FSCN1", "TPD52L1", "FKBP10"]
z = np.load(os.path.join(ND, "WP_mean_010_k3.npz"), allow_pickle=True)
score, sel_f = M.kim_park_features(z["W"])
feats = z["features"]
mg = pd.DataFrame([(feats[i], q + 1) for q, ii in enumerate(sel_f) for i in ii], columns=["protein", "metagene"])
mg["in_paper_42"] = mg.protein.isin(paper42); mg.to_csv(os.path.join(TAB, "04_WP_metagene_features.csv"), index=False)
json.dump(dict(n_features=len(mg), overlap_paper42=int(mg.in_paper_42.sum()), paper42_in_top10pct=int(np.isin(paper42, feats).sum())),
          open(os.path.join(TAB, "04_WP_metagene_summary.json"), "w"), indent=1)
print("metagenes:", len(mg), "overlap with paper 42:", mg.in_paper_42.sum())
