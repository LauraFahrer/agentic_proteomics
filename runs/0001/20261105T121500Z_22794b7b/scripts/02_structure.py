"""Global structure of the proteome: pipeline concordance with the original study,
PCA, consensus clustering and agreement with the published EMT / Metabolism subtypes,
and per-line pathway (hallmark) scores.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, spearmanr
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import NMF, PCA
from sklearn.metrics import adjusted_rand_score

from common import (CELL_META, DATA, FIG, INK2, SEED, SUBTYPE_COLORS, TAB, diverging_cmap,
                    read_gmt, savefig, style)

style()
rng = np.random.default_rng(SEED)
res = {}

ann = pd.read_csv(TAB / "sample_annotation.csv")
X = pd.read_csv(TAB / "protein_matrix_log2.csv.gz", index_col=0)
gene = X.pop("gene").fillna("")
L = pd.read_csv(TAB / "protein_matrix_cellline_mean_log2.csv.gz", index_col=0)
L.pop("gene")
lines = list(L.columns)
sub = CELL_META.paper_subtype.reindex(lines)

# --------------------------------------------------------------------------
# 1. Concordance with the authors' Spectronaut-processed WP matrix (Table S3)
# --------------------------------------------------------------------------
wp = pd.read_excel(DATA / "supp/x/MOL2-18-1437-s005.xls", sheet_name="processed WP matrix")
wp = wp.set_index("Protein Accessions")
wpg = wp.pop("Genes")
rename = {"Ocum1": "Ocum1"}
wp_line = pd.DataFrame({ln: wp[[c for c in wp.columns if c.rsplit("_", 1)[0] == ln]].mean(axis=1)
                        for ln in lines})
# map quantms protein groups to Spectronaut groups on the exact accession string first,
# then on leading accession
q_lead = pd.Series(L.index.str.split(";").str[0], index=L.index)
s_lead = pd.Series(wp_line.index.str.split(";").str[0], index=wp_line.index)
s_map = dict(zip(s_lead.values, s_lead.index))
pairs = [(q, s_map[a]) for q, a in q_lead.items() if a in s_map]
res["spectronaut_proteins"] = int(len(wp_line))
res["quantms_proteins"] = int(len(L))
res["matched_by_leading_accession"] = len(pairs)
rs = []
for q, s in pairs:
    a, b = L.loc[q], wp_line.loc[s]
    ok = a.notna() & b.notna()
    if ok.sum() >= 8 and a[ok].std() > 0 and b[ok].std() > 0:
        rs.append((q, s, pearsonr(a[ok], b[ok])[0], a[ok].std()))
rs = pd.DataFrame(rs, columns=["quantms_pg", "spectronaut_pg", "r_across_lines", "sd_quantms"])
rs["gene"] = gene.reindex(rs.quantms_pg).values
rs.to_csv(TAB / "concordance_quantms_vs_spectronaut.csv", index=False)
res["concordance_median_r"] = float(rs.r_across_lines.median())
hv = rs.sd_quantms > rs.sd_quantms.quantile(0.75)
res["concordance_median_r_top25pct_variable"] = float(rs[hv].r_across_lines.median())
res["concordance_frac_r_gt_0.8_top25pct_variable"] = float((rs[hv].r_across_lines > 0.8).mean())

fig, ax = plt.subplots(1, 2, figsize=(9, 3.2), gridspec_kw={"wspace": 0.35})
ax[0].hist(rs.r_across_lines, bins=50, color="#2a78d6")
ax[0].set_xlabel("Pearson r across 14 lines (quantms vs. Spectronaut)")
ax[0].set_ylabel("proteins")
ax[0].set_title("A  Per-protein concordance with the published matrix")
ax[1].scatter(rs.sd_quantms, rs.r_across_lines, s=3, alpha=0.3, color="#2a78d6", linewidths=0)
ax[1].set_xlabel("between-line SD (quantms, log2)")
ax[1].set_ylabel("Pearson r")
ax[1].set_title("B  Concordance rises with biological variance")
savefig(fig, "05_concordance_with_published_matrix")
plt.close(fig)

# --------------------------------------------------------------------------
# 2. PCA on proteins quantified in every run
# --------------------------------------------------------------------------
Xc = X.dropna()
res["pca_complete_proteins"] = int(len(Xc))
Z = Xc.sub(Xc.mean(axis=1), axis=0).T
pca = PCA(n_components=6, random_state=SEED).fit(Z)
pcs = pd.DataFrame(pca.transform(Z)[:, :4], index=Z.index, columns=["PC1", "PC2", "PC3", "PC4"])
ev = pca.explained_variance_ratio_
res["pca_explained_variance"] = [round(float(v), 4) for v in ev]
load = pd.DataFrame(pca.components_[:3].T, index=Xc.index, columns=["PC1", "PC2", "PC3"])
load["gene"] = gene.reindex(load.index)
load.to_csv(TAB / "pca_loadings.csv")
fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), gridspec_kw={"wspace": 0.3})
for ax, (a, b) in zip(axes, [("PC1", "PC2"), ("PC3", "PC4")]):
    for sname, c in SUBTYPE_COLORS.items():
        m = ann.paper_subtype.values == sname
        ax.scatter(pcs[a][m], pcs[b][m], s=34, color=c, label=sname, edgecolor="white", linewidth=0.8)
    for ln in lines:
        m = (ann.line == ln).values
        ax.annotate(ln, (pcs[a][m].mean(), pcs[b][m].mean()), fontsize=7, color=INK2,
                    xytext=(4, 4), textcoords="offset points")
    i, j = int(a[2]) - 1, int(b[2]) - 1
    ax.set_xlabel(f"{a} ({ev[i]*100:.1f}%)")
    ax.set_ylabel(f"{b} ({ev[j]*100:.1f}%)")
axes[0].legend(title="Hu et al. subtype", loc="best", fontsize=8)
axes[0].set_title("A  PCA, %d proteins quantified in all 42 runs" % len(Xc))
axes[1].set_title("B  Higher components")
savefig(fig, "06_pca")
plt.close(fig)

# --------------------------------------------------------------------------
# 3. Consensus clustering of cell lines (line means, proteins present in all lines)
# --------------------------------------------------------------------------
Lc = L.dropna()
mad = (Lc.sub(Lc.median(axis=1), axis=0)).abs().median(axis=1)


def cluster_labels(M, k, method="average"):
    Zl = linkage(pdist(M.T.values, "correlation"), method)
    from scipy.cluster.hierarchy import fcluster
    return fcluster(Zl, k, "maxclust")


truth = (sub == "EMT").astype(int).values
cons_res = []
cons_mats = {}
for top in (0.1, 0.2, 0.3, 0.5):
    feats = mad.sort_values(ascending=False).index[: int(len(mad) * top)]
    for k in (2, 3, 4):
        C = np.zeros((len(lines), len(lines)))
        aris = []
        B = 200
        for b in range(B):
            f = rng.choice(feats, size=int(0.8 * len(feats)), replace=False)
            M = Lc.loc[f]
            M = M.sub(M.mean(axis=1), axis=0)
            lab = cluster_labels(M, k)
            C += lab[:, None] == lab[None, :]
            aris.append(adjusted_rand_score(truth, lab))
        C /= B
        cons_mats[(top, k)] = C
        # consensus label: hierarchical on 1 - C
        lab = AgglomerativeClustering(n_clusters=k, metric="precomputed", linkage="average") \
            .fit_predict(1 - C)
        cons_res.append({"top_mad_fraction": top, "n_features": len(feats), "k": k,
                         "ARI_consensus_vs_paper": adjusted_rand_score(truth, lab),
                         "ARI_bootstrap_median": float(np.median(aris)),
                         "clusters": "; ".join(",".join(np.array(lines)[lab == c])
                                               for c in np.unique(lab))})
cons_res = pd.DataFrame(cons_res)
cons_res.to_csv(TAB / "consensus_clustering_vs_paper_subtypes.csv", index=False)
print(cons_res.to_string())

# NMF (as in the paper; data shifted to be non-negative), k = 2/3 on top 10% MAD
nmf_res = []
feats10 = mad.sort_values(ascending=False).index[: int(len(mad) * 0.1)]
M = Lc.loc[feats10]
M = M.sub(M.mean(axis=1), axis=0)
Mpos = np.concatenate([np.clip(M.values, 0, None), np.clip(-M.values, 0, None)])
for k in (2, 3):
    labs = []
    for s in range(50):
        W = NMF(n_components=k, init="random", random_state=s, max_iter=2000).fit(Mpos.T)
        labs.append(W.transform(Mpos.T).argmax(axis=1))
    C = np.mean([(l[:, None] == l[None, :]) for l in labs], axis=0)
    lab = AgglomerativeClustering(n_clusters=k, metric="precomputed", linkage="average") \
        .fit_predict(1 - C)
    nmf_res.append({"k": k, "ARI_vs_paper": adjusted_rand_score(truth, lab),
                    "clusters": "; ".join(",".join(np.array(lines)[lab == c]) for c in np.unique(lab))})
pd.DataFrame(nmf_res).to_csv(TAB / "nmf_clustering_vs_paper_subtypes.csv", index=False)
print(pd.DataFrame(nmf_res).to_string())
res["nmf"] = nmf_res

# figure: consensus heatmap (top 10%, k=2) + dendrogram of line correlation
C = cons_mats[(0.1, 2)]
Zl = linkage(pdist(M.T.values, "correlation"), "average")
order = leaves_list(Zl)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"wspace": 0.45, "width_ratios": [1, 1.1]})
ax = axes[0]
d = dendrogram(Zl, labels=lines, ax=ax, color_threshold=0, above_threshold_color=INK2,
               leaf_rotation=90)
for t in ax.get_xticklabels():
    t.set_color(SUBTYPE_COLORS[sub[t.get_text()]])
ax.set_ylabel("1 − Pearson r")
ax.set_title("A  Line clustering (top 10%% MAD, %d proteins)" % len(feats10))
ax = axes[1]
im = ax.imshow(C[np.ix_(order, order)], cmap="Blues", vmin=0, vmax=1)
ax.set_xticks(range(len(lines)), np.array(lines)[order], rotation=90)
ax.set_yticks(range(len(lines)), np.array(lines)[order])
for t in ax.get_xticklabels() + ax.get_yticklabels():
    t.set_color(SUBTYPE_COLORS[sub[t.get_text()]])
fig.colorbar(im, ax=ax, shrink=0.8, label="co-clustering frequency (k=2, 200 resamples)")
ax.set_title("B  Consensus matrix")
fig.text(0.01, -0.12, "Label colour = Hu et al. integrative subtype (orange EMT, blue Metabolism)",
         fontsize=8, color=INK2)
savefig(fig, "07_clustering_consensus")
plt.close(fig)

# --------------------------------------------------------------------------
# 4. Hallmark pathway scores per cell line (mean z-score of member proteins)
# --------------------------------------------------------------------------
sym = gene.reindex(Lc.index).str.split(";").str[0]
Zp = Lc.sub(Lc.mean(axis=1), axis=0).div(Lc.std(axis=1), axis=0)
Zp.index = sym
Zp = Zp[~Zp.index.duplicated()]
hall = read_gmt(DATA / "genesets/h.all.v2023.1.Hs.symbols.gmt", set(Zp.index), min_size=15)
score = pd.DataFrame({k.replace("HALLMARK_", ""): Zp.loc[list(v)].mean() for k, v in hall.items()}).T
nset = pd.Series({k.replace("HALLMARK_", ""): len(v) for k, v in hall.items()})
# association of each hallmark score with paper subtype (Welch t, 7 vs 7 lines)
from scipy.stats import ttest_ind
tt = []
for h in score.index:
    a, b = score.loc[h, sub == "EMT"], score.loc[h, sub == "Metabolism"]
    t, p = ttest_ind(a, b, equal_var=False)
    tt.append((h, nset[h], a.mean() - b.mean(), t, p))
tt = pd.DataFrame(tt, columns=["hallmark", "n_proteins", "diff_EMT_minus_Metab", "t", "p"])
from common import bh
tt["q"] = bh(tt.p)
tt = tt.sort_values("p")
tt.to_csv(TAB / "hallmark_scores_vs_paper_subtype.csv", index=False)
score.to_csv(TAB / "hallmark_scores_per_line.csv")
print(tt.head(15).to_string())

show = tt.head(20).hallmark.tolist()
S = score.loc[show]
col_order = list(sub.sort_values().index)
fig, ax = plt.subplots(figsize=(8.5, 6.5))
v = np.nanmax(np.abs(S.values))
im = ax.imshow(S[col_order].values, cmap=diverging_cmap(), vmin=-v, vmax=v, aspect="auto")
ax.set_xticks(range(len(col_order)), col_order, rotation=90)
for t in ax.get_xticklabels():
    t.set_color(SUBTYPE_COLORS[sub[t.get_text()]])
ax.set_yticks(range(len(show)), [f"{h}  (q={q:.2g})" for h, q in zip(show, tt.head(20).q)], fontsize=7)
fig.colorbar(im, ax=ax, shrink=0.6, label="mean protein z-score")
ax.set_title("Hallmark protein scores, 20 most subtype-associated (Welch t, 7 vs 7 lines)")
savefig(fig, "08_hallmark_scores")
plt.close(fig)

res["hallmark_top"] = tt.head(10)[["hallmark", "diff_EMT_minus_Metab", "q"]].round(4).values.tolist()
json.dump(res, open(TAB / "structure_summary.json", "w"), indent=2, default=float)
print(json.dumps({k: v for k, v in res.items() if k not in ("hallmark_top",)}, indent=2, default=float))
