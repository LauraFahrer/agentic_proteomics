"""Co-abundance modules (candidate complexes / functional units), their annotation,
subtype association, and identification of strongly co-regulated pairs not in known complexes.
"""
import json
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom, pointbiserialr, ttest_ind

from common import CELL_META, DATA, INK2, SUBTYPE_COLORS, TAB, diverging_cmap, read_gmt, savefig, style

style()
out = {}
G = pd.read_csv(TAB / "gene_matrix_cellline_mean_log2.csv.gz", index_col=0)
C = pd.read_parquet(TAB / "coabundance_bicor.parquet")
lines = list(G.columns)
sub = CELL_META.paper_subtype.reindex(lines)

# most variable genes carry the biology
mad = (G.sub(G.median(axis=1), axis=0)).abs().median(axis=1)
topg = mad.sort_values(ascending=False).head(2000).index
Cm = C.loc[topg, topg].copy()
np.fill_diagonal(Cm.values, 1.0)
Cm = Cm.fillna(0)
D = 1 - Cm.values
np.fill_diagonal(D, 0)
D = (D + D.T) / 2
Z = linkage(squareform(D, checks=False), method="average")
labels = fcluster(Z, t=0.75, criterion="distance")   # merge genes with bicor>~0.25 avg
modules = {}
for m in np.unique(labels):
    members = list(topg[labels == m])
    if len(members) >= 8:
        modules[int(m)] = members
out["n_modules"] = len(modules)
out["genes_in_modules"] = int(sum(len(v) for v in modules.values()))

# gene-set libraries for annotation
universe = set(G.index)
libs = {}
for f, tag in [("c5.go.bp.v7.5.1.symbols.gmt", "GOBP"), ("c5.go.cc.v7.5.1.symbols.gmt", "GOCC"),
               ("c2.cp.reactome.v7.5.1.symbols.gmt", "REACTOME"), ("c2.cp.kegg.v7.5.1.symbols.gmt", "KEGG")]:
    for k, v in read_gmt(DATA / "genesets" / f, universe, 3, 400).items():
        libs[f"{tag}:{k}"] = v
complex_sets = {k: v for k, v in libs.items() if k.startswith("GOCC:") and "COMPLEX" in k}
Nuni = len(universe)


def enrich(members, sets, topn=3):
    m = set(members) & universe
    res = []
    for name, gs in sets.items():
        k = len(m & gs)
        if k < 3:
            continue
        p = hypergeom.sf(k - 1, Nuni, len(gs), len(m))
        res.append((name, k, len(gs), p))
    res.sort(key=lambda x: x[3])
    return res[:topn]


# module eigengene = PC1; subtype association
rows = []
mod_eigen = {}
for m, members in modules.items():
    X = G.loc[members]
    Xz = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1) + 1e-9, axis=0)
    u, s, vt = np.linalg.svd(Xz.values - Xz.values.mean(0), full_matrices=False)
    eig = pd.Series(vt[0], index=lines)
    if np.corrcoef(eig, Xz.mean(0))[0, 1] < 0:
        eig = -eig
    mod_eigen[m] = eig
    a, b = eig[sub == "EMT"], eig[sub == "Metabolism"]
    t, p = ttest_ind(a, b, equal_var=False)
    top_go = enrich(members, libs, 3)
    top_cx = enrich(members, complex_sets, 1)
    # within-module mean co-abundance (cohesion)
    sub_c = C.loc[members, members].values
    cohesion = np.nanmean(sub_c[np.triu_indices(len(members), 1)])
    rows.append({"module": m, "size": len(members), "cohesion_bicor": float(cohesion),
                 "subtype_diff_EMT_minus_Metab": float(a.mean() - b.mean()), "subtype_p": float(p),
                 "top_annotation": top_go[0][0] if top_go else "", "top_annotation_p": top_go[0][3] if top_go else np.nan,
                 "annotations": "; ".join(f"{n}(p={p_:.1g})" for n, k, s_, p_ in top_go),
                 "known_complex": top_cx[0][0].replace("GOCC:", "") if top_cx else "",
                 "members": ",".join(members)})
modtab = pd.DataFrame(rows).sort_values("cohesion_bicor", ascending=False)
from common import bh
modtab["subtype_q"] = bh(modtab.subtype_p.values)
modtab.to_csv(TAB / "coabundance_modules.csv", index=False)
out["modules_subtype_assoc_q<0.1"] = int((modtab.subtype_q < 0.1).sum())
print(modtab[["module", "size", "cohesion_bicor", "top_annotation", "known_complex",
              "subtype_diff_EMT_minus_Metab", "subtype_q"]].head(25).to_string())

# --- strongly co-regulated pairs NOT in any known GO-CC complex (candidate novel) ---
known_pairs = set()
for gs in complex_sets.values():
    for a, b in itertools.combinations(sorted(gs), 2):
        known_pairs.add((a, b))
iu = np.triu_indices(len(G.index), 1)
gi = np.array(G.index)
strong = []
Cv = C.values
for i, j in zip(*np.triu_indices(len(gi), 1)):
    pass  # (vectorised below)
mask = np.triu(np.ones(Cv.shape, bool), 1) & (Cv >= 0.9)
ii, jj = np.where(mask)
for i, j in zip(ii, jj):
    a, b = gi[i], gi[j]
    key = (a, b) if a < b else (b, a)
    strong.append((a, b, float(Cv[i, j]), key in known_pairs))
strong = pd.DataFrame(strong, columns=["geneA", "geneB", "bicor", "in_known_complex"])
strong.to_csv(TAB / "strong_coabundance_pairs.csv.gz", index=False)
out["strong_pairs_bicor>=0.9"] = int(len(strong))
out["strong_pairs_novel_frac"] = float((~strong.in_known_complex).mean())
json.dump(out, open(TAB / "modules_summary.json", "w"), indent=2, default=float)
print(json.dumps(out, indent=2, default=float))

# --- figure: module eigengene x line heatmap, modules annotated & subtype-sorted ---
show = modtab[modtab["size"] >= 10].copy()
E = pd.DataFrame({f"M{r.module} ({r['size']})": mod_eigen[r.module] for _, r in show.iterrows()}).T
col_order = list(sub.sort_values().index)
fig, ax = plt.subplots(figsize=(8.5, max(4, 0.32 * len(E))))
im = ax.imshow(E[col_order].values, cmap=diverging_cmap(), vmin=-0.6, vmax=0.6, aspect="auto")
ax.set_xticks(range(len(col_order)), col_order, rotation=90)
for t in ax.get_xticklabels():
    t.set_color(SUBTYPE_COLORS[sub[t.get_text()]])
lab = [f"M{r.module}: {(r.known_complex or r.top_annotation).replace('GOBP_','').replace('REACTOME_','').replace('HALLMARK_','')[:34]}"
       + (" *" if r.subtype_q < 0.1 else "") for _, r in show.iterrows()]
ax.set_yticks(range(len(E)), lab, fontsize=6.5)
fig.colorbar(im, ax=ax, shrink=0.5, label="module eigengene")
ax.set_title("Co-abundance modules across lines (* = subtype-associated q<0.1)\n"
             "label colour of lines: orange EMT, blue Metabolism", fontsize=9)
savefig(fig, "02_module_landscape")
plt.close(fig)
