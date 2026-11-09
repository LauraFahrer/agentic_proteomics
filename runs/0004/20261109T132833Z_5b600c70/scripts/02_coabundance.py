"""Protein co-abundance network and benchmark against known complexes (prompt_0004, PPI).

Rationale: obligate complex subunits are co-regulated, so their abundances co-vary across
conditions. We compute protein-protein co-abundance across the 14 cell lines and test whether
it recovers known complex membership (GO cellular-component 'complex' terms) better than chance.
This validates the co-abundance signal before any novel-interaction claim.
"""
import json
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CAT, DATA, INK2, TAB, read_gmt, savefig, style

style()
rng = np.random.default_rng(20261109)
out = {}
L = pd.read_csv(TAB / "protein_matrix_cellline_mean_log2.csv.gz", index_col=0)
gene = L.pop("gene").fillna("")
sym = gene.str.split(";").str[0]
# proteins quantified in all 14 lines and with a gene symbol (for complex mapping)
Lc = L[(L.notna().sum(axis=1) == L.shape[1]) & (sym.reindex(L.index).values != "")]
symc = sym.reindex(Lc.index)
# collapse to gene level (mean over protein groups of the same leading gene)
G = Lc.groupby(symc.values).mean()
out["genes_complete"] = int(G.shape[0])


def biweight_midcorr(X):
    """Rows = genes, columns = samples. Returns gene x gene bicor matrix (robust to outliers)."""
    X = X.values.astype(float)
    med = np.median(X, axis=1, keepdims=True)
    mad = np.median(np.abs(X - med), axis=1, keepdims=True)
    mad[mad == 0] = 1e-9
    u = (X - med) / (9 * mad)
    w = (1 - u ** 2) ** 2
    w[np.abs(u) >= 1] = 0
    num = (X - med) * w
    denom = np.sqrt((num ** 2).sum(axis=1, keepdims=True))
    denom[denom == 0] = 1e-9
    nx = num / denom
    return nx @ nx.T


C = biweight_midcorr(G)
Cdf = pd.DataFrame(C, index=G.index, columns=G.index)
np.fill_diagonal(Cdf.values, np.nan)

# --- benchmark: within-complex vs random pairs (GO-CC complex terms) ---
gocc = read_gmt(DATA / "genesets/c5.go.cc.v7.5.1.symbols.gmt", set(G.index), min_size=3, max_size=100)
complexes = {k: v for k, v in gocc.items() if "COMPLEX" in k}
out["gocc_complex_terms_used"] = len(complexes)
pos_pairs = set()
for members in complexes.values():
    ms = sorted(members)
    for a, b in itertools.combinations(ms, 2):
        pos_pairs.add((a, b) if a < b else (b, a))
genes = list(G.index)
gi = {g: i for i, g in enumerate(genes)}
pos = np.array([Cdf.values[gi[a], gi[b]] for a, b in pos_pairs if a in gi and b in gi])
pos = pos[~np.isnan(pos)]
# random background pairs
N = len(genes)
ridx = rng.integers(0, N, size=(200000, 2))
ridx = ridx[ridx[:, 0] != ridx[:, 1]]
posset = pos_pairs
rand = []
for i, j in ridx[:120000]:
    a, b = genes[i], genes[j]
    key = (a, b) if a < b else (b, a)
    if key not in posset:
        rand.append(Cdf.values[i, j])
rand = np.array(rand); rand = rand[~np.isnan(rand)]
out["within_complex_pairs"] = int(len(pos))
out["within_complex_median_bicor"] = float(np.median(pos))
out["random_pairs_median_bicor"] = float(np.median(rand))
# AUC: probability a within-complex pair outranks a random pair
sample_pos = rng.choice(pos, 40000); sample_rand = rng.choice(rand, 40000)
auc = float(np.mean(sample_pos > sample_rand) + 0.5 * np.mean(sample_pos == sample_rand))
out["coabundance_AUC_complex_vs_random"] = auc
from scipy.stats import mannwhitneyu
out["mannwhitney_p"] = float(mannwhitneyu(pos, rand, alternative="greater").pvalue)
print(json.dumps(out, indent=2))

Cdf.to_parquet(TAB / "coabundance_bicor.parquet")
G.to_csv(TAB / "gene_matrix_cellline_mean_log2.csv.gz", float_format="%.4f")
json.dump(out, open(TAB / "coabundance_benchmark.json", "w"), indent=2)

# --- figure: benchmark + example recovered complexes ---
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3), gridspec_kw={"wspace": 0.33})
ax = axes[0]
bins = np.linspace(-1, 1, 60)
ax.hist(rand, bins=bins, density=True, color="#c3c2b7", alpha=0.8, label="random pairs")
ax.hist(pos, bins=bins, density=True, color=CAT[0], alpha=0.65, label="same GO-CC complex")
ax.axvline(np.median(rand), color="#898781", lw=1, ls="--")
ax.axvline(np.median(pos), color=CAT[0], lw=1, ls="--")
ax.set_xlabel("co-abundance (biweight midcorr, 14 lines)")
ax.set_ylabel("density")
ax.set_title(f"A  Co-abundance recovers complexes\nAUC = {auc:.2f} (n={len(pos):,} vs {len(rand):,} pairs)", fontsize=9)
ax.legend(fontsize=7)

# AUC curve across correlation threshold (precision of complex pairs among top edges)
ax = axes[1]
thr = np.linspace(0.3, 0.95, 40)
prec = []
allpairs_pos = {(a, b) if a < b else (b, a) for a, b in pos_pairs}
iu = np.triu_indices(N, 1)
vals = Cdf.values[iu]
pr_a = np.array(genes)[iu[0]]; pr_b = np.array(genes)[iu[1]]
is_pos = np.array([((a, b) if a < b else (b, a)) in allpairs_pos for a, b in zip(pr_a, pr_b)])
for t in thr:
    m = vals >= t
    prec.append(is_pos[m].mean() if m.sum() > 50 else np.nan)
base = is_pos.mean()
ax.plot(thr, prec, color=CAT[0], lw=1.8)
ax.axhline(base, color="#898781", ls="--", lw=1, label=f"baseline {base:.3f}")
ax.set_xlabel("co-abundance threshold")
ax.set_ylabel("fraction of edges within a known complex")
ax.set_title("B  High co-abundance enriches for\nknown complex pairs", fontsize=9)
ax.legend(fontsize=7)
out["baseline_complex_pair_fraction"] = float(base)
out["precision_at_bicor_0.9"] = float(prec[np.argmin(np.abs(thr - 0.9))])

# example: recovered proteasome 20S/19S co-abundance submatrix
ax = axes[2]
ex = [g for g in G.index if g.startswith("PSM")][:22]
if len(ex) >= 6:
    sub = Cdf.loc[ex, ex].fillna(1)
    im = ax.imshow(sub.values, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(ex)), ex, rotation=90, fontsize=5)
    ax.set_yticks(range(len(ex)), ex, fontsize=5)
    fig.colorbar(im, ax=ax, shrink=0.7, label="bicor")
    ax.set_title("C  Example: proteasome (PSM*)\nsubunit co-abundance", fontsize=9)
savefig(fig, "01_coabundance_benchmark")
plt.close(fig)
json.dump(out, open(TAB / "coabundance_benchmark.json", "w"), indent=2)
print("precision@0.9 =", out["precision_at_bicor_0.9"], "baseline =", base)
