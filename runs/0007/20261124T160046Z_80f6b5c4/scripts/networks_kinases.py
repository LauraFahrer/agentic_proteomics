"""Self-implemented co-expression network and kinase-activity methods (run 0007).

WGCNA-type modules (Zhang & Horvath 2005; Langfelder & Horvath 2008):
  soft threshold (scale-free topology fit), unsigned adjacency |cor|^beta, topological overlap (TOM),
  average-linkage clustering of 1-TOM, module detection by a *simplified* recursive tree cut (own rule,
  documented below; NOT the dynamicTreeCut algorithm), module eigengenes (PC1), merging of modules whose
  eigengenes correlate > 1 - merge_cut.
Kinase activity (PhosMap / Casado 2013 definitions, re-implemented with exact symbol+site matching):
  mean      - mean centred log2 value of a kinase's substrate sites
  KSEA      - PhosMap-style KS enrichment of substrates in the per-sample ranked list; activity = -log10(P)*sign
  MLR       - ridge regression (no intercept) of site values on the binary kinase-substrate matrix, lambda by CV
INKA (Beekhof et al. 2019): kinase-centric = own phosphorylation + activation-loop phosphorylation,
  substrate-centric = PhosphoSitePlus + NetworKIN substrate phosphorylation; INKA = sqrt(KC * SC).
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform


# ------------------------------------------------------------------ WGCNA-type
def scale_free_fit(X, powers=range(1, 21), nbreaks=10):
    """X: samples x features. Returns table of signed R^2 of log10 p(k) ~ log10 k (WGCNA pickSoftThreshold)."""
    R = np.abs(np.corrcoef(X.T)).astype(np.float32)
    np.fill_diagonal(R, 0)
    out = []
    for b in powers:
        k = (R ** b).sum(axis=0)
        h, edges = np.histogram(k, bins=nbreaks)
        mids = np.array([k[(k >= edges[i]) & (k <= edges[i + 1])].mean() if h[i] else np.nan for i in range(nbreaks)])
        ok = (h > 0) & np.isfinite(mids) & (mids > 0)
        x, y = np.log10(mids[ok]), np.log10(h[ok] / h.sum())
        sl, ic, r, _, _ = stats.linregress(x, y)
        out.append((b, -np.sign(sl) * r ** 2, sl, np.mean(k)))
    return pd.DataFrame(out, columns=["power", "signed_R2", "slope", "mean_k"])


def tom_dissimilarity(X, beta):
    A = np.abs(np.corrcoef(X.T)).astype(np.float64) ** beta
    np.fill_diagonal(A, 0)
    L = A @ A
    k = A.sum(axis=0)
    T = (L + A) / (np.minimum.outer(k, k) + 1 - A)
    np.fill_diagonal(T, 1)
    return 1 - T


def tree_cut(Z, n, min_size=30, gap=0.1):
    """Simplified recursive tree cut (own rule). Start from branches below 0.99*max height (static cut used by
    WGCNA's hybrid method). A branch is split into its two children when both children contain >= min_size
    objects and the branch's merge height exceeds both children's merge heights by more than `gap`
    (in 1-TOM units); otherwise the branch is a module (if >= min_size) or unassigned (label 0)."""
    root = hierarchy.to_tree(Z)
    hmax = Z[:, 2].max()
    labels = np.zeros(n, int); nxt = [1]

    def assign(node):
        ids = node.pre_order()
        if len(ids) >= min_size:
            labels[ids] = nxt[0]; nxt[0] += 1

    def rec(node):
        if node.is_leaf():
            return
        if node.dist > 0.99 * hmax:
            rec(node.left); rec(node.right); return
        l, r = node.left, node.right
        if l.count >= min_size and r.count >= min_size and \
           node.dist - max(l.dist, r.dist) > gap:
            rec(l); rec(r)
        else:
            assign(node)

    rec(root)
    return labels


def eigengenes(X, labels):
    """X: samples x features. PC1 of scaled features per module (sign aligned with mean expression)."""
    out = {}
    Xs = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)
    for m in sorted(set(labels) - {0}):
        sub = Xs[:, labels == m]
        u, s, vt = np.linalg.svd(sub, full_matrices=False)
        pc = u[:, 0]
        if np.corrcoef(pc, sub.mean(axis=1))[0, 1] < 0:
            pc = -pc
        out[m] = pc
    return pd.DataFrame(out)


def merge_modules(X, labels, merge_cut=0.25):
    labels = labels.copy()
    while True:
        ME = eigengenes(X, labels)
        if ME.shape[1] < 2:
            return labels
        C = np.corrcoef(ME.values.T); D = 1 - C
        np.fill_diagonal(D, np.inf)
        i, j = np.unravel_index(np.argmin(D), D.shape)
        if D[i, j] >= merge_cut:
            return labels
        a, b = ME.columns[i], ME.columns[j]
        labels[labels == b] = a


def wgcna_modules(X, beta, min_size=30, merge_cut=0.25, gap=0.1):
    D = tom_dissimilarity(X, beta)
    Z = hierarchy.linkage(squareform(D, checks=False), method="average")
    lab = tree_cut(Z, X.shape[1], min_size=min_size, gap=gap)
    lab = merge_modules(X, lab, merge_cut)
    # relabel by size
    sizes = pd.Series(lab[lab > 0]).value_counts()
    remap = {m: i + 1 for i, m in enumerate(sizes.index)}
    return np.array([remap.get(l, 0) for l in lab])


# ------------------------------------------------------------------ kinase-substrate mapping
def build_regulons(site_ids, ks):
    """site_ids like GENE_S123; ks: DataFrame kinase, substrate, site. Exact symbol+site matching."""
    key = pd.Series(range(len(site_ids)), index=site_ids)
    reg = {}
    pairs = ks.assign(sid=ks["substrate"].astype(str) + "_" + ks["site"].astype(str))
    pairs = pairs[pairs["sid"].isin(key.index)]
    for kin, g in pairs.groupby("kinase"):
        reg[kin] = np.unique(key[g["sid"]].values)
    return reg


def activity_mean(M, reg):
    """M: sites x samples (centred log2). Returns kinases x samples mean log2 (PhosMap 'mean', log scale)."""
    return pd.DataFrame({k: M[idx].mean(axis=0) for k, idx in reg.items()}).T


def _kses(sorted_vals, member_mask):
    """PhosMap compute_kses: KS-type deviation with |value|-weighted hits (vectorised over rows of mask)."""
    m = np.atleast_2d(member_mask)
    hit = np.abs(sorted_vals)[None, :] * m
    ph = np.cumsum(hit, axis=1) / np.maximum(hit.sum(axis=1, keepdims=True), 1e-300)
    miss = ~m
    pm = np.cumsum(miss, axis=1) / np.maximum(miss.sum(axis=1, keepdims=True), 1)
    d = ph - pm
    i = np.argmax(np.abs(d), axis=1)
    return d[np.arange(d.shape[0]), i]


def activity_ksea(M, reg, nperm=1000, seed=2024):
    """Per sample: sites ranked by value (descending); ES for each kinase's substrates; P from nperm random
    site sets of the same size (null shared across kinases with equal set size); P = 0 replaced by
    0.1 * minimum non-zero P (as PhosMap); activity = -log10(P) * sign(mean substrate value)."""
    rng = np.random.default_rng(seed)
    G, S = M.shape
    kin = list(reg); out = np.full((len(kin), S), np.nan)
    sizes = sorted({len(reg[k]) for k in kin})
    for s in range(S):
        v = M[:, s]; o = np.argsort(-v, kind="stable"); vs = v[o]
        pos = np.empty(G, int); pos[o] = np.arange(G)
        nulls = {}
        permall = np.argsort(rng.random((nperm, G)), axis=1)[:, :max(sizes)]  # one random ordering per sample;
        for n in sizes:                                                        # size-n sets = first n positions
            mask = np.zeros((nperm, G), bool); mask[np.arange(nperm)[:, None], permall[:, :n]] = True
            nulls[n] = _kses(vs, mask)
        pv = np.zeros(len(kin)); sign = np.zeros(len(kin))
        for i, k in enumerate(kin):
            mask = np.zeros(G, bool); mask[pos[reg[k]]] = True
            es = _kses(vs, mask)[0]; nl = nulls[len(reg[k])]
            pv[i] = np.mean(nl >= es) if es >= 0 else np.mean(nl < es)
            sign[i] = 1 if v[reg[k]].mean() > 0 else -1
        nz = pv[pv > 0]
        pv[pv == 0] = (nz.min() if nz.size else 1 / nperm) * 0.1
        out[:, s] = -np.log10(pv) * sign
    return pd.DataFrame(out, index=kin, columns=None)


def activity_mlr(M, reg, nfold=10, seed=2024, nlambda=100):
    """Ridge regression without intercept, per sample: y (sites with >=1 kinase; values averaged) ~ X (binary
    site x kinase), predictors standardised, lambda chosen by nfold CV (min MSE) - glmnet(alpha=0) analogue."""
    rng = np.random.default_rng(seed)
    kin = sorted(reg)
    sites = np.unique(np.concatenate([reg[k] for k in kin]))
    X = np.zeros((len(sites), len(kin)))
    sidx = {s: i for i, s in enumerate(sites)}
    for j, k in enumerate(kin):
        X[[sidx[s] for s in reg[k]], j] = 1
    mu, sd = X.mean(axis=0), X.std(axis=0); sd[sd == 0] = 1
    Xs = (X - mu) / sd
    folds = rng.integers(0, nfold, len(sites))
    out = np.zeros((len(kin), M.shape[1]))
    for s in range(M.shape[1]):
        y = M[sites, s]; ysd = y.std(); ys = y / ysd
        lmax = np.max(np.abs(Xs.T @ ys)) / len(ys) * 1000  # ridge: glmnet uses lambda_max of alpha=0.001
        lams = np.exp(np.linspace(np.log(lmax), np.log(lmax * 1e-4), nlambda))
        cv = np.zeros(len(lams))
        for f in range(nfold):
            tr, te = folds != f, folds == f
            A = Xs[tr]; U, d, Vt = np.linalg.svd(A, full_matrices=False); uy = U.T @ ys[tr]
            for li, lam in enumerate(lams):
                b = Vt.T @ (d / (d ** 2 + lam * tr.sum()) * uy)
                cv[li] += np.sum((ys[te] - Xs[te] @ b) ** 2)
        lam = lams[np.argmin(cv)]
        U, d, Vt = np.linalg.svd(Xs, full_matrices=False)
        b = Vt.T @ (d / (d ** 2 + lam * len(ys)) * (U.T @ ys))
        out[:, s] = b / sd * ysd
    return pd.DataFrame(out, index=kin)


# ------------------------------------------------------------------ INKA
def inka_scores(intens, site_gene, site_acc_pos, kinases, act_loop, psp, nwk):
    """intens: sites x samples (linear intensities). site_gene: gene symbol per site; site_acc_pos:
    'ACC-Ypos' keys; kinases: set of kinase symbols; act_loop: set of (gene, residue-position) activation-loop
    sites; psp / nwk: DataFrames (kinase, site_key) with site_key 'ACC-Ypos'.
    Returns dict of DataFrames (kinase x sample) for the four components and the INKA score."""
    intens = np.asarray(intens, float)
    S = intens.shape[1]
    gene = np.asarray(site_gene); keys = pd.Index(site_acc_pos)
    kin = sorted(set(kinases))
    comp = {c: pd.DataFrame(0.0, index=kin, columns=range(S)) for c in ["kinome", "actloop", "psp", "nwk"]}
    pos_res = np.array([k.split("-")[-1] for k in site_acc_pos])
    for k in kin:
        own = gene == k
        if own.any():
            comp["kinome"].loc[k] = intens[own].sum(axis=0)
            loop = own & np.isin(pos_res, [p for (g, p) in act_loop if g == k])
            comp["actloop"].loc[k] = intens[loop].sum(axis=0)
    for name, tab in [("psp", psp), ("nwk", nwk)]:
        t = tab[tab["site_key"].isin(keys) & tab["kinase"].isin(kin)]
        loc = keys.get_indexer(t["site_key"])
        for k, g in pd.DataFrame({"k": t["kinase"].values, "i": loc}).groupby("k"):
            comp[name].loc[k] = intens[np.unique(g["i"].values)].sum(axis=0)
    kc = comp["kinome"] + comp["actloop"]; sc = comp["psp"] + comp["nwk"]
    comp["kinase_centric"] = kc; comp["substrate_centric"] = sc
    comp["inka"] = np.sqrt(kc * sc)
    return comp
