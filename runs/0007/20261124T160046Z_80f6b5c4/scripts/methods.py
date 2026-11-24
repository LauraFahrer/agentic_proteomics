"""Self-implemented analytical methods for run 0007 (independent reproduction of Hu et al. 2024).

Only general-purpose libraries are used: numpy (linear algebra, random numbers), scipy (special functions,
distributions, hierarchical clustering primitives), pandas (tables). Every *specialised* method used in the
paper is re-implemented here from its original publication:

  quantile normalisation ............ Bolstad et al. 2003
  NMF (Brunet / KL updates) ......... Brunet et al. 2004 PNAS; stopping rule as in Gaujoux & Seoighe 2010 (R 'NMF')
  consensus, cophenetic, silhouette . Brunet 2004; Rousseeuw 1987
  metagene feature extraction ....... Kim & Park 2007 (score + selection rule)
  consensus clustering .............. Monti et al. 2003 (as in ConsensusClusterPlus, Wilkerson & Hayes 2010)
  moderated t-statistics ............ Smyth 2004 (limma lmFit/eBayes, incl. fitFDist)
  GSEA (S2N, weighted KS, NES, FDR) . Subramanian et al. 2005 PNAS
  ssGSEA / PTM-SEA .................. Barbie et al. 2009; Krug et al. 2019 MCP
  WGCNA-type modules ................ Zhang & Horvath 2005; Langfelder & Horvath 2008 (simplified tree cut)
  KSEA / mean / MLR kinase activity . Casado et al. 2013; PhosMap (Zhan et al. 2023) definitions
  INKA score ........................ Beekhof et al. 2019 Mol Syst Biol
"""
import numpy as np
import pandas as pd
from scipy import special, stats
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform


# ============================================================ preprocessing
def quantile_normalize(df):
    """Missing-value-aware quantile normalisation (preprocessCore-like: mean of rank-interpolated columns)."""
    a = df.values.astype(float)
    L = a.shape[0]
    grid = np.linspace(0, 1, L)
    ref = np.zeros(L)
    for j in range(a.shape[1]):
        v = np.sort(a[~np.isnan(a[:, j]), j])
        ref += np.interp(grid, np.linspace(0, 1, len(v)), v)
    ref /= a.shape[1]
    out = np.full_like(a, np.nan)
    for j in range(a.shape[1]):
        ok = ~np.isnan(a[:, j])
        r = stats.rankdata(a[ok, j])
        out[ok, j] = np.interp((r - 1) / max(ok.sum() - 1, 1), grid, ref)
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def impute_downshift(df, shift=1.8, width=0.3, seed=1):
    rng = np.random.default_rng(seed)
    out = df.copy()
    for c in out.columns:
        v = out[c]; na = v.isna()
        out.loc[na, c] = rng.normal(v.mean() - shift * v.std(), width * v.std(), na.sum())
    return out


def mean_abs_dev(m):
    return np.mean(np.abs(m - m.mean(axis=1, keepdims=True)), axis=1)


def median_abs_dev(m):
    med = np.median(m, axis=1, keepdims=True)
    return 1.4826 * np.median(np.abs(m - med), axis=1)


def bh(p):
    p = np.asarray(p, float); n = np.sum(~np.isnan(p))
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p); o = np.argsort(p[ok]); ranked = p[ok][o] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty(n); tmp[o] = np.minimum(ranked, 1); out[ok] = tmp
    return out


# ============================================================ NMF (Brunet, KL divergence)
EPS = np.finfo(float).eps


def nmf_brunet(V, k, rng, max_iter=2000, check_every=10, stationary=40):
    """Brunet et al. multiplicative updates minimising KL divergence. Random uniform init in [0, max(V)]
    (as the R 'NMF' package 'random' seeding). Stops when the sample connectivity (argmax of H) is unchanged
    for `stationary` consecutive checks, or at max_iter (R NMF default stopping criterion for brunet)."""
    n, m = V.shape
    vmax = V.max()
    W = rng.uniform(0, vmax, (n, k)); H = rng.uniform(0, vmax, (k, m))
    last = None; same = 0
    for it in range(1, max_iter + 1):
        WH = W @ H + EPS
        H *= (W.T @ (V / WH)) / (W.sum(axis=0)[:, None] + EPS)
        WH = W @ H + EPS
        W *= ((V / WH) @ H.T) / (H.sum(axis=1)[None, :] + EPS)
        if it % check_every == 0:
            cl = H.argmax(axis=0)
            if last is not None and np.array_equal(cl, last):
                same += 1
                if same >= stationary:
                    break
            else:
                same = 0
            last = cl
    WH = W @ H + EPS
    dev = np.sum(V * np.log((V + EPS) / WH) - V + WH)
    return W, H, dev


def nmf_consensus(V, k, nrun=100, seed=123456):
    rng = np.random.default_rng(seed)
    m = V.shape[1]
    C = np.zeros((m, m)); best = None
    for _ in range(nrun):
        W, H, dev = nmf_brunet(V, k, rng)
        cl = H.argmax(axis=0)
        C += (cl[:, None] == cl[None, :])
        if best is None or dev < best[2]:
            best = (W, H, dev)
    C /= nrun
    return C, best


def consensus_hclust(C, k, method="average"):
    D = 1 - C; np.fill_diagonal(D, 0)
    Z = hierarchy.linkage(squareform(D, checks=False), method=method)
    return Z, hierarchy.fcluster(Z, k, criterion="maxclust")


def cophenetic_coef(C):
    D = 1 - C; np.fill_diagonal(D, 0)
    d = squareform(D, checks=False)
    Z = hierarchy.linkage(d, method="average")
    return np.corrcoef(hierarchy.cophenet(Z), d)[0, 1]


def silhouette(D, labels):
    """Rousseeuw silhouette widths from a distance matrix."""
    labels = np.asarray(labels); s = np.zeros(len(labels))
    for i in range(len(labels)):
        own = labels == labels[i]; own[i] = False
        if own.sum() == 0:
            continue
        a = D[i, own].mean()
        b = min(D[i, labels == c].mean() for c in np.unique(labels) if c != labels[i])
        s[i] = (b - a) / max(a, b)
    return s


def kim_park_features(W):
    """Kim & Park (2007) feature score and selection, as used by NMF::extractFeatures(method='kim')."""
    P = W / (W.sum(axis=1, keepdims=True) + EPS)
    k = W.shape[1]
    score = 1 + np.sum(P * np.log2(P + EPS), axis=1) / np.log2(k)
    thr = np.median(score) + 3 * 1.4826 * np.median(np.abs(score - np.median(score)))
    sel = []
    for q in range(k):
        keep = (score > thr) & (W.argmax(axis=1) == q) & (W[:, q] > np.median(W[:, q]))
        sel.append(np.where(keep)[0])
    return score, sel


# ============================================================ consensus clustering (Monti / CCP)
def consensus_clustering(X, kmax=4, reps=500, p_item=0.8, seed=123456, linkage="average"):
    """X: features x samples. Inner clustering = hierarchical (average) on 1 - Pearson correlation between
    samples, item subsampling p_item, all features (ConsensusClusterPlus defaults: clusterAlg='hc',
    distance='pearson', innerLinkage='average', finalLinkage='average')."""
    rng = np.random.default_rng(seed)
    m = X.shape[1]
    res = {}
    Ms = {k: np.zeros((m, m)) for k in range(2, kmax + 1)}
    Is = np.zeros((m, m))
    for _ in range(reps):
        idx = np.sort(rng.choice(m, int(np.floor(m * p_item)), replace=False))
        sub = X[:, idx]
        sd = sub.std(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            R = np.corrcoef(sub.T)
        R = np.nan_to_num(R, nan=0.0)
        D = 1 - R; np.fill_diagonal(D, 0); D = np.clip(D, 0, None)
        Z = hierarchy.linkage(squareform((D + D.T) / 2, checks=False), method=linkage)
        Is[np.ix_(idx, idx)] += 1
        for k in range(2, kmax + 1):
            cl = hierarchy.fcluster(Z, k, criterion="maxclust")
            Ms[k][np.ix_(idx, idx)] += (cl[:, None] == cl[None, :])
    for k in range(2, kmax + 1):
        C = Ms[k] / np.maximum(Is, 1); np.fill_diagonal(C, 1)
        Z, cl = consensus_hclust(C, k, linkage)
        res[k] = {"consensus": C, "classes": cl}
    return res


# ============================================================ limma-type moderated statistics
def _trigamma_inverse(x):
    """Smyth's Newton iteration for the inverse trigamma function."""
    x = np.asarray(x, float)
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.all(np.abs(-dif / y) < 1e-8):
            break
    return y


def fit_f_dist(s2, df1):
    """limma::fitFDist - moment estimation of prior d0, s0^2 from log residual variances."""
    s2 = np.asarray(s2, float); ok = (s2 > 0) & np.isfinite(s2)
    z = np.log(s2[ok]); e = z - special.digamma(df1 / 2) + np.log(df1 / 2)
    emean = e.mean(); evar = np.sum((e - emean) ** 2) / (ok.sum() - 1) - special.polygamma(1, df1 / 2)
    if evar > 0:
        d0 = 2 * _trigamma_inverse(evar)
        s02 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0 = np.inf; s02 = np.exp(emean)
    return d0, s02


def moderated_t(Y, design, contrast):
    """lmFit + contrasts.fit + eBayes for a full-rank design (no missing values). Returns logFC, t, P, adjP."""
    X = np.asarray(design, float); Y = np.asarray(Y, float)
    XtX_inv = np.linalg.inv(X.T @ X)
    B = Y @ X @ XtX_inv
    resid = Y - B @ X.T
    df = X.shape[0] - X.shape[1]
    s2 = np.sum(resid ** 2, axis=1) / df
    c = np.asarray(contrast, float)
    est = B @ c; unscaled = np.sqrt(c @ XtX_inv @ c)
    d0, s02 = fit_f_dist(s2, df)
    s2post = (d0 * s02 + df * s2) / (d0 + df) if np.isfinite(d0) else np.full_like(s2, s02)
    dft = df + d0 if np.isfinite(d0) else np.inf
    t = est / (np.sqrt(s2post) * unscaled)
    p = 2 * stats.t.sf(np.abs(t), dft) if np.isfinite(dft) else 2 * stats.norm.sf(np.abs(t))
    return pd.DataFrame({"logFC": est, "t": t, "P": p, "adjP": bh(p)}), (d0, s02)


def oneway_anova(Y, groups):
    groups = np.asarray(groups); levels = np.unique(groups)
    p = np.array([stats.f_oneway(*[y[groups == g] for g in levels]).pvalue for y in Y])
    return p


# ============================================================ gene-set statistics
def signal2noise(X, lab):
    """GSEA Signal2Noise with the GSEA sigma floor (sigma >= 0.2*|mu|, min 0.2)."""
    a = X[:, lab]; b = X[:, ~lab]
    ma, mb = a.mean(axis=1), b.mean(axis=1)
    sa, sb = a.std(axis=1, ddof=1), b.std(axis=1, ddof=1)
    sa = np.maximum(sa, 0.2 * np.abs(ma)); sb = np.maximum(sb, 0.2 * np.abs(mb))
    sa = np.where(ma == 0, np.maximum(sa, 0.2), sa); sb = np.where(mb == 0, np.maximum(sb, 0.2), sb)
    return (ma - mb) / (sa + sb)


def _es_matrix(metric, M):
    """Weighted (p=1) KS enrichment scores for all sets. metric: (G,), M: sets x G boolean membership."""
    o = np.argsort(-metric, kind="stable")
    r = metric[o]; Ms = M[:, o]
    w = np.abs(r)[None, :] * Ms
    nr = w.sum(axis=1, keepdims=True); nh = Ms.shape[1] - Ms.sum(axis=1, keepdims=True)
    run = np.cumsum(w / np.maximum(nr, EPS) - (~Ms) / nh, axis=1)
    mx, mn = run.max(axis=1), run.min(axis=1)
    return np.where(np.abs(mx) > np.abs(mn), mx, mn)


def gsea(X, genes, lab, gene_sets, nperm=1000, permute="phenotype", min_size=15, max_size=500, seed=149):
    """GSEA (Subramanian et al. 2005). X: genes x samples (already collapsed to unique genes),
    lab: boolean vector (True = phenotype of interest). Returns ES, NES, nominal P and FDR q per set."""
    rng = np.random.default_rng(seed)
    genes = np.asarray(genes); gidx = {g: i for i, g in enumerate(genes)}
    names, rows = [], []
    for n, gs in gene_sets.items():
        ii = [gidx[g] for g in set(gs) if g in gidx]
        if min_size <= len(ii) <= max_size:
            r = np.zeros(len(genes), bool); r[ii] = True; names.append(n); rows.append(r)
    M = np.array(rows)
    metric = signal2noise(X, lab)
    es = _es_matrix(metric, M)
    null = np.zeros((len(names), nperm))
    for p in range(nperm):
        if permute == "phenotype":
            null[:, p] = _es_matrix(signal2noise(X, rng.permutation(lab)), M)
        else:  # gene_set permutation: random sets of equal size on the observed ranking
            Mp = np.zeros_like(M)
            for i in range(M.shape[0]):
                Mp[i, rng.choice(M.shape[1], M[i].sum(), replace=False)] = True
            null[:, p] = _es_matrix(metric, Mp)
    pos_mean = np.array([null[i][null[i] >= 0].mean() if np.any(null[i] >= 0) else np.nan for i in range(len(names))])
    neg_mean = np.array([np.abs(null[i][null[i] < 0].mean()) if np.any(null[i] < 0) else np.nan for i in range(len(names))])
    nes = np.where(es >= 0, es / pos_mean, es / neg_mean)
    nnull = np.where(null >= 0, null / pos_mean[:, None], null / neg_mean[:, None])
    pnom = np.array([(np.sum(null[i] >= es[i]) / max(np.sum(null[i] >= 0), 1)) if es[i] >= 0 else
                     (np.sum(null[i] <= es[i]) / max(np.sum(null[i] < 0), 1)) for i in range(len(names))])
    fdr = np.zeros(len(names))
    allpos = nnull[nnull >= 0]; allneg = nnull[nnull < 0]
    for i in range(len(names)):
        if nes[i] >= 0:
            a = np.mean(allpos >= nes[i]) if allpos.size else 0; b = np.mean(nes[nes >= 0] >= nes[i])
        else:
            a = np.mean(allneg <= nes[i]) if allneg.size else 0; b = np.mean(nes[nes < 0] <= nes[i])
        fdr[i] = min(1.0, a / b) if b > 0 else 1.0
    return pd.DataFrame({"set": names, "size": M.sum(axis=1), "ES": es, "NES": nes, "NOM_p": pnom, "FDR_q": fdr})


def collapse_max_probe(df, gene_of):
    """GSEA 'Max_probe' collapse: per gene and sample, the maximum over its probes/sites."""
    g = pd.Series(gene_of, index=df.index)
    return df.groupby(g.values).max()


def ora(hits, universe, gene_sets, min_size=10, max_size=500):
    """Hypergeometric over-representation test with BH adjustment."""
    hits = set(hits) & set(universe); U = set(universe); N = len(U); n = len(hits)
    out = []
    for name, gs in gene_sets.items():
        s = set(gs) & U
        if not (min_size <= len(s) <= max_size):
            continue
        x = len(s & hits)
        out.append((name, len(s), x, stats.hypergeom.sf(x - 1, N, len(s), n) if x > 0 else 1.0))
    r = pd.DataFrame(out, columns=["set", "size", "overlap", "P"])
    r["adjP"] = bh(r["P"].values); r["n_hits"] = n
    return r.sort_values("P")


def read_gmt(path):
    gs = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        gs[f[0]] = [x for x in f[2:] if x]
    return gs


def ssgsea_scores(X, ids, signatures, weight=0.75, nperm=1000, min_overlap=10, seed=1):
    """PTM-SEA/ssGSEA2 (Krug 2019): per sample, features ranked by value, weighted running sum
    (weight 0.75), enrichment = area under the running-sum difference (ssGSEA, Barbie 2009).
    Directional signatures: '<id>;u' sites counted in the up set, '<id>;d' in the down set; the score is
    ES(up) - ES(down). NES = ES / sd(null ES) from feature-label permutations; P from permutations."""
    rng = np.random.default_rng(seed)
    idx = {f: i for i, f in enumerate(ids)}
    sets = {}
    for name, members in signatures.items():
        up = [idx[m[:-2]] for m in members if m.endswith(";u") and m[:-2] in idx]
        dn = [idx[m[:-2]] for m in members if m.endswith(";d") and m[:-2] in idx]
        if len(up) + len(dn) >= min_overlap:
            sets[name] = (np.array(up, int), np.array(dn, int))
    G, S = X.shape

    def es_area(vals, members):
        # members: (G,) or (n, G) boolean in descending-rank order; returns area under P_hit - P_miss
        mem = np.atleast_2d(members)
        w = np.abs(vals)[None, :] ** weight
        hw = w * mem
        tot = hw.sum(axis=1, keepdims=True)
        ph = np.cumsum(hw, axis=1) / np.where(tot > 0, tot, 1)
        miss = ~mem
        pm = np.cumsum(miss, axis=1) / np.maximum(miss.sum(axis=1, keepdims=True), 1)
        out = np.sum(ph - pm, axis=1)
        out[tot[:, 0] == 0] = 0
        return out if members.ndim == 2 else out[0]

    names = list(sets)
    ES = np.zeros((len(names), S)); P = np.ones((len(names), S)); NES = np.zeros((len(names), S))
    for s in range(S):
        v = X[:, s]; ok = ~np.isnan(v)
        rk = np.zeros(G); rk[ok] = stats.rankdata(v[ok])  # sample normalisation: rank (NaN -> bottom)
        o = np.argsort(-rk, kind="stable"); vals = rk[o]
        pos = np.empty(G, int); pos[o] = np.arange(G)
        kmax = max(len(u) + len(d) for u, d in sets.values())
        permall = np.argsort(rng.random((nperm, G)), axis=1)[:, :kmax]
        for i, n in enumerate(names):
            up, dn = sets[n]
            mu = np.zeros(G, bool); mu[pos[up]] = True
            md = np.zeros(G, bool); md[pos[dn]] = True
            e = es_area(vals, mu) - es_area(vals, md)
            k = len(up) + len(dn)
            perm = permall[:, :k]   # one random ordering per sample, prefixes give random sets of size k
            a = np.zeros((nperm, G), bool); b = np.zeros((nperm, G), bool)
            rows = np.arange(nperm)[:, None]
            if len(up): a[rows, perm[:, :len(up)]] = True
            if len(dn): b[rows, perm[:, len(up):]] = True
            nulls = es_area(vals, a) - es_area(vals, b)
            ES[i, s] = e; NES[i, s] = e / nulls.std()
            P[i, s] = (np.sum(np.abs(nulls) >= abs(e)) + 1) / (nperm + 1)
    return names, ES, NES, P
