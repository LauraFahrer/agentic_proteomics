"""Small statistics helpers: moderated t (limma-style empirical Bayes) and rank-based
gene-set tests. Implemented here to avoid an R dependency."""
import numpy as np
import pandas as pd
from scipy import special, stats

from common import bh


def _fit_prior(s2, d):
    """Estimate prior df d0 and prior variance s0^2 (Smyth 2004, method of moments)."""
    s2 = s2[np.isfinite(s2) & (s2 > 0)]
    z = np.log(s2)
    e = z - special.digamma(d / 2) + np.log(d / 2)
    emean = e.mean()
    evar = np.var(e, ddof=1) - special.polygamma(1, d / 2)
    if evar <= 0:
        return np.inf, np.exp(emean)
    # solve trigamma(d0/2) = evar
    lo, hi = 1e-3, 1e4
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if special.polygamma(1, mid / 2) > evar:
            lo = mid
        else:
            hi = mid
    d0 = mid
    s02 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    return d0, s02


def moderated_t(M, g1, g2, min_per_group=2):
    """Two-group moderated t-test on rows of M (log2). g1/g2: column lists.
    Returns DataFrame with logFC (g1 - g2), t, p, q, n1, n2."""
    A, B = M[g1], M[g2]
    n1, n2 = A.notna().sum(axis=1), B.notna().sum(axis=1)
    ok = (n1 >= min_per_group) & (n2 >= min_per_group)
    A, B, n1, n2 = A[ok], B[ok], n1[ok], n2[ok]
    m1, m2 = A.mean(axis=1), B.mean(axis=1)
    ss = ((A.sub(m1, axis=0)) ** 2).sum(axis=1) + ((B.sub(m2, axis=0)) ** 2).sum(axis=1)
    d = (n1 + n2 - 2).astype(float)
    s2 = ss / d
    dmode = float(d.mode().iloc[0])
    d0, s02 = _fit_prior(s2[d == dmode].values, dmode)
    if np.isinf(d0):
        s2post = pd.Series(s02, index=s2.index)
        dtot = np.full(len(s2), 1e6)
    else:
        s2post = (d0 * s02 + d * s2) / (d0 + d)
        dtot = d0 + d
    se = np.sqrt(s2post * (1 / n1 + 1 / n2))
    t = (m1 - m2) / se
    p = 2 * stats.t.sf(np.abs(t), dtot)
    out = pd.DataFrame({"logFC": m1 - m2, "mean_g1": m1, "mean_g2": m2, "t": t, "p": p,
                        "n1": n1, "n2": n2})
    out["q"] = bh(out.p.values)
    out.attrs["d0"] = d0
    return out


def rank_set_test(stat, sets, min_size=10):
    """For each set: Wilcoxon rank-sum of member statistics vs non-members (two-sided),
    plus mean statistic. stat: Series indexed by gene symbol."""
    stat = stat.dropna()
    stat = stat[~stat.index.duplicated()]
    ranks = stat.rank()
    N = len(stat)
    rows = []
    for name, genes in sets.items():
        mem = stat.index.isin(list(genes))
        n = int(mem.sum())
        if n < min_size or n > N - min_size:
            continue
        r1 = ranks[mem].sum()
        u = r1 - n * (n + 1) / 2
        mu = n * (N - n) / 2
        sd = np.sqrt(n * (N - n) * (N + 1) / 12)
        z = (u - mu) / sd
        rows.append((name, n, float(stat[mem].mean()), float(stat[mem].median()), z,
                     2 * stats.norm.sf(abs(z))))
    df = pd.DataFrame(rows, columns=["set", "n", "mean_stat", "median_stat", "z", "p"])
    df["q"] = bh(df.p.values)
    return df.sort_values("p")
