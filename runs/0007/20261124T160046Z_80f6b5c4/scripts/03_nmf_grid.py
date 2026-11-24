#!/usr/bin/env python3
"""Step 03 - NMF subtype grid with the self-implemented Brunet NMF (methods.nmf_consensus).
Grid: top 10/20/30/50/100 % features by mean absolute deviation (Methods) x k = 2..5, 100 runs each.
Usage: python 03_nmf_grid.py [WP PP YP WPq] [--mad median]"""
import os
import sys
import numpy as np
import pandas as pd
from multiprocessing import Pool
from common import *
import methods as M

ND = os.path.join(DC, "nmf"); os.makedirs(ND, exist_ok=True)
FRACS = [0.1, 0.2, 0.3, 0.5, 1.0]; KS = [2, 3, 4, 5]
madtype = "median" if "--mad" in sys.argv and sys.argv[sys.argv.index("--mad") + 1] == "median" else "mean"
dsets = [a for a in sys.argv[1:] if a in ("WP", "PP", "YP", "WPq")] or ["WP", "PP", "YP", "WPq"]
loaders = {"WP": load_WP, "PP": load_PP, "YP": load_YP, "WPq": load_WPq}


def job(args):
    d, fr, k = args
    m, _ = loaders[d]()
    V = m.values
    dev = M.mean_abs_dev(V) if madtype == "mean" else M.median_abs_dev(V)
    n = int(np.ceil(fr * len(V)))
    idx = np.argsort(-dev, kind="stable")[:n]
    C, (W, H, devi) = M.nmf_consensus(V[idx], k, nrun=100, seed=123456)
    Z, cl = M.consensus_hclust(C, k)
    D = 1 - C; np.fill_diagonal(D, 0)
    sil = M.silhouette(D, cl)
    np.savez_compressed(os.path.join(ND, f"{d}_{madtype}_{int(fr*100):03d}_k{k}.npz"), C=C, W=W, H=H,
                        features=np.array(m.index[idx]), samples=np.array(m.columns), cl=cl, sil=sil)
    return dict(dataset=d, mad=madtype, frac=fr, n_features=n, k=k, cophenetic=M.cophenetic_coef(C),
                mean_silhouette=sil.mean(), ARI_vs_paper=ari(cl, paper_label(d, list(m.columns))))


if __name__ == "__main__":
    jobs = [(d, fr, k) for d in dsets for fr in FRACS for k in KS]
    jobs.sort(key=lambda j: -j[1] * (8 if j[0] == "PP" else 5))   # big jobs first
    with Pool(int(os.environ.get("NCORES", 4))) as p:
        res = p.map(job, jobs, chunksize=1)
    r = pd.DataFrame(res).sort_values(["dataset", "frac", "k"])
    for d in dsets:
        r[r.dataset == d].to_csv(os.path.join(TAB, f"03_nmf_grid_{d}_{madtype}.csv"), index=False)
    print(r.to_string())
