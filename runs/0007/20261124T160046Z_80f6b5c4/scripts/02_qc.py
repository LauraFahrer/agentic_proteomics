#!/usr/bin/env python3
"""Step 02 - replicate correlation (Fig. S1C) and GO cellular-component distribution (Fig. 1B).
Pearson correlation: numpy (general-purpose). GO annotations: QuickGO (external database)."""
import json
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from common import *

WP, wa = load_WP(); PP, pa = load_PP(); YP, ya = load_YP(); WPq, _ = load_WPq()
rows = []
fig, axes = plt.subplots(2, 2, figsize=(14, 13))
for ax, (name, m, ref) in zip(axes.flat, [("WP (Table S3)", WP, 0.93), ("PP (Table S3)", PP, 0.82),
                                          ("YP (Table S3)", YP, 0.72), ("WP (quantms, reprocessed)", WPq, np.nan)]):
    C = np.corrcoef(m.values.T)
    ln = np.array([line_of(c) for c in m.columns])
    same = (ln[:, None] == ln[None, :]) & np.triu(np.ones_like(C, bool), 1)
    rows.append(dict(dataset=name, mean_all_pairs=C[np.triu_indices_from(C, 1)].mean(),
                     mean_within_line_replicates=C[same].mean(), min_within_line_replicates=C[same].min(),
                     paper_mean_pearson=ref))
    im = ax.imshow(C, cmap="RdYlBu_r"); ax.set_title(name)
    ax.set_xticks(range(len(m.columns))); ax.set_xticklabels(m.columns, rotation=90, fontsize=5)
    ax.set_yticks(range(len(m.columns))); ax.set_yticklabels(m.columns, fontsize=5)
    fig.colorbar(im, ax=ax, shrink=0.7)
savefig(fig, "02_FigS1C_replicate_correlation")
cs = pd.DataFrame(rows); cs.to_csv(os.path.join(TAB, "02_replicate_correlation.csv"), index=False); print(cs)

terms = {"Nucleus": "GO:0005634", "Nucleolus": "GO:0005730", "Extracellular matrix": "GO:0031012",
         "Ribosome": "GO:0005840", "Endoplasmic reticulum": "GO:0005783", "Golgi apparatus": "GO:0005794",
         "Cytosol": "GO:0005829", "Membrane": "GO:0016020", "Cytoskeleton": "GO:0005856",
         "Centrosome": "GO:0005813", "Extracellular exosome": "GO:0070062", "Lysosome": "GO:0005764",
         "Peroxisome": "GO:0005777", "Mitochondrion": "GO:0005739"}
paper = {"Nucleus": (2307, 2206), "Nucleolus": (493, 464), "Extracellular matrix": (129, 100), "Ribosome": (138, 99),
         "Endoplasmic reticulum": (751, 15), "Golgi apparatus": (581, 350), "Cytosol": (1786, 1768),
         "Membrane": (1785, 1627), "Cytoskeleton": (434, 512), "Centrosome": (172, 231),
         "Extracellular exosome": (1523, 1172), "Lysosome": (229, 80), "Peroxisome": (88, 64), "Mitochondrion": (757, 476)}
wpg = set(wa["Gene"]); ppg = set(pa["GeneSymbol"]) | set(ya["GeneSymbol"])
out = []
for t, go in terms.items():
    sym = set(pd.read_csv(os.path.join(DC, "go", go.replace(":", "_") + "_exact.tsv"), sep="\t")["SYMBOL"])
    out.append(dict(term=t, GO=go, proteins=len(wpg & sym), pproteins=len(ppg & sym),
                    paper_proteins=paper[t][0], paper_pproteins=paper[t][1]))
go = pd.DataFrame(out); go.to_csv(os.path.join(TAB, "02_Fig1B_GO_CC_distribution.csv"), index=False)
from scipy.stats import spearmanr
rho = (spearmanr(go.proteins, go.paper_proteins)[0], spearmanr(go.pproteins, go.paper_pproteins)[0])
print(go, "\nSpearman", rho)
json.dump({"spearman_proteins": rho[0], "spearman_pproteins": rho[1]}, open(os.path.join(TAB, "02_Fig1B_spearman.json"), "w"), indent=1)
fig, axs = plt.subplots(1, 2, figsize=(12, 6))
o = go.sort_values("paper_proteins")
for ax, (a, b, col, ttl) in zip(axs, [("paper_proteins", "proteins", "#3b7dd8", "Proteins (WP)"),
                                      ("paper_pproteins", "pproteins", "#d83b3b", "Phosphoproteins (PP+YP)")]):
    y = np.arange(len(o)); ax.barh(y - 0.2, o[a], 0.4, color="grey", label="paper")
    ax.barh(y + 0.2, o[b], 0.4, color=col, label="this run (QuickGO exact term)")
    ax.set_yticks(y); ax.set_yticklabels(o.term); ax.set_title(ttl); ax.legend(fontsize=7)
savefig(fig, "02_Fig1B_GO_CC_distribution")
