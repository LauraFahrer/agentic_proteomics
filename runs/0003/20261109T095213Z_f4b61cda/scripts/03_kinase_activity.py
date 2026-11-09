"""Differential kinase-activity inference with phosX (motif-based, Kinase Library PSSMs).

Reads the phosX outputs produced by run_phosx.sh for each contrast and layer:
  results/tables/kinase/<contrast>_{ST,Y}_kinases.tsv
phosX columns: KS (signed enrichment), 'p value', 'FDR q value', 'Activity Score' (= -log10 p * sign).
Builds a summary table and a figure of the most differentially active kinases.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CAT, INK2, TAB, savefig, style

style()
KIN = TAB / "kinase"
out = {}


def read_k(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "kinase"
    # phosX Activity Score (-log10 p) saturates at the permutation floor; rank by the
    # continuous KS enrichment statistic instead (sign = direction, |KS| = strength).
    return df


def top(df, n=12, direction=1):
    d = df[df["KS"] * direction > 0].copy()
    return d.sort_values("KS", ascending=(direction < 0)).head(n)


panels = []
for contrast in ["EMT_vs_Metabolism", "EBV_vs_nonEBV"]:
    for layer in ["ST", "Y"]:
        f = KIN / f"{contrast}_{layer}_kinases.tsv"
        if not f.exists():
            continue
        df = read_k(f)
        df.to_csv(TAB / f"kinase_activity_{contrast}_{layer}.csv")
        sig = df[df["FDR q value"] < 0.1]
        out[f"{contrast}_{layer}"] = {
            "kinases_scored": int(len(df)), "q<0.1": int(len(sig)),
            "up_g1_byKS": top(df, 12, 1).index.tolist(), "up_g2_byKS": top(df, 12, -1).index.tolist(),
            "sig_up_g1": sig[sig.KS > 0].sort_values("KS", ascending=False).index.tolist()[:12],
            "sig_up_g2": sig[sig.KS < 0].sort_values("KS").index.tolist()[:12]}
        panels.append((contrast, layer, df))
print(json.dumps(out, indent=2))
json.dump(out, open(TAB / "kinase_activity_summary.json", "w"), indent=2)

# figure: for the primary contrast, S/T and Y kinase activities (top differential)
fig, axes = plt.subplots(1, 2, figsize=(12.5, 7.5), gridspec_kw={"wspace": 0.55})
for ax, layer, lab in [(axes[0], "ST", "Ser/Thr kinome"), (axes[1], "Y", "Tyr kinome")]:
    f = KIN / f"EMT_vs_Metabolism_{layer}_kinases.tsv"
    if not f.exists():
        ax.set_visible(False); continue
    df = read_k(f)
    sel = pd.concat([top(df, 14, 1), top(df, 14, -1)]).sort_values("KS")
    colors = [CAT[1] if v > 0 else CAT[0] for v in sel["KS"]]
    ax.barh(range(len(sel)), sel["KS"], color=colors, height=0.74)
    ax.set_yticks(range(len(sel)), sel.index, fontsize=7)
    ax.axvline(0, color="#c3c2b7", lw=0.6)
    for i, (v, q) in enumerate(zip(sel["KS"], sel["FDR q value"])):
        if q < 0.1:
            ax.text(v + (0.01 if v > 0 else -0.01), i, "*", va="center", ha="left" if v > 0 else "right",
                    color=INK2, fontsize=11)
    ax.set_xlabel("phosX kinase enrichment (KS)\n(+ = higher in EMT, − = higher in Metabolism)", fontsize=8)
    ax.set_title(f"{lab}  (* = FDR<0.1)", fontsize=10)
fig.suptitle("Differential kinase activity between GC proteome subtypes (phosX / Kinase Library)",
             fontsize=11, y=1.0)
savefig(fig, "03_kinase_activity_subtype")
plt.close(fig)
