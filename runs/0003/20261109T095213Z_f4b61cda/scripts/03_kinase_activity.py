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
    return df


def top(df, n=12, direction=1):
    d = df[df["Activity Score"] * direction > 0].copy()
    d["rank"] = d["Activity Score"].abs()
    return d.sort_values("rank", ascending=False).head(n)


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
            "up_g1": top(df, 10, 1).index.tolist(), "up_g2": top(df, 10, -1).index.tolist(),
            "sig_up_g1": sig[sig["Activity Score"] > 0].sort_values("Activity Score", ascending=False).index.tolist()[:12],
            "sig_up_g2": sig[sig["Activity Score"] < 0].sort_values("Activity Score").index.tolist()[:12]}
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
    sel = pd.concat([top(df, 14, 1), top(df, 14, -1)]).sort_values("Activity Score")
    colors = [CAT[1] if v > 0 else CAT[0] for v in sel["Activity Score"]]
    ax.barh(range(len(sel)), sel["Activity Score"], color=colors, height=0.74)
    ax.set_yticks(range(len(sel)), sel.index, fontsize=7)
    ax.axvline(0, color="#c3c2b7", lw=0.6)
    for i, (v, q) in enumerate(zip(sel["Activity Score"], sel["FDR q value"])):
        if q < 0.1:
            ax.text(v + (0.03 if v > 0 else -0.03), i, "*", va="center", ha="left" if v > 0 else "right",
                    color=INK2, fontsize=10)
    ax.set_xlabel("phosX activity score\n(+ = higher in EMT, − = higher in Metabolism)", fontsize=8)
    ax.set_title(f"{lab}  (* = FDR<0.1)", fontsize=10)
fig.suptitle("Differential kinase activity between GC proteome subtypes (phosX / Kinase Library)",
             fontsize=11, y=1.0)
savefig(fig, "03_kinase_activity_subtype")
plt.close(fig)
