"""Differential abundance + pathway enrichment (prompt_0002 core deliverable).

Primary biological contrast = the two proteome subtypes (Hu et al.): EMT vs Metabolism
(7 vs 7 cell lines; cell-line means are the independent unit). Differential proteins are
found with a moderated t-test; the pathways associated with the changes are characterised
by (a) preranked GSEA on the t-statistic across Hallmark/KEGG/Reactome/GO-BP, and
(b) over-representation (ORA) of the significant proteins.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CAT, CELL_META, DATA, INK2, SEED, TAB, bh, read_gmt, savefig, style
from stats import gsea_preranked, moderated_t

style()
out = {}
L = pd.read_csv(TAB / "protein_matrix_cellline_mean_log2.csv.gz", index_col=0)
gene = L.pop("gene").fillna("")
sym = gene.str.split(";").str[0]
lines = list(L.columns)
EMT = [l for l in lines if CELL_META.paper_subtype[l] == "EMT"]
MET = [l for l in lines if CELL_META.paper_subtype[l] == "Metabolism"]

# --- differential proteins ------------------------------------------------
d = moderated_t(L, EMT, MET, min_per_group=3)   # logFC = EMT - Metabolism
d.insert(0, "gene", gene.reindex(d.index))
d["sym"] = sym.reindex(d.index)
d = d.sort_values("p")
d.to_csv(TAB / "diff_EMT_vs_Metabolism.csv.gz", float_format="%.5g")
out["diff"] = {"proteins_tested": int(len(d)), "q<0.05": int((d.q < 0.05).sum()),
               "q<0.05_up_EMT": int(((d.q < 0.05) & (d.logFC > 0)).sum()),
               "q<0.05_up_Metab": int(((d.q < 0.05) & (d.logFC < 0)).sum())}

# gene-level rank statistic = mean moderated t per leading gene symbol
rank = d.query("sym != ''").groupby("sym").t.mean()

# --- preranked GSEA across collections ------------------------------------
universe = set(rank.index)
collections = {"H": "h.all.v2023.1.Hs.symbols.gmt", "KEGG": "c2.cp.kegg.v7.5.1.symbols.gmt",
               "REACTOME": "c2.cp.reactome.v7.5.1.symbols.gmt", "GOBP": "c5.go.bp.v7.5.1.symbols.gmt"}
all_sets = {}
for tag, f in collections.items():
    for k, v in read_gmt(DATA / "genesets" / f, universe, 12, 300).items():
        all_sets[f"{tag}:{k}"] = v
g = gsea_preranked(rank, all_sets, min_size=12, max_size=300, n_perm=2000, seed=SEED)
g.to_csv(TAB / "gsea_EMT_vs_Metabolism.csv", index=False, float_format="%.4g")
sigpos = g[(g.q < 0.1) & (g.NES > 0)]
signeg = g[(g.q < 0.1) & (g.NES < 0)]
out["gsea"] = {"sets_tested": int(len(g)), "q<0.1_up_EMT": int(len(sigpos)), "q<0.1_up_Metab": int(len(signeg)),
               "top_EMT": sigpos.head(12).set.tolist(), "top_Metab": signeg.head(12).set.tolist()}
print(json.dumps(out, indent=2, default=float))

# --- figure: top enriched pathways (barbell NES) --------------------------
def clean(s):
    return s.split(":", 1)[1].replace("HALLMARK_", "").replace("REACTOME_", "R:").replace("KEGG_", "K:") \
            .replace("GOBP_", "GO:")[:46]


topN = pd.concat([signeg.sort_values("NES").head(14), sigpos.sort_values("NES", ascending=False).head(14)])
topN = topN.sort_values("NES")
fig, ax = plt.subplots(figsize=(8.8, 9))
colors = [CAT[1] if v > 0 else CAT[0] for v in topN.NES]
ax.barh(range(len(topN)), topN.NES, color=colors, height=0.72)
ax.set_yticks(range(len(topN)), [clean(s) for s in topN.set], fontsize=7)
ax.axvline(0, color="#c3c2b7", lw=0.6)
ax.set_xlabel("GSEA NES  (positive = higher in EMT, negative = higher in Metabolism)")
ax.set_title("Pathways enriched between GC proteome subtypes\n(preranked GSEA, q < 0.1)", fontsize=11)
for i, (nes, q) in enumerate(zip(topN.NES, topN.q)):
    ax.text(nes + (0.03 if nes > 0 else -0.03), i, f"q={q:.2g}", va="center",
            ha="left" if nes > 0 else "right", fontsize=6, color=INK2)
ax.set_xlim(topN.NES.min() - 0.9, topN.NES.max() + 0.9)
savefig(fig, "03_gsea_subtype")
plt.close(fig)

# --- figure: leading-edge enrichment plot for two flagship sets -----------
def running_es(rank, members, weight=1.0):
    r = rank.sort_values(ascending=False)
    genes = r.index.to_numpy(); sv = np.abs(r.to_numpy()) ** weight
    hit = np.isin(genes, list(members))
    Nhit = sv[hit].sum(); Nmiss = (~hit).sum()
    run = np.cumsum(np.where(hit, sv, 0)) / Nhit - np.cumsum(np.where(hit, 0, 1)) / Nmiss
    return run, hit, genes


flag = {}
for key in ["H:HALLMARK_OXIDATIVE_PHOSPHORYLATION", "H:HALLMARK_MYC_TARGETS_V1"]:
    if key in all_sets:
        flag[key] = all_sets[key]
fig, axes = plt.subplots(len(flag), 1, figsize=(8, 2.6 * len(flag)), gridspec_kw={"hspace": 0.55})
for ax, (key, mem) in zip(np.atleast_1d(axes), flag.items()):
    run, hit, genes = running_es(rank, mem)
    ax.plot(run, color=CAT[0] if run[np.argmax(np.abs(run))] < 0 else CAT[1], lw=1.5)
    ax.axhline(0, color="#c3c2b7", lw=0.5)
    for i in np.where(hit)[0]:
        ax.vlines(i, min(run) - 0.05, min(run) - 0.01, color=INK2, lw=0.3)
    row = g[g.set == key].iloc[0]
    ax.set_title(f"{clean(key)}  (NES {row.NES:.2f}, q {row.q:.2g})", fontsize=9)
    ax.set_ylabel("running ES")
axes[-1].set_xlabel("proteins ranked by EMT−Metabolism moderated t  (left = up in EMT)")
savefig(fig, "04_gsea_leadingedge")
plt.close(fig)

json.dump(out, open(TAB / "enrichment_summary.json", "w"), indent=2, default=float)
