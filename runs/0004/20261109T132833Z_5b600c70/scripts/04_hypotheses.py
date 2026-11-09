"""Validate strong co-abundance pairs against STRING (to separate 'known but missing from the
GO-CC reference' from genuinely novel), test complex-level subtype association, and assemble
the three PPI hypotheses.
"""
import io
import json
import time
import urllib.parse
import urllib.request

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from common import CAT, CELL_META, DATA, INK2, SUBTYPE_COLORS, TAB, bh, diverging_cmap, read_gmt, savefig, style

style()
out = {}
G = pd.read_csv(TAB / "gene_matrix_cellline_mean_log2.csv.gz", index_col=0)
C = pd.read_parquet(TAB / "coabundance_bicor.parquet")
strong = pd.read_csv(TAB / "strong_coabundance_pairs.csv.gz")
lines = list(G.columns)
sub = CELL_META.paper_subtype.reindex(lines)
genes = sorted(set(strong.geneA) | set(strong.geneB))


def string_id_map(ids):
    """Map input gene symbols to STRING preferred names (handles aliases, e.g. PHB1->PHB)."""
    m = {}
    for i in range(0, len(ids), 500):
        chunk = ids[i:i + 500]
        data = urllib.parse.urlencode({"identifiers": "\r".join(chunk), "species": 9606,
                                       "caller_identity": "agentic_proteomics", "limit": 1}).encode()
        try:
            req = urllib.request.Request("https://string-db.org/api/tsv/get_string_ids", data=data)
            df = pd.read_csv(io.StringIO(urllib.request.urlopen(req, timeout=90).read().decode()), sep="\t")
            for qi, p in zip(df.queryIndex, df.preferredName):   # queryIndex is 0-based into chunk
                m[chunk[int(qi)]] = p
            time.sleep(1)
        except Exception as e:
            out.setdefault("string_errors", []).append("idmap:" + str(e))
    return m


def string_network(ids, network_type="functional"):
    """Query STRING for ALL edges among `ids` in one call (network API returns only edges among
    the submitted identifiers, so chunking would miss cross-chunk pairs)."""
    data = urllib.parse.urlencode({"identifiers": "\r".join(ids), "species": 9606,
                                   "caller_identity": "agentic_proteomics",
                                   "network_type": network_type}).encode()
    try:
        req = urllib.request.Request("https://string-db.org/api/tsv/network", data=data)
        txt = urllib.request.urlopen(req, timeout=150).read().decode()
        time.sleep(1)
        return pd.read_csv(io.StringIO(txt), sep="\t")
    except Exception as e:
        out.setdefault("string_errors", []).append(str(e))
        return pd.DataFrame(columns=["preferredName_A", "preferredName_B", "score", "escore"])


# STRING functional (combined) and physical edges among the strong-pair genes
sn = string_network(genes, "functional")
sp = string_network(genes, "physical")


def edge_map(df, col="score"):
    m = {}
    for a, b, s in zip(df.preferredName_A, df.preferredName_B, df[col]):
        m[frozenset((a, b))] = max(m.get(frozenset((a, b)), 0), float(s))
    return m


comb = edge_map(sn, "score"); phys = edge_map(sp, "score"); expc = edge_map(sn, "escore")
idm = string_id_map(genes)   # my symbol -> STRING preferred name
def gp(a, b, m):
    return m.get(frozenset((idm.get(a, a), idm.get(b, b))), 0.0)
strong["string_combined"] = [gp(a, b, comb) for a, b in zip(strong.geneA, strong.geneB)]
strong["string_physical"] = [gp(a, b, phys) for a, b in zip(strong.geneA, strong.geneB)]
strong["string_experimental"] = [gp(a, b, expc) for a, b in zip(strong.geneA, strong.geneB)]
strong["support"] = np.where(strong.in_known_complex, "GO-CC complex",
                     np.where(strong.string_combined >= 0.4, "STRING-supported",
                     np.where(strong.string_combined < 0.15, "candidate novel", "weak")))
strong.sort_values("bicor", ascending=False).to_csv(TAB / "strong_pairs_string_annotated.csv.gz", index=False)
vc = strong.support.value_counts().to_dict()
out["strong_pairs_total"] = int(len(strong))
out["support_breakdown"] = vc
out["pct_known_any"] = float((strong.in_known_complex | (strong.string_combined >= 0.4)).mean() * 100)
novel = strong[strong.support == "candidate novel"].sort_values("bicor", ascending=False)
out["candidate_novel_pairs"] = int(len(novel))
out["top_candidate_novel"] = novel.head(20)[["geneA", "geneB", "bicor", "string_combined"]].values.tolist()

# are the candidate-novel-pair genes enriched (vs all strong-pair genes) for a cell compartment?
from scipy.stats import fisher_exact
bg = sorted(set(strong.geneA) | set(strong.geneB))
novg = set(novel.geneA) | set(novel.geneB)
gocc_all = read_gmt(DATA / "genesets/c5.go.cc.v7.5.1.symbols.gmt", set(bg), 3, 2000)
erows = []
for k, v in gocc_all.items():
    a = len(novg & v); c = len((set(bg) - novg) & v)
    if a >= 5:
        orr, p = fisher_exact([[a, len(novg) - a], [c, len(bg) - len(novg) - c]], alternative="greater")
        erows.append((k.replace("GOCC_", ""), a, len(v & set(bg)), float(orr), float(p)))
enr = pd.DataFrame(erows, columns=["term", "novel_genes", "bg_genes", "OR", "p"]).sort_values("p")
enr["q"] = bh(enr.p.values)
enr.to_csv(TAB / "candidate_novel_gene_enrichment.csv", index=False)
out["novel_gene_top_compartments"] = enr.head(6)[["term", "OR", "p", "q"]].values.tolist()
print("support:", vc, "| pct known(any):", round(out["pct_known_any"], 1))
print("top candidate-novel pairs:"); print(novel.head(15)[["geneA", "geneB", "bicor", "string_combined", "string_physical"]].to_string())

# --- complex-level subtype association: mean abundance of recovered complexes vs subtype ---
gocc = read_gmt(DATA / "genesets/c5.go.cc.v7.5.1.symbols.gmt", set(G.index), 4, 80)
complexes = {k.replace("GOCC_", ""): v for k, v in gocc.items() if "COMPLEX" in k}
Z = G.sub(G.mean(axis=1), axis=0).div(G.std(axis=1) + 1e-9, axis=0)
rows = []
for name, mem in complexes.items():
    mm = [g for g in mem if g in Z.index]
    if len(mm) < 4:
        continue
    score = Z.loc[mm].mean()
    coh = C.loc[mm, mm].values
    cohesion = np.nanmean(coh[np.triu_indices(len(mm), 1)])
    a, b = score[sub == "EMT"], score[sub == "Metabolism"]
    t, p = ttest_ind(a, b, equal_var=False)
    rows.append({"complex": name, "n": len(mm), "cohesion_bicor": float(cohesion),
                 "diff_EMT_minus_Metab": float(a.mean() - b.mean()), "p": float(p)})
cx = pd.DataFrame(rows)
cx["q"] = bh(cx.p.values)
cx = cx.sort_values("p")
cx.to_csv(TAB / "complex_subtype_association.csv", index=False)
out["complexes_tested"] = int(len(cx))
out["complexes_subtype_q<0.1"] = int((cx.q < 0.1).sum())
out["top_subtype_complexes"] = cx.head(12)[["complex", "n", "diff_EMT_minus_Metab", "q", "cohesion_bicor"]].round(3).values.tolist()
print("\ncomplexes subtype-associated q<0.1:", out["complexes_subtype_q<0.1"])
print(cx.head(12).to_string())
json.dump(out, open(TAB / "hypotheses_summary.json", "w"), indent=2, default=float)

# ---- figure ----
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), gridspec_kw={"wspace": 0.42, "width_ratios": [1, 1.2, 1.1]})
# (A) support breakdown of strong pairs
ax = axes[0]
order_s = ["GO-CC complex", "STRING-supported", "weak", "candidate novel"]
cols = {"GO-CC complex": CAT[0], "STRING-supported": CAT[2], "weak": "#c3c2b7", "candidate novel": CAT[1]}
vals = [vc.get(k, 0) for k in order_s]
ax.bar(range(len(order_s)), vals, color=[cols[k] for k in order_s])
ax.set_xticks(range(len(order_s)), order_s, rotation=25, ha="right", fontsize=7)
ax.set_ylabel("strong co-abundance pairs (bicor≥0.9)")
for i, v in enumerate(vals):
    ax.text(i, v, str(v), ha="center", va="bottom", fontsize=8)
ax.set_title(f"A  Prior support of strong pairs\n{out['pct_known_any']:.0f}% known (complex or STRING)", fontsize=9)
# (B) bicor vs STRING combined
ax = axes[1]
for lab_, c in [("GO-CC complex", CAT[0]), ("STRING-supported", CAT[2]), ("candidate novel", CAT[1])]:
    d = strong[strong.support == lab_]
    ax.scatter(d.string_combined, d.bicor, s=10, color=c, label=lab_, alpha=0.7, linewidths=0)
for _, r in novel.head(10).iterrows():
    ax.annotate(f"{r.geneA}-{r.geneB}", (r.string_combined, r.bicor), fontsize=5.5, color=INK2,
                xytext=(2, 2), textcoords="offset points")
ax.set_xlabel("STRING combined score"); ax.set_ylabel("co-abundance (bicor)")
ax.set_title("B  Strong pairs vs STRING evidence", fontsize=9); ax.legend(fontsize=6)
# (C) complexes most subtype-associated
ax = axes[2]
cc = cx.head(12).iloc[::-1]
colors = [CAT[1] if v > 0 else CAT[0] for v in cc.diff_EMT_minus_Metab]
ax.barh(range(len(cc)), cc.diff_EMT_minus_Metab, color=colors, height=0.72)
ax.set_yticks(range(len(cc)), [c[:30] for c in cc.complex], fontsize=6)
ax.axvline(0, color="#c3c2b7", lw=0.6)
for i, (v, q) in enumerate(zip(cc.diff_EMT_minus_Metab, cc.q)):
    if q < 0.1:
        ax.text(v + (0.02 if v > 0 else -0.02), i, "*", va="center", ha="left" if v > 0 else "right", fontsize=11, color=INK2)
ax.set_xlabel("complex abundance EMT − Metab")
ax.set_title("C  Complex abundance vs subtype\n(* q<0.1; + EMT, − Metabolism)", fontsize=9)
savefig(fig, "03_ppi_validation_subtype")
plt.close(fig)
json.dump(out, open(TAB / "hypotheses_summary.json", "w"), indent=2, default=float)
