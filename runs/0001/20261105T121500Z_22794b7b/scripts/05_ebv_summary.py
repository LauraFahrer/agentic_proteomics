"""Summarise the human+EBV Sage search (05_ebv_search.sh).

EBV evidence = target PSMs at 1% spectrum- and peptide-level q-value whose peptide maps
ONLY to EBV proteins (not to any human protein, I/L treated as equal). Because EBV is a tiny
part of the database, the global FDR does not guarantee a low FDR within the EBV subset, so
the number of EBV-*decoy* PSMs passing the same thresholds is reported as a subset-specific
false-hit estimate.
Usage: python 05_ebv_summary.py <work_dir>
"""
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import CAT, FIG, TAB, savefig, style

style()
W = Path(sys.argv[1])
psm = pd.read_csv(W / "sage/out/results.sage.tsv", sep="\t")
human = "".join(l.strip() for l in open(W / "db/human_sp.fasta") if not l.startswith(">")).replace("I", "L")
ebv_acc = {l.split("|")[1] for l in open(W / "db/ebv.fasta") if l.startswith(">")}
ebv_name = {l.split("|")[1]: l.split("|")[2].split(" OS=")[0] for l in open(W / "db/ebv.fasta") if l.startswith(">")}


def accs(p):
    return [x.split("|")[1] if "|" in x else x for x in str(p).split(";")]


psm["accs"] = psm.proteins.map(lambda p: [a.replace("rev_", "") for a in accs(p)])
psm["all_ebv"] = psm.accs.map(lambda a: all(x in ebv_acc for x in a))
psm["stripped"] = psm.peptide.map(lambda s: re.sub(r"\[.*?\]|[^A-Z]", "", s))
psm["in_human"] = psm.stripped.map(lambda s: s.replace("I", "L") in human)
psm["line"] = psm.filename.str.extract(r"DDAlibrary_(\w+)\.mzML")[0]
ok = (psm.spectrum_q <= 0.01) & (psm.peptide_q <= 0.01)
tgt = psm[ok & (psm.label == 1)]
dec = psm[ok & (psm.label == -1)]
summ = {"psms_target_1pct": int(len(tgt)),
        "per_line_target_psms": tgt.groupby("line").size().to_dict(),
        "per_line_target_peptides": tgt.groupby("line").stripped.nunique().to_dict()}
ebv_t = tgt[tgt.all_ebv & ~tgt.in_human]
ebv_d = dec[dec.all_ebv]
summ["EBV_target_psms_per_line"] = ebv_t.groupby("line").size().to_dict()
summ["EBV_target_peptides_per_line"] = ebv_t.groupby("line").stripped.nunique().to_dict()
summ["EBV_decoy_psms_per_line"] = ebv_d.groupby("line").size().to_dict()
rows = ebv_t.assign(protein=ebv_t.accs.map(lambda a: "/".join(sorted({ebv_name.get(x, x) for x in a})))) \
    [["line", "stripped", "peptide", "protein", "hyperscore", "spectrum_q", "peptide_q", "protein_q", "charge"]] \
    .sort_values(["line", "protein", "hyperscore"], ascending=[True, True, False])
rows.to_csv(TAB / "ebv_peptide_hits.csv", index=False)
summ["EBV_proteins_per_line"] = rows.groupby("line").protein.apply(lambda s: sorted(set(s))).to_dict()
json.dump(summ, open(TAB / "ebv_search_summary.json", "w"), indent=2)
print(json.dumps(summ, indent=2))

lines = sorted(tgt.line.unique())
order = [l for l in ["AGS", "eAGS", "TMK1", "SNU719", "YCCEL1", "NCC24", "N87"] if l in lines]
fig, ax = plt.subplots(figsize=(6, 3.2))
vals = [summ["EBV_target_peptides_per_line"].get(l, 0) for l in order]
dvals = [summ["EBV_decoy_psms_per_line"].get(l, 0) for l in order]
x = range(len(order))
ax.bar([i - 0.2 for i in x], vals, width=0.38, color=CAT[6], label="EBV-unique target peptides")
ax.bar([i + 0.2 for i in x], dvals, width=0.38, color="#c3c2b7", label="EBV decoy PSMs (false-hit estimate)")
for i, v in zip(x, vals):
    ax.text(i - 0.2, v, str(v), ha="center", va="bottom", fontsize=8)
ax.set_xticks(list(x), order)
ax.set_ylabel("count (1% PSM & peptide FDR)")
ax.set_title("EBV peptides in whole-proteome DDA runs (human + EBV search)")
ax.legend(fontsize=7)
savefig(fig, "15_ebv_peptide_search")
plt.close(fig)
