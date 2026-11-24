#!/usr/bin/env python3
"""Step 04a - 15-aa flanking windows for PP / YP sites (PTMsigDB v1.9.1 'flanking' IDs).

PP: window built from UniProt (Swiss-Prot incl. isoforms, downloaded 2026) using
    UniprotID + residue position from Table S3. Sites whose residue does not match the
    sequence are flagged (window = NA).
YP: window = central 15 aa of MaxQuant's 31-aa 'Sequence window'.
Output: data_cache/processed/{PP,YP}_windows.tsv (id, window15)
"""
import os
import re
import pandas as pd

RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DC = os.path.join(RUN, "data_cache")
P = os.path.join(DC, "processed")

seqs, cur = {}, None
for line in open(os.path.join(DC, "uniprot", "human_sp_isoforms.fasta")):
    if line.startswith(">"):
        cur = line.split("|")[1]
        seqs[cur] = []
    else:
        seqs[cur].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}


def window(acc, res):
    m = re.match(r"([STY])(\d+)$", str(res))
    if not m:
        return None
    aa, pos = m.group(1), int(m.group(2))
    for a in [acc, str(acc).split("-")[0]]:
        s = seqs.get(a)
        if s and pos <= len(s) and s[pos - 1] == aa:
            pad = "_" * 7 + s + "_" * 7
            return pad[pos - 1: pos + 14]
    return None


pp = pd.read_csv(os.path.join(P, "PP_S3.tsv"), sep="\t", usecols=["id", "UniprotID", "Residue"])
pp["window15"] = [window(a, r) for a, r in zip(pp.UniprotID, pp.Residue)]
print("PP windows found:", pp.window15.notna().sum(), "/", len(pp))
pp[["id", "window15"]].to_csv(os.path.join(P, "PP_windows.tsv"), sep="\t", index=False)

yp = pd.read_csv(os.path.join(P, "YP_S3.tsv"), sep="\t", usecols=["id", "SeqWindow"])
yp["window15"] = yp.SeqWindow.str[8:23].str.replace("_", "_")
print("YP windows:", yp.window15.str.len().eq(15).sum(), "/", len(yp))
yp[["id", "window15"]].to_csv(os.path.join(P, "YP_windows.tsv"), sep="\t", index=False)
