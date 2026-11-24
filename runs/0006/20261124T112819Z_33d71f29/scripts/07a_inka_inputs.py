#!/usr/bin/env python3
"""Step 07a - build MaxQuant-like INKA inputs from the phosphosite-centric YP matrix (Table S3).

The paper uploaded the YP site matrix to the (now unreachable) InKA web server. We run the
original InKA R code (v1.6.0 / InKAscore 1.5.0, intensity mode) locally, which expects
MaxQuant 'Phospho (STY)Sites.txt' + a 'pp-Peptide-report.txt'. We therefore create:
  * one pseudo phosphopeptide per pY site = fully tryptic peptide (cleave C-term to K/R,
    not before P) containing the site, from UniProt; fallback = central 15-aa window.
  * intensity columns = 2^log2 (Table S3 values are log2, normalised, imputed).
Output: data_cache/inka/work/{pp-Peptide-report.txt, Phospho (STY)Sites.txt}
"""
import os
import re
import numpy as np
import pandas as pd

RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DC = os.path.join(RUN, "data_cache")
W = os.path.join(DC, "inka", "work")
os.makedirs(W, exist_ok=True)

seqs, cur = {}, None
for line in open(os.path.join(DC, "uniprot", "human_sp_isoforms.fasta")):
    if line.startswith(">"):
        cur = line.split("|")[1]; seqs[cur] = []
    else:
        seqs[cur].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}


def tryptic_peptide(seq, pos):  # pos 1-based
    i = pos - 1
    s = i
    while s > 0 and not (seq[s - 1] in "KR" and seq[s] != "P"):
        s -= 1
    e = i
    while e < len(seq) - 1 and not (seq[e] in "KR" and seq[e + 1] != "P"):
        e += 1
    return seq[s:e + 1]


raw = pd.read_csv(os.path.join(DC, "processed", "YP_S3_raw.tsv"), sep="\t")
ann = pd.read_csv(os.path.join(DC, "processed", "YP_S3.tsv"), sep="\t", usecols=["id", "GeneSymbol"])
samples = [c for c in raw.columns if re.match(r".+_r\d$", c)]


def accs(s):
    return [m.group(1) for m in re.finditer(r"(?:sp|tr)\|([^|]+)\|", str(s))]


rows_pp, rows_ps = [], []
n_tryp = 0
for i, r in raw.iterrows():
    acc_list = accs(r["Proteins"])
    pos_list = str(r["Positions within proteins"]).split(";")
    # de-duplicate accession/position pairs (headers are repeated across databases)
    seen, A, P = set(), [], []
    for a, p in zip(acc_list, pos_list):
        if (a, p) not in seen:
            seen.add((a, p)); A.append(a); P.append(p)
    lead = r["LeadAcc"]; lpos = int(r["LeadPos"])
    s = seqs.get(lead) or seqs.get(lead.split("-")[0])
    if s and lpos <= len(s) and s[lpos - 1] == "Y":
        pep = tryptic_peptide(s, lpos); n_tryp += 1
    else:
        pep = r["SeqWindow"][8:23].replace("_", "")
    pid = i + 1
    inten = {f"Intensity.{c}": 2 ** r[c] for c in samples}
    rows_pp.append({"ppModPeptide.ID": pid, "Sequence": pep, "Modifications": "Phospho (STY)", "Mass": 0,
                    "Proteins": ";".join(A), "Gene.names": ann.GeneSymbol.iloc[i],
                    "Protein.Names": "", "Phospho..STY.": 1, "MS.MS.Count": 1, "Score": 100,
                    "Phospho..STY..site.IDs": pid, **inten,
                    **{f"Norm.Intensity.{c}": 2 ** r[c] for c in samples},
                    **{f"OPL.ppCount.{c}": 2 ** r[c] for c in samples}})
    rows_ps.append({"Proteins": ";".join(A), "Positions within proteins": ";".join(P),
                    "Localization prob": 1.0, "Amino acid": "Y", "Sequence window": r["SeqWindow"],
                    "Reverse": "", "Potential contaminant": "", "id": pid, "Mod. peptide IDs": pid})

pp = pd.DataFrame(rows_pp)
pp.to_csv(os.path.join(W, "pp-Peptide-report.txt"), sep="\t", index=False)
# one dummy reverse/contaminant row (as in any MaxQuant output) so that fread types the
# Reverse / Potential contaminant columns as character; it is removed by InKA's own filter
rows_ps.append({"Proteins": "REV__P00000", "Positions within proteins": "1", "Localization prob": 1.0,
                "Amino acid": "Y", "Sequence window": "_" * 31, "Reverse": "+", "Potential contaminant": "+",
                "id": 999999, "Mod. peptide IDs": 999999})
pd.DataFrame(rows_ps).to_csv(os.path.join(W, "Phospho (STY)Sites.txt"), sep="\t", index=False)
print("sites:", len(pp), "tryptic peptides derived:", n_tryp)
