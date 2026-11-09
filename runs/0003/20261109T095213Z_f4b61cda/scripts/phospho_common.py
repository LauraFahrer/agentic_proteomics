"""Loaders for the (phospho)proteome matrices of PXD048538.

IMPORTANT: the quantms dataset version e60b9e9aa4c8 contains ONLY the whole proteome.
The phosphoproteome (PP, Ser/Thr/Tyr) and tyrosine-phosphorylation (YP) layers used here
come from the associated publication's Supplementary Table S3 (file MOL2-18-1437-s005.xls,
sheets 'processed PP/YP/WP matrix'): author-processed, log2, quantile-normalised, imputed and
batch-corrected. This is stated throughout the report.
"""
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(os.environ.get("PXD048538_DATA", "data_cache"))
SUPP = DATA / "supp/x/MOL2-18-1437-s005.xls"
FASTA = DATA.parent / "ebv/db/human_sp.fasta"   # cached UniProt human Swiss-Prot (run 0001)

# phospho sample-name -> canonical cell-line name (common.CELL_META index)
LINE_FIX = {"M28": "MKN28", "Ocum": "Ocum1", "TMK": "TMK1"}


def _fix(name):
    return LINE_FIX.get(name, name)


def load_fasta(path=FASTA):
    seqs, acc, cur = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if acc:
                seqs[acc] = "".join(cur)
            acc = line.split("|")[1] if "|" in line else line[1:].split()[0]
            cur = []
        else:
            cur.append(line.strip())
    if acc:
        seqs[acc] = "".join(cur)
    return seqs


def acc2gene(path=FASTA):
    m = {}
    for line in open(path):
        if line.startswith(">") and "|" in line:
            a = line.split("|")[1]
            g = re.search(r"GN=([A-Za-z0-9\-]+)", line)
            if g:
                m[a] = g.group(1)
    return m


def window10(seq, pos):
    """10-mer with the residue at 1-based `pos` placed at position 6 (index 5); '_' padded."""
    if seq is None or pos < 1 or pos > len(seq):
        return None
    out = []
    for i in range(pos - 6, pos + 4):          # protein indices (0-based) pos-6 .. pos+3
        out.append(seq[i] if 0 <= i < len(seq) else "_")
    return "".join(out)


def load_PP(line_order):
    df = pd.read_excel(SUPP, sheet_name="processed PP matrix")
    meta = df[["GeneSymbol", "isoforms", "Residue", "UniprotID", "Sequence"]].copy()
    for c in ["GeneSymbol", "isoforms", "Residue", "UniprotID", "Sequence"]:
        meta[c] = meta[c].astype(str)   # fix Excel gene->date corruption (SEPT/MARCH)
    samp = [c for c in df.columns if re.match(r".+_r\d$", c)]
    X = df[samp].copy()
    cols = pd.MultiIndex.from_tuples([(_fix(c.rsplit("_r", 1)[0]), "r" + c.rsplit("_r", 1)[1]) for c in samp])
    X.columns = cols
    lines = [l for l in line_order if l in X.columns.get_level_values(0)]
    linemean = pd.DataFrame({l: X[l].mean(axis=1) for l in lines})
    meta["restype"] = meta.Residue.str[0]
    meta["pos"] = pd.to_numeric(meta.Residue.str[1:], errors="coerce")
    meta["site"] = meta.GeneSymbol + "_" + meta.Residue
    fasta = load_fasta()
    meta["win10"] = [window10(fasta.get(u), int(p)) if pd.notna(p) else None
                     for u, p in zip(meta.UniprotID, meta.pos)]
    return meta.reset_index(drop=True), linemean.reset_index(drop=True)


def load_YP(line_order):
    df = pd.read_excel(SUPP, sheet_name="processed YP matrix")
    samp = [c for c in df.columns if re.match(r".+_R\d$", c)]
    X = df[samp].copy()
    cols = pd.MultiIndex.from_tuples([(_fix(c.rsplit("_R", 1)[0]), "R" + c.rsplit("_R", 1)[1]) for c in samp])
    X.columns = cols
    lines = [l for l in line_order if l in X.columns.get_level_values(0)]
    linemean = pd.DataFrame({l: X[l].mean(axis=1) for l in lines})
    win = df["Sequence window"].astype(str).str.split(";").str[0]
    # central residue of the (odd-length) window is the phosphoacceptor; take the central 10-mer
    def centre10(w):
        if not w or w == "nan":
            return None
        c = len(w) // 2
        s = w[c - 5:c + 5]
        return s if len(s) == 10 else None
    lead = df["Leading proteins"].astype(str).str.extract(r"\|([A-Z0-9\-]+)\|")[0]
    gene = df["Leading proteins"].astype(str).str.extract(r"GN=([A-Za-z0-9\-]+?)(?:PE=|OX=|SV=)")[0]
    pos = df["Positions within proteins"].astype(str).str.split(";").str[0]
    a2g = acc2gene()
    gene = gene.fillna(lead.str.replace(r"-\d+$", "", regex=True).map(a2g))
    meta = pd.DataFrame({"UniprotID": lead, "GeneSymbol": gene, "pos": pd.to_numeric(pos, errors="coerce"),
                         "win10": [centre10(w) for w in win], "restype": "Y"})
    meta["site"] = meta.GeneSymbol.fillna(meta.UniprotID) + "_Y" + meta.pos.astype("Int64").astype(str)
    return meta.reset_index(drop=True), linemean.reset_index(drop=True)


def load_WP_linemean(line_order):
    """Author WP matrix (same workbook) as cell-line means, keyed by leading UniProt accession + gene."""
    df = pd.read_excel(SUPP, sheet_name="processed WP matrix")
    acc = df["Protein Accessions"].astype(str).str.split(";").str[0]
    gene = df["Genes"].astype(str).str.split(";").str[0]
    samp = [c for c in df.columns if re.match(r".+_R\d$", c)]
    X = df[samp].copy()
    X.columns = pd.MultiIndex.from_tuples([(_fix(c.rsplit("_R", 1)[0]), c) for c in samp])
    lines = [l for l in line_order if l in X.columns.get_level_values(0)]
    lm = pd.DataFrame({l: X[l].mean(axis=1) for l in lines})
    lm.index = acc.values
    byacc = lm.groupby(level=0).mean()
    lm2 = lm.copy(); lm2.index = gene.values
    bygene = lm2.groupby(level=0).mean()
    return byacc, bygene


def write_seqrnk(win_scores, path):
    """win_scores: iterable of (win10, score) with valid 10-mers; writes a phosX seqrnk."""
    rows = [(w, s) for w, s in win_scores if isinstance(w, str) and len(w) == 10 and "_" * 10 != w]
    df = pd.DataFrame(rows, columns=["seq", "score"]).dropna()
    df = df.sort_values("score", ascending=False)
    df.to_csv(path, sep="\t", header=False, index=False, float_format="%.5f")
    return len(df)
