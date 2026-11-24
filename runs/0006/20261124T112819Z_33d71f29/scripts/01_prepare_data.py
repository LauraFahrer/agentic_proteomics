#!/usr/bin/env python3
"""Step 01 - build analysis matrices.

Two data sources are prepared:

A. Author-processed matrices from the publication (Table S3, MOL2-18-1437-s005.xls):
   WP (proteins x 42), PP (phosphosites x 42), YP (pTyr sites x 28). These were produced
   by the authors with Spectronaut 13 / MaxQuant 1.6.0.16, quantile normalisation,
   imputation (Spectronaut / PhosR / Perseus) and proBatch. They are the closest
   available input to the paper's downstream analyses (exact upstream reproduction is
   impossible: Spectronaut is commercial and PP/YP raw data were not reprocessed by quantms).

B. quantms reprocessing e60b9e9aa4c8 (DIA-NN 2.5.1, MaxLFQ, WP only). We apply the
   paper's *data-processing* steps (log2 -> quantile normalisation -> imputation) so the
   downstream WP subtyping can be repeated on an independently processed matrix.

Outputs (data_cache/processed/*.tsv) are consumed by the R scripts.
"""
import os
import re
import sys

import numpy as np
import pandas as pd

RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DC = os.path.join(RUN, "data_cache")
SUPP = os.path.join(DC, "supp", "x", "MOL2-18-1437-s005.xls")
OUT = os.path.join(DC, "processed")
TAB = os.path.join(RUN, "results", "tables")
os.makedirs(OUT, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

LINES = ["AGS", "AZ521", "eAGS", "KatoIII", "MKN28", "N87", "NCC24", "Ocum1",
         "SNU1750", "SNU601", "SNU719", "TMK1", "YCC3", "YCCEL1"]
FIX = {"M28": "MKN28", "Ocum": "Ocum1", "TMK": "TMK1", "YCCEL": "YCCEL1", "YCCE1": "YCCEL1"}


def std_col(c):
    m = re.match(r"^(.+?)_[rR](\d)$", str(c))
    if not m:
        return None
    line = FIX.get(m.group(1), m.group(1))
    assert line in LINES, c
    return f"{line}_r{m.group(2)}"


def order_cols(cols):
    return sorted(cols, key=lambda c: (LINES.index(c.split("_r")[0]), int(c.split("_r")[1])))


# ---------------------------------------------------------------- A. Table S3
xl = pd.read_excel(SUPP, sheet_name=None)
wp = xl["processed WP matrix"]
pp = xl["processed PP matrix"]
yp = xl["processed YP matrix"]

# --- WP
wp = wp.rename(columns={c: std_col(c) for c in wp.columns if std_col(c)})
samp = order_cols([c for c in wp.columns if re.match(r".+_r\d$", str(c))])
wp["Accession"] = wp["Protein Accessions"].astype(str).str.split(";").str[0]
# Table S3 WP contains the iRT spike-in (Biognosys HRM peptides) and porcine trypsin -> removed
n_wp_nonhuman = int(wp["Accession"].isin(["iRT", "TRYP_PIG"]).sum())
wp = wp[~wp["Accession"].isin(["iRT", "TRYP_PIG"])].copy()
wp["Gene"] = wp["Genes"].fillna("").astype(str).str.split(";").str[0]
wp.loc[wp["Gene"] == "", "Gene"] = wp.loc[wp["Gene"] == "", "Accession"]
wp["id"] = wp["Gene"]
# duplicated first gene symbols -> make unique as GENE|ACCESSION (documented)
wp["_mean"] = wp[samp].mean(axis=1)
wp = wp.sort_values("_mean", ascending=False)
dup = wp["id"].duplicated(keep=False)
wp.loc[dup, "id"] = wp.loc[dup, "Gene"] + "|" + wp.loc[dup, "Accession"]
wp = wp.sort_values("id")
wp_out = wp[["id", "Gene", "Accession", "Protein Accessions", "Genes"] + samp]
wp_out.to_csv(os.path.join(OUT, "WP_S3.tsv"), sep="\t", index=False)

# --- PP
pp = pp.rename(columns={c: std_col(c) for c in pp.columns if std_col(c)})
psamp = order_cols([c for c in pp.columns if re.match(r".+_r\d$", c)])
pp["id"] = pp["GeneSymbol"].astype(str) + "_" + pp["Residue"].astype(str)
# a site id can appear on several rows (isoforms); keep first, report count
n_dup_pp = int(pp["id"].duplicated().sum())
pp = pp.drop_duplicates("id")
pp_out = pp[["id", "GeneSymbol", "isoforms", "Residue", "UniprotID", "Sequence"] + psamp]
pp_out.to_csv(os.path.join(OUT, "PP_S3.tsv"), sep="\t", index=False)

# --- YP (MaxQuant Phospho(STY)Sites-derived)
yp = yp.rename(columns={c: std_col(c) for c in yp.columns if std_col(c)})
ysamp = order_cols([c for c in yp.columns if re.match(r".+_r\d$", c)])


def first_acc(s):
    s = str(s).split(";")[0]
    m = re.match(r"(?:sp|tr)\|([^|]+)\|", s)
    return m.group(1) if m else s


def first_name(s):
    s = str(s).split(";")[0]
    m = re.match(r"(?:sp|tr)\|[^|]+\|(\S+?)_HUMAN", s)
    return m.group(1) if m else ""


yp["LeadAcc"] = yp["Leading proteins"].map(first_acc)
lead_list = yp["Leading proteins"].astype(str).str.split(";")
prot_list = yp["Proteins"].astype(str).str.split(";")
pos_list = yp["Positions within proteins"].astype(str).str.split(";")


def lead_pos(i):
    lead = lead_list.iloc[i][0]
    try:
        j = prot_list.iloc[i].index(lead)
        return pos_list.iloc[i][j]
    except ValueError:
        return pos_list.iloc[i][0]


yp["LeadPos"] = [lead_pos(i) for i in range(len(yp))]
yp["SeqWindow"] = yp["Sequence window"].astype(str).str.split(";").str[0]
yp["Residue_aa"] = yp["SeqWindow"].str[15]
yp["UniProtName"] = yp["Leading proteins"].map(first_name)
yp.to_csv(os.path.join(OUT, "YP_S3_raw.tsv"), sep="\t", index=False)

# --------------------------------------------------- B. quantms WP reprocessing
pg = pd.read_parquet(os.path.join(DC, "quantms", "PXD048538.pg.parquet"))
pg = pg[(~pg["is_decoy"]) & (~pg["contaminant"])]
pg = pg[(pg["pg_qvalue"] <= 0.01) & (pg["global_qvalue"] <= 0.01)]
pg["run"] = pg["grouped_runs"].map(lambda x: x[0])
pg["sample"] = pg["run"].str.replace("Gastric_WP_DIA_", "", regex=False).map(std_col)
pg["gene"] = pg["gg_names"].map(lambda x: ";".join(x) if x is not None else "")
qm = pg.pivot_table(index=["anchor_protein", "gene"], columns="sample",
                    values="intensity", aggfunc="max")
qm = qm[order_cols(list(qm.columns))]
n_ids = qm.notna().sum()
qm = np.log2(qm.replace(0, np.nan))
qm.to_csv(os.path.join(OUT, "WP_quantms_log2_raw.tsv"), sep="\t")

# protein FDR-filtered matrix; keep proteins quantified in >=50% of samples
# (the paper does not state a completeness filter: Spectronaut's default q-value
#  filtering + "imputation by Spectronaut"; we document this choice)
keep = qm.notna().mean(axis=1) >= 0.5
qk = qm[keep].astype("float64").copy()

# quantile normalisation (missing-value aware, preprocessCore-style mean of sorted columns)
def quantile_normalize(df):
    arr = df.values.copy()
    srt = np.sort(arr, axis=0)  # NaNs last
    n_obs = (~np.isnan(arr)).sum(axis=0)
    # reference: mean of rank-interpolated columns
    L = arr.shape[0]
    grid = np.linspace(0, 1, L)
    ref = np.zeros(L)
    for j in range(arr.shape[1]):
        v = srt[: n_obs[j], j]
        ref += np.interp(grid, np.linspace(0, 1, len(v)), v)
    ref /= arr.shape[1]
    out = np.full_like(arr, np.nan)
    for j in range(arr.shape[1]):
        ok = ~np.isnan(arr[:, j])
        ranks = pd.Series(arr[ok, j]).rank(method="average").values
        q = (ranks - 1) / max(ok.sum() - 1, 1)
        out[ok, j] = np.interp(q, grid, ref)
    return pd.DataFrame(out, index=df.index, columns=df.columns)


qn = quantile_normalize(qk)
# imputation: Spectronaut's imputation is proprietary. Approximation: per-sample
# down-shifted normal (Perseus-like; width 0.3, shift 1.8 SD), seeded.
rng = np.random.default_rng(20240101)
qi = qn.copy()
for c in qi.columns:
    v = qi[c]
    mu, sd = v.mean(), v.std()
    na = v.isna()
    qi.loc[na, c] = rng.normal(mu - 1.8 * sd, 0.3 * sd, na.sum())
qi = qi.reset_index()
qi["Gene"] = qi["gene"].str.split(";").str[0]
qi.loc[qi["Gene"] == "", "Gene"] = qi["anchor_protein"]
qi["_mean"] = qi[order_cols([c for c in qi.columns if "_r" in c])].mean(axis=1)
qi = qi.sort_values("_mean", ascending=False)
qi["id"] = qi["Gene"]
d = qi["id"].duplicated(keep=False)
qi.loc[d, "id"] = qi.loc[d, "Gene"] + "|" + qi.loc[d, "anchor_protein"]
cols = order_cols([c for c in qi.columns if re.match(r".+_r\d$", str(c))])
qi = qi.sort_values("id")
qi[["id", "Gene", "anchor_protein"] + cols].to_csv(os.path.join(OUT, "WP_quantms_processed.tsv"), sep="\t", index=False)

summary = {
    "S3_WP_proteins": int(len(wp_out)), "S3_WP_nonhuman_rows_removed": n_wp_nonhuman, "S3_PP_sites": int(len(pp_out)),
    "S3_PP_duplicated_site_ids_dropped": n_dup_pp, "S3_YP_sites": int(len(yp)),
    "quantms_proteins_any_sample_1pctFDR": int(len(qm)),
    "quantms_proteins_ge50pct_samples": int(keep.sum()),
    "quantms_ids_per_sample_median": float(n_ids.median()),
    "quantms_missing_fraction_after_filter": float(qk.isna().mean().mean()),
}
pd.Series(summary).to_json(os.path.join(TAB, "01_data_summary.json"), indent=1)
print(pd.Series(summary))
