"""Build the protein-group x sample matrix for PXD048538 (quantms e60b9e9aa4c8).
Same pipeline as run 0001 step 01 (kept self-contained for this run).
Output: tables/protein_matrix_log2.csv.gz, protein_matrix_cellline_mean_log2.csv.gz, sample_annotation.csv
"""
import json

import numpy as np
import pandas as pd

from common import CELL_META, QDIR, TAB, short_line

ann_s = pd.read_parquet(QDIR / "PXD048538.sample.parquet")
runs = pd.read_parquet(QDIR / "PXD048538.run.parquet")
sdrf = pd.read_csv(QDIR / "PXD048538.sdrf.tsv", sep="\t")
ann = sdrf[["source name", "characteristics[cell line]", "characteristics[disease]",
            "characteristics[biological replicate]", "comment[data file]"]].copy()
ann.columns = ["sample", "cell_line_sdrf", "disease_sdrf", "bio_rep", "raw_file"]
ann["run"] = ann.raw_file.str.replace(".raw", "", regex=False)
ann["line"] = ann.cell_line_sdrf.map(short_line)
ann = ann.join(CELL_META, on="line")
ann["sample_id"] = ann.line + "_R" + ann.bio_rep.astype(str)
order = list(CELL_META.sort_values(["paper_subtype"]).index)
ann["line"] = pd.Categorical(ann.line, categories=order, ordered=True)
ann = ann.sort_values(["line", "bio_rep"]).reset_index(drop=True)
ann.to_csv(TAB / "sample_annotation.csv", index=False)

pg = pd.read_parquet(QDIR / "PXD048538.pg.parquet",
                     columns=["pg_accessions", "gg_names", "grouped_runs", "pg_qvalue",
                              "global_qvalue", "intensity", "is_decoy", "contaminant", "peptide_counts"])
pg["run"] = pg.grouped_runs.str[0]
pg["protein"] = pg.pg_accessions.map(lambda a: ";".join(a))
pg["gene"] = pg.gg_names.map(lambda a: ";".join(a) if a is not None and len(a) else "")
keep = (~pg.is_decoy) & (~pg.contaminant) & (pg.pg_qvalue <= 0.01) & (pg.global_qvalue <= 0.01) & (pg.intensity > 0)
pgf = pg[keep]
mat = pgf.pivot_table(index="protein", columns="run", values="intensity", aggfunc="max")
gene = pgf.groupby("protein").gene.first()
mat = np.log2(mat)
mat.columns = [ann.set_index("run").sample_id[c] for c in mat.columns]
mat = mat[ann.sample_id]
mat = mat - mat.median() + mat.median().mean()  # median alignment (shift only)

pres = pd.DataFrame({line: mat[list(cols)].notna().sum(axis=1) >= 2
                     for line, cols in ann.groupby("line", observed=True).sample_id})
line_mean = pd.DataFrame({line: mat[list(cols)].mean(axis=1).where(pres[line])
                          for line, cols in ann.groupby("line", observed=True).sample_id})

out = mat.copy(); out.insert(0, "gene", gene.reindex(out.index))
out.to_csv(TAB / "protein_matrix_log2.csv.gz", float_format="%.4f")
lm = line_mean.copy(); lm.insert(0, "gene", gene.reindex(lm.index))
lm.to_csv(TAB / "protein_matrix_cellline_mean_log2.csv.gz", float_format="%.4f")
summary = {"protein_groups": int(mat.shape[0]), "runs": int(mat.shape[1]), "lines": int(line_mean.shape[1]),
           "present_in_all_lines": int(pres.all(axis=1).sum())}
json.dump(summary, open(TAB / "matrix_summary.json", "w"), indent=2)
print(json.dumps(summary, indent=2))
