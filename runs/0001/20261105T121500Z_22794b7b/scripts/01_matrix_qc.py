"""Build the protein-group x sample matrix from the quantms pg.parquet and run QC.

Input : PXD048538.pg.parquet (DIA-NN PG.MaxLFQ per run), sample/run parquet
Output: tables/protein_matrix_log2.csv.gz, tables/sample_annotation.csv, QC figures
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (CELL_META, FIG, INK2, MUTED, QDIR, SUBTYPE_COLORS, TAB, savefig,
                    short_line, style)

style()
qc = {}

# ---- samples -------------------------------------------------------------
samp = pd.read_parquet(QDIR / "PXD048538.sample.parquet")
runs = pd.read_parquet(QDIR / "PXD048538.run.parquet")
sdrf = pd.read_csv(QDIR / "PXD048538.sdrf.tsv", sep="\t")
run2sample = {r.run_accession: r.samples[0]["sample_accession"] for r in runs.itertuples()}
ann = sdrf[["source name", "characteristics[cell line]", "characteristics[disease]",
            "characteristics[sex]", "characteristics[age]",
            "characteristics[biological replicate]", "comment[data file]"]].copy()
ann.columns = ["sample", "cell_line_sdrf", "disease_sdrf", "sex", "age", "bio_rep", "raw_file"]
ann["run"] = ann.raw_file.str.replace(".raw", "", regex=False)
assert set(ann.run) == set(run2sample)
assert all(run2sample[r] == s for r, s in zip(ann.run, ann["sample"]))
ann["line"] = ann.cell_line_sdrf.map(short_line)
ann = ann.join(CELL_META, on="line")
ann["sample_id"] = ann.line + "_R" + ann.bio_rep.astype(str)
order = list(CELL_META.sort_values(["paper_subtype"]).index)
ann["line"] = pd.Categorical(ann.line, categories=order, ordered=True)
ann = ann.sort_values(["line", "bio_rep"]).reset_index(drop=True)
ann.to_csv(TAB / "sample_annotation.csv", index=False)
qc["n_runs"] = len(ann)
qc["n_lines"] = ann.line.nunique()

# ---- protein groups ------------------------------------------------------
pg = pd.read_parquet(QDIR / "PXD048538.pg.parquet",
                     columns=["pg_accessions", "gg_names", "grouped_runs", "pg_qvalue",
                              "global_qvalue", "intensity", "is_decoy", "contaminant",
                              "peptide_counts"])
qc["pg_rows_total"] = len(pg)
pg["run"] = pg.grouped_runs.str[0]
pg["protein"] = pg.pg_accessions.map(lambda a: ";".join(a))
pg["gene"] = pg.gg_names.map(lambda a: ";".join(a) if a is not None and len(a) else "")
pg["n_pep"] = pg.peptide_counts.map(lambda d: d["unique_sequences"])
qc["pg_rows_contaminant"] = int(pg.contaminant.sum())
keep = (~pg.is_decoy) & (~pg.contaminant) & (pg.pg_qvalue <= 0.01) & \
       (pg.global_qvalue <= 0.01) & (pg.intensity > 0)
qc["pg_rows_after_filter"] = int(keep.sum())
pgf = pg[keep]
mat = pgf.pivot_table(index="protein", columns="run", values="intensity", aggfunc="max")
gene = pgf.groupby("protein").gene.first()
mat = np.log2(mat)
mat.columns = [ann.set_index("run").sample_id[c] for c in mat.columns]
mat = mat[ann.sample_id]
qc["protein_groups_quantified"] = int(mat.shape[0])

# ---- QC: IDs per run, missingness, intensity distribution ---------------
ids = mat.notna().sum()
qc["ids_per_run_min_median_max"] = [int(ids.min()), int(ids.median()), int(ids.max())]
med = mat.median()
qc["run_median_log2_range"] = [float(med.min()), float(med.max())]

# Median-centring: DIA-NN MaxLFQ is already cross-run normalised; we check and apply
# an additional median alignment (shift only) so that any residual loading offset is removed.
norm = mat - med + med.mean()

# replicate reproducibility
corr = norm.corr(method="pearson", min_periods=500)
same, diff = [], []
for i, a in enumerate(norm.columns):
    for b in norm.columns[i + 1:]:
        (same if a.split("_R")[0] == b.split("_R")[0] else diff).append(corr.loc[a, b])
qc["pearson_within_line_median"] = float(np.median(same))
qc["pearson_between_line_median"] = float(np.median(diff))
# CV on linear scale within cell line
cvs = []
for line, cols in ann.groupby("line", observed=True).sample_id:
    lin = 2 ** norm[list(cols)]
    ok = lin.notna().sum(axis=1) == 3
    cvs.append(pd.Series((lin[ok].std(axis=1, ddof=1) / lin[ok].mean(axis=1)).values, name=line))
qc["median_cv_within_line_percent"] = {s.name: round(float(s.median() * 100), 1) for s in cvs}

# ---- presence per cell line (>=2 of 3 replicates) -------------------------
pres = pd.DataFrame({line: norm[list(cols)].notna().sum(axis=1) >= 2
                     for line, cols in ann.groupby("line", observed=True).sample_id})
qc["proteins_present_in_all_lines"] = int(pres.all(axis=1).sum())
qc["proteins_present_in_some_lines_only"] = int(((pres.sum(axis=1) > 0) & ~pres.all(axis=1)).sum())

# cell-line-level mean matrix (mean of available replicates, require >=2 values)
line_mean = pd.DataFrame({line: norm[list(cols)].mean(axis=1).where(pres[line])
                          for line, cols in ann.groupby("line", observed=True).sample_id})

out = norm.copy()
out.insert(0, "gene", gene.reindex(out.index))
out.to_csv(TAB / "protein_matrix_log2.csv.gz", float_format="%.4f")
lm = line_mean.copy()
lm.insert(0, "gene", gene.reindex(lm.index))
lm.to_csv(TAB / "protein_matrix_cellline_mean_log2.csv.gz", float_format="%.4f")
json.dump(qc, open(TAB / "qc_summary.json", "w"), indent=2)
print(json.dumps(qc, indent=2))

# ---- figures ---------------------------------------------------------------
cols = [SUBTYPE_COLORS[s] for s in ann.paper_subtype]
fig, axes = plt.subplots(2, 1, figsize=(10, 6.2), gridspec_kw={"hspace": 0.55})
ax = axes[0]
ax.bar(range(len(ids)), ids.values, color=cols, width=0.8)
ax.set_xticks(range(len(ids)), ids.index, rotation=90, fontsize=7)
ax.set_ylabel("protein groups\n(1% FDR)")
ax.set_title("A  Protein groups quantified per run")
ax.set_ylim(0, ids.max() * 1.1)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=c, label=f"{s} subtype (Hu et al.)") for s, c in SUBTYPE_COLORS.items()],
          loc="upper right", ncol=2, fontsize=8, bbox_to_anchor=(1, 1.25))
ax = axes[1]
bp = ax.boxplot([mat[c].dropna() for c in mat.columns], showfliers=False, widths=0.6,
                patch_artist=True, medianprops={"color": "#0b0b0b"})
for p, c in zip(bp["boxes"], cols):
    p.set_facecolor(c); p.set_alpha(0.55); p.set_edgecolor(c)
ax.set_xticks(range(1, len(ids) + 1), ids.index, rotation=90, fontsize=7)
ax.set_ylabel("log2 MaxLFQ")
ax.set_title("B  Intensity distribution per run (before median alignment)")
savefig(fig, "01_qc_ids_intensity")
plt.close(fig)

# missingness vs abundance
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw={"wspace": 0.35})
nlines = pres.sum(axis=1)
ax = axes[0]
vc = nlines.value_counts().sort_index()
ax.bar(vc.index, vc.values, color="#2a78d6", width=0.8)
ax.set_xlabel("number of cell lines in which the protein is detected (≥2/3 reps)")
ax.set_ylabel("protein groups")
ax.set_title("A  Detection breadth across 14 lines")
ax = axes[1]
mean_int = norm.mean(axis=1)
frac_missing = norm.isna().mean(axis=1)
ax.scatter(mean_int, frac_missing, s=3, alpha=0.3, color="#2a78d6", linewidths=0)
ax.set_xlabel("mean log2 intensity (where observed)")
ax.set_ylabel("fraction of runs missing")
ax.set_title("B  Missingness is intensity-dependent")
savefig(fig, "02_qc_missingness")
plt.close(fig)

# correlation heatmap
fig, ax = plt.subplots(figsize=(8.2, 7))
im = ax.imshow(corr.values, cmap="Blues", vmin=np.nanpercentile(corr.values, 2), vmax=1)
ax.set_xticks(range(len(corr)), corr.columns, rotation=90, fontsize=6)
ax.set_yticks(range(len(corr)), corr.index, fontsize=6)
for i, c in enumerate(cols):
    ax.add_patch(plt.Rectangle((-2.2, i - 0.5), 1, 1, color=c, clip_on=False))
cb = fig.colorbar(im, ax=ax, shrink=0.7)
cb.set_label("Pearson r (log2 intensities)")
ax.set_title("Run-to-run correlation (median within-line r = %.3f, between-line r = %.3f)"
             % (qc["pearson_within_line_median"], qc["pearson_between_line_median"]))
savefig(fig, "03_qc_correlation")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 3))
bp = ax.boxplot([s.dropna() * 100 for s in cvs], showfliers=False, patch_artist=True,
                medianprops={"color": "#0b0b0b"})
for p, line in zip(bp["boxes"], [s.name for s in cvs]):
    c = SUBTYPE_COLORS[CELL_META.paper_subtype[line]]
    p.set_facecolor(c); p.set_alpha(0.55); p.set_edgecolor(c)
ax.set_xticks(range(1, len(cvs) + 1), [s.name for s in cvs], rotation=45)
ax.set_ylabel("CV (%) across 3 biological replicates")
ax.set_title("Within-line coefficient of variation (proteins quantified in all 3 replicates)")
savefig(fig, "04_qc_cv")
plt.close(fig)
