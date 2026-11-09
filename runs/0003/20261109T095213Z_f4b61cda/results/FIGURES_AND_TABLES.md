# Figure and table index — run 0003 (signaling / PTM)

PTM data source: publication Supplementary Table S3 (quantms e60b9e9aa4c8 has no PTM layer).

## Figures
| File | Content |
|---|---|
| 01_phospho_coverage | Phosphosite coverage per line; modified-residue composition |
| 02_diff_phospho_volcano | Differential phosphorylation EMT vs Metabolism: A protein-normalised S/T (0 sig), B tyrosine (25 sig) |
| 03_kinase_activity_subtype | phosX differential kinase activity (Ser/Thr + Tyr kinomes), ranked by KS |
| 04_signaling_modules | A differential pY sites across lines (RTK/adaptor module); B Tyr kinome; C AGC vs MAP2K/MAP3K |

## Tables
| File | Content |
|---|---|
| phospho_site_meta_{PP,YP}.csv.gz | site metadata incl. 10-mer flanking windows |
| pp_linemean.csv.gz, yp_linemean.csv.gz, pp_protnorm_linemean.csv.gz | line-mean phospho matrices |
| phospho_coverage.json | coverage counts |
| diff_PP_*.csv.gz, diff_YP_*.csv.gz | differential phosphorylation (raw + protein-normalised) |
| diff_phospho_summary.json | differential summary |
| kinase/*_kinases.tsv | phosX kinase-activity tables per contrast × layer |
| kinase_activity_summary.json, hypotheses_summary.json | summaries |
| seqrnk/*.seqrnk | ranked flanking-sequence inputs to phosX |
