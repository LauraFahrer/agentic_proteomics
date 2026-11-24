# Scripts — run 0007 (independent re-implementation, Python)

Self-implemented method libraries:
- `methods.py`: normalisation, imputation, MAD, Brunet NMF + consensus/cophenetic/silhouette, Kim–Park features,
  consensus clustering, moderated t (limma-type), ANOVA wrapper, BH, GSEA, Max_probe collapse, ORA, ssGSEA/PTM-SEA.
- `networks_kinases.py`: WGCNA-type network modules, kinase–substrate regulons, kinase activity (mean, KSEA, MLR),
  INKA.
- `common.py`: paths, loaders, the paper's published cluster labels.

Inputs: `data_cache/` (git-ignored) holds the quantms files, the paper supplements, QuickGO, UniProt and the
`db/` resources. They are obtained as in `runs/0006/.../scripts/00_setup.sh`. The `db/` files are copied from the
data folders of the repositories cloned there; the Rdata kinase/UniProt tables were exported to CSV with
`load()` + `write.csv()`.

    export OMP_NUM_THREADS=1 PY=/opt/mm/envs/repro/bin/python   # any Python 3.12 + numpy/scipy/pandas/matplotlib
    cd scripts
    $PY 01_prepare_data.py          # Table S3 + quantms WP (quantile normalisation, imputation)
    $PY 02_qc.py                    # Fig. S1C, Fig. 1B
    $PY 03_nmf_grid.py              # 4 datasets x 5 fractions x k=2..5, 100 runs (~40 min on 4 cores)
    $PY 04_subtypes.py              # final subtypes, integrative consensus clustering, metagenes
    $PY 05_functional.py gsea; $PY 05_functional.py wgcna; $PY 05_functional.py diff
    $PY 06_kinase.py kinome; $PY 06_kinase.py inka; $PY 06_kinase.py ptmsea; $PY 06_kinase.py activity
    $PY 06b_activity_sensitivity.py --ksea   # Fig. 8A sensitivity

Validation: `methods.moderated_t` vs R limma on simulated data, max |Δt| = 0.001 (REPORT.md section 1).
