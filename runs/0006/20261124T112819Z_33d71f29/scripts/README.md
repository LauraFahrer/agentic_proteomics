# Scripts — run 0006 (reproduction of Hu et al. 2024, PXD048538)

Run from the run directory. `00_setup.sh` fetches the data into `data_cache/` (git-ignored), clones the original
tools to `/opt/src` (pinned commits) and builds the R 4.3/Bioconductor + GSEA 4.3.2 environment with micromamba.

    bash scripts/00_setup.sh
    export PATH=/opt/mm/envs/repro/bin:$PATH
    python  scripts/01_prepare_data.py          # Table S3 matrices + quantms WP (quantile norm, imputation)
    Rscript scripts/02_qc.R                      # YP annotation, replicate correlation (S1C), GO CC (1B)
    for d in WP PP YP WPq; do Rscript scripts/03_nmf_grid.R $d mean; done   # NMF grid (~1.5 h on 4 cores)
    Rscript scripts/03_nmf_grid.R WP median      # median-AD sensitivity
    Rscript scripts/04_subtypes.R                # final subtypes, ConsensusClusterPlus integration
    python  scripts/04a_ptm_windows.py           # flanking windows for PTM-SEA
    R_LIBS=/opt/stubs Rscript scripts/05_functional.R   # GSEA, WGCNA, ANOVA/limma + GO, Fig. 5B
    Rscript scripts/05b_gsea_permutation_sensitivity.R
    Rscript scripts/05c_compare_gsea.R
    Rscript scripts/08_kinome.R                  # kinome, outliers, kinome NMF, 40 kinases

Prepared but NOT executed in this run (running cloned GitHub code was refused by the session permission policy):

    Rscript scripts/06_pp_kinase_activity.R      # PhosMap mean/KSEA/MLR x raw/normalised PP (Fig. 8A)
    python  scripts/07a_inka_inputs.py && Rscript scripts/07b_inka_run.R   # InKA on YP (Fig. 7A/D)
    # PTM-SEA: ssGSEA2.0 ssgsea-cli.R with PTMsigDB v1.9.1 flanking gmt on PP_windows (Fig. S5D) - no driver yet

`/opt/stubs` holds empty placeholder packages GO.db and GenomeInfoDbData that only satisfy WGCNA's namespace
imports (bioconductor.org unreachable); no GO.db functionality is used.
`common.R` holds paths, loaders, the paper's cluster labels (read from the figures) and helpers.
