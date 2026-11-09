# Analysis scripts — PXD048538 run 0004 (protein-protein interaction / complexes)

Co-abundance (guilt-by-co-regulation) analysis on the whole proteome. Run in order with
PXD048538_DATA set to the data dir (populated by 00_download.sh).

    export PXD048538_DATA=/path/to/data
    bash   00_download.sh "$PXD048538_DATA"   # quantms data + gene sets
    python 01_build_matrix.py                  # protein matrix
    python 02_coabundance.py                   # bicor matrix + known-complex benchmark (fig 01)
    python 03_modules.py                       # co-abundance modules + strong pairs (fig 02)
    python 04_hypotheses.py                    # STRING annotation + complex-subtype + enrichment (fig 03)

References: GO-CC complex terms (MSigDB, cached) for benchmarking; STRING REST API
(string-db.org, queried by 04) for PPI evidence; ComplexPortal (EBI) for spot checks.
Co-abundance = biweight midcorrelation across the 14 cell-line means (robust to outlier lines).
NOTE n=14: only strong edges (bicor>=0.9) are reliable (45x complex enrichment); co-abundance
is co-regulation, not proven physical interaction.
- common.py, stats.py : shared helpers.
