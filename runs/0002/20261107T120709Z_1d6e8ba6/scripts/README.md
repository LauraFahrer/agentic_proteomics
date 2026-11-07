# Analysis scripts — PXD048538 run 0002 (pathway analysis)

Run in order; set PXD048538_DATA to the data directory populated by 00_download.sh
(shared with run 0001). Fixed seed 20261105.

    export PXD048538_DATA=/path/to/data
    bash   00_download.sh "$PXD048538_DATA"   # data + gene sets (checks sha256)
    python 01_build_matrix.py                  # protein matrix
    python 02_pathway_landscape.py             # pathway modules (fig 01-02)
    python 03_enrichment.py                    # differential + preranked GSEA (fig 03-04; ~9 min)
    python 04_hypotheses.py                    # H1-H3 + STRING enrichment (fig 05-06)

- common.py : paths, cell-line annotation, palette, read_gmt, helpers.
- stats.py  : moderated t-test, rank-based set test, preranked GSEA (gsea_preranked).

Deps: Python 3.13 + pandas, numpy, scipy, scikit-learn, matplotlib, statsmodels, pyarrow
(see ../run_info/python_packages.txt). 04 queries string-db.org (needs network).
Inputs are not committed; 00_download.sh fetches and checksum-verifies them.
