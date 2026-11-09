# Synthesis scripts — run 0005 (integrates runs 0001-0004)

    export PXD048538_DATA=/path/to/data   # for the gene-set files (Hallmark)
    python 01_integrate.py             # per-line axis scorecard + convergence (fig 01)
    python 02_hypothesis_register.py   # hypothesis register + prioritisation (fig 02)

01 reads the committed cell-line-mean protein matrix from run 0002 and the per-line/
per-contrast result tables of runs 0001-0004 (paths in common.py). 02 encodes the
curated hypothesis evaluation (author judgement, documented for audit). No new primary
data; this run synthesises the prior four. common.py holds paths, palette and helpers.
