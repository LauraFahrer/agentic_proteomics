# Analysis scripts — PXD048538 run 0001 / 20261105T121500Z_22794b7b

Run in order from this directory. Set `PXD048538_DATA` to the directory that
`00_download.sh` populates (default `data_cache`). A fixed seed (20261105) is used.

```bash
export PXD048538_DATA=/path/to/data
bash   00_download.sh "$PXD048538_DATA"     # data + paper + gene sets (checks sha256)
python 01_matrix_qc.py                       # protein matrix + QC (fig 01-04)
python 02_structure.py                       # concordance, PCA, clustering, hallmarks (fig 05-08)
python 03_differential.py                    # differential + gene-set tests (fig 09-10)
python 04_targeted.py                        # amplicons, lineage TFs, EBV module, markers (fig 11-14)
bash   05_ebv_search.sh "$PXD048538_DATA/ebv" <tools_env> TMK1 SNU719 YCCEL1 AGS eAGS
python 05_ebv_summary.py "$PXD048538_DATA/ebv"   # EBV peptide search summary (fig 15)
```

Dependencies: Python 3.13 with pandas, numpy, scipy, scikit-learn, matplotlib,
statsmodels, pyarrow, mudata, anndata, openpyxl, xlrd (see ../run_info/python_packages.txt).
For step 05: thermorawfileparser and sage-proteomics (bioconda; installed via micromamba in this run).

- `common.py` — paths, cell-line annotation (incl. Hu et al. Table S1), palette, helpers.
- `stats.py`  — moderated t-test (empirical Bayes) and rank-based gene-set test.

Inputs are not committed (size + redistribution terms); `00_download.sh` fetches them.
