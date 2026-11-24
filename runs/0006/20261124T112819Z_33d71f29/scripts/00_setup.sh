#!/usr/bin/env bash
# Run 0006 - fetch inputs and tools for the reproduction of Hu et al. 2024 (Mol Oncol 18:1437).
# Usage: bash scripts/00_setup.sh            (run from the run directory)
# Creates: data_cache/ (git-ignored), /opt/mm/envs/repro (R 4.3 + Bioconductor), /opt/src (tools)
set -euo pipefail
RUN=$(cd "$(dirname "$0")/.." && pwd)
DC=$RUN/data_cache; mkdir -p "$DC"/{quantms,lit,supp,go,uniprot}
SRC=${REPRO_SRC:-/opt/src}; mkdir -p "$SRC"

# ---- 1. quantms reprocessing (WP only, DIA-NN 2.5.1) ----
B=https://browse.quantms.org/quantms/datasets/PXD048538/e60b9e9aa4c8
for f in PXD048538.dataset.parquet PXD048538.ontology.parquet PXD048538.provenance.parquet \
         PXD048538.run.parquet PXD048538.sample.parquet PXD048538.sdrf.tsv provenance.json \
         PXD048538.pg.parquet PXD048538_msstats_in.csv.gz; do
  curl -sSL --retry 3 -o "$DC/quantms/$f" "$B/$f"
done
(cd "$DC/quantms" && python3 -c "
import json,hashlib,os
for f in json.load(open('provenance.json'))['files']:
    if os.path.exists(f['path']):
        h=hashlib.sha256(open(f['path'],'rb').read()).hexdigest(); print('OK ' if h==f['sha256'] else 'BAD', f['path'])")

# ---- 2. publication full text + supplementary files (Table S3 = author-processed WP/PP/YP) ----
curl -sS -o "$DC/lit/PMC11161736.xml" "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11161736/fullTextXML"
curl -sSL -o "$DC/supp/supp.zip" "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11161736/supplementaryFiles"
unzip -o -q "$DC/supp/supp.zip" -d "$DC/supp/x"

# ---- 3. GO cellular-component annotations (QuickGO; geneontology.org/bioconductor are blocked) ----
for go in GO:0005634 GO:0005730 GO:0031012 GO:0005840 GO:0005783 GO:0005794 GO:0005829 GO:0016020 \
          GO:0005856 GO:0005813 GO:0070062 GO:0005764 GO:0005777 GO:0005739; do
  curl -s -H "Accept: text/tsv" -o "$DC/go/${go/:/_}_exact.tsv" \
    "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?goId=$go&taxonId=9606&goUsage=exact&geneProductType=protein&geneProductSubset=Swiss-Prot&downloadLimit=2000000&selectedFields=symbol&selectedFields=goId"
done

# ---- 4. UniProt kinase family annotation (Manning groups via UniProt 'family' comments) ----
curl -sS -o "$DC/uniprot/human_kinase_families.tsv" \
  "https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession,gene_primary,protein_families&query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28keyword%3AKW-0418%29"

# ---- 5. tools (original implementations, pinned commits) ----
clone() { [ -d "$SRC/$2" ] || git clone -q https://github.com/$1/$2 "$SRC/$2"; git -C "$SRC/$2" checkout -q "$3"; }
clone broadinstitute ssGSEA2.0      f6c76af4977322b1708ebdc6dfc077d8a9bcdaf8  # PTM-SEA + MSigDB v7.0 + PTMsigDB
clone ecnuzdd        PhosMap        354539869e3dcdb7eeaffa2fd9bdb0c05e4a490c  # mean / KSEA / MLR kinase activity
clone digitalbiologylab inka-intensity 782c1aba10430e3c76a0dddf2b8a6b8e95998b18 # INKA 1.6.0 (InKAscore 1.5.0)
clone symbioticMe    proBatch       524249ee9db539ca0e93442ca698bb14cf827d90  # not installed (see report)

# ---- 6. R / Bioconductor / GSEA environment (conda-forge + bioconda) ----
if [ ! -x /opt/mm/envs/repro/bin/Rscript ]; then
  mkdir -p /opt/mm && cd /opt/mm
  curl -sSL -o mm.tar.bz2 https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.3.3-0.tar.bz2
  tar xjf mm.tar.bz2 bin/micromamba
  MAMBA_ROOT_PREFIX=/opt/mm/root ./bin/micromamba create -y -p /opt/mm/envs/repro -c conda-forge -c bioconda \
    --channel-priority strict r-base=4.3 r-nmf bioconductor-consensusclusterplus r-wgcna bioconductor-limma \
    bioconductor-fgsea bioconductor-clusterprofiler r-cluster r-data.table r-ggplot2 r-pheatmap r-readxl \
    r-matrixstats r-optparse r-doparallel r-verification r-pls r-glmnet r-mass r-entropy r-gtools \
    bioconductor-sva bioconductor-complexheatmap bioconductor-phosr bioconductor-cmapr bioconductor-impute \
    bioconductor-preprocesscore r-gplots r-ggrepel r-reshape2 r-stringr r-dplyr r-tidyr r-jsonlite r-openxlsx \
    r-arrow r-splitstackshape r-network r-lme4 bioconductor-pvca r-gridextra gsea \
    python=3.12 pandas numpy scipy scikit-learn matplotlib seaborn pyarrow openpyxl xlrd statsmodels
fi
echo "setup done"
