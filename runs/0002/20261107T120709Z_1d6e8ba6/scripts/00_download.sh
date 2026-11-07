#!/usr/bin/env bash
# Download the quantms output (version e60b9e9aa4c8) of PXD048538, the paper's
# full text (Europe PMC) and MSigDB gene-set files used for enrichment.
# Usage: bash 00_download.sh <data_dir>
set -euo pipefail
OUT=${1:?usage: 00_download.sh <data_dir>}
mkdir -p "$OUT/quantms" "$OUT/lit" "$OUT/genesets"
B=https://browse.quantms.org/quantms/datasets/PXD048538/e60b9e9aa4c8
for f in PXD048538.dataset.parquet PXD048538.ontology.parquet PXD048538.provenance.parquet \
         PXD048538.run.parquet PXD048538.sample.parquet PXD048538.sdrf.tsv provenance.json \
         PXD048538.pg.parquet PXD048538.h5mu PXD048538_msstats_in.csv.gz; do
  curl -sSL --retry 3 -o "$OUT/quantms/$f" "$B/$f"
done
# verify checksums against provenance.json
python3 - "$OUT/quantms" <<'EOF'
import hashlib, json, sys, os
d = sys.argv[1]
prov = json.load(open(os.path.join(d, "provenance.json")))
for f in prov["files"]:
    p = os.path.join(d, f["path"])
    if not os.path.exists(p):
        print("missing (not downloaded):", f["path"]); continue
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    print(("OK  " if h == f["sha256"] else "BAD ") + f["path"])
EOF
curl -sS -o "$OUT/lit/PMC11161736.xml" \
  "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11161736/fullTextXML"
# MSigDB gene sets (mirrored copies on GitHub; msigdb.org is not reachable from the sandbox)
curl -sSL -o "$OUT/genesets/h.all.v2023.1.Hs.symbols.gmt" \
  https://raw.githubusercontent.com/saezlab/VisiumMS/HEAD/config/h.all.v2023.1.Hs.symbols.gmt
for f in c2.cp.kegg.v7.5.1.symbols.gmt c2.cp.reactome.v7.5.1.symbols.gmt c5.go.bp.v7.5.1.symbols.gmt c5.go.cc.v7.5.1.symbols.gmt; do
  curl -sSL -o "$OUT/genesets/$f" "https://raw.githubusercontent.com/MLO-lab/MuVI/HEAD/msigdb/$f"
done
# HGNC gene annotation (cytogenetic location, used for the amplicon analysis)
curl -sSL -o "$OUT/genesets/hgnc_complete_set.txt" \
  https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
# Paper supplementary tables (Table S1 cell-line features, Table S3 processed WP matrix)
mkdir -p "$OUT/supp"
curl -sSL -o "$OUT/supp/supp.zip" "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11161736/supplementaryFiles"
unzip -o -q "$OUT/supp/supp.zip" -d "$OUT/supp/x"
