#!/usr/bin/env bash
# Direct test for Epstein-Barr virus (EBV) proteins in the raw MS data.
# The quantms reprocessing searched a human-only database, so viral proteins cannot appear
# in its output. Here the whole-proteome DDA spectral-library runs of selected lines are
# re-searched against human (UniProt Swiss-Prot, UP000005640) + all reviewed EBV proteins.
#
# Requirements (installed into a conda env via micromamba in this run):
#   thermorawfileparser (bioconda), sage-proteomics 0.14.x (bioconda)
# Usage: bash 05_ebv_search.sh <work_dir> <tools_env_dir> <line> [<line> ...]
set -euo pipefail
WORK=${1:?work dir}; ENV=${2:?tools env}; shift 2
mkdir -p "$WORK/db" "$WORK/raw" "$WORK/mzml" "$WORK/sage"
if [ ! -s "$WORK/db/human_ebv.fasta" ]; then
  curl -sS -o "$WORK/db/human_sp.fasta" \
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000005640%29%20AND%20%28reviewed%3Atrue%29"
  curl -sS -o "$WORK/db/ebv.fasta" \
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A10376%29%20AND%20%28reviewed%3Atrue%29"
  cat "$WORK/db/human_sp.fasta" "$WORK/db/ebv.fasta" > "$WORK/db/human_ebv.fasta"
fi
for L in "$@"; do
  RAW="$WORK/raw/Gastric_WP_DDAlibrary_${L}.raw"
  if [ ! -s "$RAW" ]; then
    curl -sS -C - --retry 20 -o "$RAW" \
      "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/12/PXD048538/Gastric_WP_DDAlibrary_${L}.raw"
  fi
  MZ="$WORK/mzml/Gastric_WP_DDAlibrary_${L}.mzML"
  [ -s "$MZ" ] || "$ENV/bin/thermorawfileparser" -i "$RAW" -o "$WORK/mzml" -f 2 -L 2 >/dev/null
done
cat > "$WORK/sage/config.json" <<EOF
{
  "database": {
    "bucket_size": 32768,
    "enzyme": {"missed_cleavages": 2, "min_len": 7, "max_len": 40, "cleave_at": "KR", "restrict": "P"},
    "fragment_min_mz": 150, "fragment_max_mz": 2000,
    "peptide_min_mass": 500, "peptide_max_mass": 5000,
    "ion_kinds": ["b", "y"], "min_ion_index": 2,
    "static_mods": {"C": 57.021464},
    "variable_mods": {"M": [15.994915], "^": [42.010565]},
    "max_variable_mods": 2,
    "decoy_tag": "rev_", "generate_decoys": true,
    "fasta": "$WORK/db/human_ebv.fasta"
  },
  "precursor_tol": {"ppm": [-10, 10]},
  "fragment_tol": {"ppm": [-20, 20]},
  "isotope_errors": [0, 1],
  "deisotope": true, "chimera": false, "wide_window": false,
  "predict_rt": true, "min_peaks": 15, "max_peaks": 150, "report_psms": 1,
  "mzml_paths": [$(for L in "$@"; do printf '"%s",' "$WORK/mzml/Gastric_WP_DDAlibrary_${L}.mzML"; done | sed 's/,$//')]
}
EOF
"$ENV/bin/sage" "$WORK/sage/config.json" -o "$WORK/sage/out" --write-pin
