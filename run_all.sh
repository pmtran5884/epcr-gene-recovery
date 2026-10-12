#!/usr/bin/env bash
# End-to-end run. Expect ~1-2 h wall clock per WGS run, dominated by the
# download (15-39 GB per run) and the single streaming pass over all reads.
#   conda env create -f environment.yml && conda activate procr-miniprot
#   ./run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
source config/params.sh

scripts/00_fetch_inputs.sh
scripts/01_extract_reference.sh
for acc in $RUNS; do
  scripts/02_fetch_sra.sh "$acc"
  scripts/03_reconstruct_locus.sh "$acc"
done
scripts/04_align_and_report.py
echo "[run_all] done -- see $RESULTS/report.md"
