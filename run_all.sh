#!/usr/bin/env bash
# End-to-end run. Expect ~1-2 h wall clock per WGS run, dominated by the
# download (15-39 GB per run) and the single streaming pass over all reads.
#   conda env create -f environment.yml && conda activate procr-miniprot
#   ./run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
source config/params.sh

# 1. the scaffold species: assembly + annotation, also supplies the bait
scripts/01_extract_from_assembly.sh "$SCAFFOLD_LABEL"

# 2. the two shrews: unannotated assemblies, miniprot finds the gene directly
for label in $GENOME_TARGETS; do
  scripts/01_extract_from_assembly.sh "$label"
done

# 3. the taxon with no assembly: reconstruct the locus from WGS reads
for acc in $RUNS; do
  scripts/02_fetch_sra.sh "$acc"
  scripts/03_reconstruct_locus.sh "$acc"
done

# 4. align everything and write the tables
scripts/04_align_and_report.py
echo "[run_all] done -- see $RESULTS/report.md"
