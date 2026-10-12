#!/usr/bin/env bash
# Fetch the two inputs that come straight from NCBI: the query protein and the
# reference genome assembly. Both steps are skipped if the output already exists.
set -euo pipefail
cd "$(dirname "$0")/.."
source config/params.sh
mkdir -p "$WORK" "$RESULTS/$REF_LABEL"

QFA="$RESULTS/query_${QUERY_ACC}.faa"
if [[ ! -s "$QFA" ]]; then
  echo "[00] fetching query protein $QUERY_ACC from NCBI"
  curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${QUERY_ACC}&rettype=fasta&retmode=text" \
    -o "$QFA"
  grep -q "^>" "$QFA" || { echo "[00] ERROR: query protein fetch failed"; exit 1; }
fi
echo "[00] query protein: $(grep -v '^>' "$QFA" | tr -d '\n' | wc -c | tr -d ' ') aa"

GZ="$WORK/${REF_LABEL}_genomic.fna.gz"
FNA="$WORK/${REF_LABEL}_genomic.fna"
if [[ ! -s "$FNA" ]]; then
  [[ -s "$GZ" ]] || { echo "[00] downloading $ASM_ACC ($ASM_NAME)"; curl -sSL -o "$GZ" "$ASM_URL"; }
  echo "[00] decompressing assembly"
  gunzip -c "$GZ" > "$FNA"
fi
[[ -s "$FNA.fai" ]] || samtools faidx "$FNA"
echo "[00] assembly: $(awk '{s+=$2} END{print s}' "$FNA.fai") bp in $(wc -l < "$FNA.fai" | tr -d ' ') sequences"
