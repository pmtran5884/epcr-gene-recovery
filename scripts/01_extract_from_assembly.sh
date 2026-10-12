#!/usr/bin/env bash
# Recover EPCR from one genome assembly by splice-aware protein-to-genome
# alignment. Works on an assembly with no annotation of any kind -- miniprot
# infers the exon structure from the query protein, so no GFF/GTF is needed.
#
# Usage: scripts/01_extract_from_assembly.sh <label>
#        label must appear in config/targets.tsv with mode 'genome' or 'scaffold'
#
# Outputs (results/<label>/): miniprot.gff, protein.faa, cds.fna, cds_exons.bed,
# locus.txt, qc.tsv  (+ bait.fa for the 'scaffold' species, consumed by step 03)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/params.sh
LABEL="${1:?usage: 01_extract_from_assembly.sh <label from config/targets.tsv>}"

read -r MODE ASM_ACC QUERY_ACC ASM_URL < <(
  awk -F'\t' -v l="$LABEL" '$1==l {print $2"\t"$3"\t"$4"\t"$5; exit}' config/targets.tsv)
[[ -n "${MODE:-}" ]] || { echo "[01] $LABEL not found in config/targets.tsv"; exit 1; }
[[ "$MODE" == "genome" || "$MODE" == "scaffold" ]] || {
  echo "[01] $LABEL has mode '$MODE'; this script handles genome/scaffold only"; exit 1; }

OUT="$RESULTS/$LABEL"
FNA="$WORK/${LABEL}.fna"
mkdir -p "$OUT" "$WORK"

# --- assembly ---------------------------------------------------------------
if [[ ! -s "$FNA" ]]; then
  echo "[01] $LABEL: downloading $ASM_ACC"
  curl -sSL -o "$FNA.gz" "$ASM_URL"
  gunzip -f "$FNA.gz"
fi
[[ -s "$FNA.fai" ]] || samtools faidx "$FNA"

# --- query protein ----------------------------------------------------------
QFA="$RESULTS/query_${QUERY_ACC}.faa"
if [[ ! -s "$QFA" ]]; then
  echo "[01] fetching query protein $QUERY_ACC"
  curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${QUERY_ACC}&rettype=fasta&retmode=text" -o "$QFA"
  grep -q "^>" "$QFA" || { echo "[01] ERROR: query fetch failed"; exit 1; }
fi

# --- splice-aware protein-to-genome alignment -------------------------------
echo "[01] $LABEL: miniprot $(miniprot --version) -t $THREADS -j $SPLICE_MODEL"
miniprot -t "$THREADS" -j "$SPLICE_MODEL" --gff --trans --aln --outn=5 \
  "$FNA" "$QFA" > "$OUT/miniprot.gff" 2> "$OUT/miniprot.log"

# --- CDS, protein, QC (handles the N-terminal over-extension) ---------------
scripts/extract_cds_from_gff.py --genome "$FNA" --gff "$OUT/miniprot.gff" \
  --label "$LABEL" --out-dir "$OUT" --query-acc "$QUERY_ACC" --assembly-acc "$ASM_ACC"

# --- bait for the reads route ----------------------------------------------
if [[ "$MODE" == "scaffold" ]]; then
  CTG=$(awk '$1=="contig"{print $2}' "$OUT/locus.txt")
  GS=$(awk '$1=="gene_start"{print $2}' "$OUT/locus.txt")
  GE=$(awk '$1=="gene_end"{print $2}' "$OUT/locus.txt")
  BS=$(( GS - FLANK )); [[ $BS -lt 1 ]] && BS=1
  BE=$(( GE + FLANK ))
  printf 'bait_start\t%s\nbait_end\t%s\n' "$BS" "$BE" >> "$OUT/locus.txt"
  samtools faidx "$FNA" "${CTG}:${BS}-${BE}" > "$OUT/bait.fa"
  samtools faidx "$OUT/bait.fa"
  echo "[01] $LABEL: bait $(( BE - BS + 1 )) bp written for the reads route"
fi
