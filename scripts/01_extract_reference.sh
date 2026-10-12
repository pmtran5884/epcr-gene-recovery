#!/usr/bin/env bash
# miniprot the query protein against the reference assembly, then cut a bait
# locus (gene +/- FLANK) that later steps use to capture reads.
# Outputs (results/$REF_LABEL/):
#   miniprot.gff   GFF3 in genome coordinates, with ##PAF and ##STA records
#   protein.faa    translated CDS of the best hit
#   cds_exons.bed  CDS exons of the best hit (for depth QC)
#   locus.txt      contig, gene start/end, bait start/end
#   bait.fa        the bait sequence
set -euo pipefail
cd "$(dirname "$0")/.."
source config/params.sh
FNA="$WORK/${REF_LABEL}_genomic.fna"
OUT="$RESULTS/$REF_LABEL"
mkdir -p "$OUT"

echo "[01] miniprot $(miniprot --version) -t $THREADS against $ASM_ACC"
/usr/bin/time -l miniprot -t "$THREADS" -j "$SPLICE_MODEL" --gff --trans --aln --outn=5 \
    "$FNA" "$RESULTS/query_${QUERY_ACC}.faa" > "$OUT/miniprot.gff" 2> "$OUT/miniprot.log" \
  || miniprot -t "$THREADS" -j "$SPLICE_MODEL" --gff --trans --aln --outn=5 \
    "$FNA" "$RESULTS/query_${QUERY_ACC}.faa" > "$OUT/miniprot.gff" 2> "$OUT/miniprot.log"

# Best hit = the Rank=1 mRNA; miniprot prints hits in descending score order.
awk -F'\t' '$3=="mRNA" && $9 ~ /Rank=1/ {print; exit}' "$OUT/miniprot.gff" > "$OUT/best_mrna.gff"
CTG=$(cut -f1 "$OUT/best_mrna.gff")
GS=$(cut -f4 "$OUT/best_mrna.gff")
GE=$(cut -f5 "$OUT/best_mrna.gff")
IDENT=$(sed 's/.*Identity=\([0-9.]*\).*/\1/' "$OUT/best_mrna.gff")

awk -F'\t' -v c="$CTG" '$1==c && $3=="CDS" && $9 ~ /Rank=1/ {print $1"\t"$4-1"\t"$5}' \
  "$OUT/miniprot.gff" | sort -k2,2n > "$OUT/cds_exons.bed"

# Translated protein of the first (best) alignment.
awk '/^##STA/{print ">'"${REF_LABEL}"'_EPCR_miniprot '"${ASM_ACC}"' '"${CTG}:${GS}-${GE}"'"; print $2; exit}' \
  "$OUT/miniprot.gff" | sed 's/\*$//' > "$OUT/protein.faa"

BS=$(( GS - FLANK )); [[ $BS -lt 1 ]] && BS=1
BE=$(( GE + FLANK ))
printf 'contig\t%s\ngene_start\t%s\ngene_end\t%s\nbait_start\t%s\nbait_end\t%s\nidentity_to_query\t%s\n' \
  "$CTG" "$GS" "$GE" "$BS" "$BE" "$IDENT" > "$OUT/locus.txt"

samtools faidx "$FNA" "${CTG}:${BS}-${BE}" > "$OUT/bait.fa"
samtools faidx "$OUT/bait.fa"

echo "[01] best hit $CTG:$GS-$GE  identity=$IDENT  exons=$(wc -l < "$OUT/cds_exons.bed" | tr -d ' ')"
echo "[01] protein $(grep -v '^>' "$OUT/protein.faa" | tr -d '\n' | wc -c | tr -d ' ') aa; bait $(( BE - BS + 1 )) bp"
