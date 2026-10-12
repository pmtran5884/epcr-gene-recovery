#!/usr/bin/env bash
# Reconstruct the gene locus of a taxon that has reads but no assembly, then run
# miniprot on the reconstruction.
#
# Two-pass read selection, which is the part that matters:
#   pass 1  stream every read of the run against the small bait locus. This is
#           deliberately promiscuous: a 50 kb reference has no competing repeat
#           copies, so off-target repeat reads pile onto it and MAPQ stays high
#           (we measured 825x over PROCR against an expected ~14x). Pass 1 is a
#           sensitivity filter only -- never call a consensus from it.
#   pass 2  re-map the pass-1 survivors COMPETITIVELY against the whole genome
#           and keep only MAPQ >= MIN_MAPQ alignments that land back in the
#           locus. Depth then drops to the run's genome-wide coverage, which is
#           the sanity check that the capture is specific.
#
# Usage: scripts/03_reconstruct_locus.sh SRR3948161
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source config/params.sh
ACC="${1:?usage: 03_reconstruct_locus.sh <SRR accession>}"
FNA="$WORK/${REF_LABEL}.fna"
REF="$RESULTS/$REF_LABEL"
OUT="$RESULTS/$ACC"
TMP="$WORK/$ACC"
mkdir -p "$OUT" "$TMP/fqtmp"

CTG=$(awk '$1=="contig"{print $2}' "$REF/locus.txt")
GS=$(awk '$1=="gene_start"{print $2}' "$REF/locus.txt")
GE=$(awk '$1=="gene_end"{print $2}' "$REF/locus.txt")
BS=$(awk '$1=="bait_start"{print $2}' "$REF/locus.txt")
BE=$(awk '$1=="bait_end"{print $2}' "$REF/locus.txt")
BAIT_REGION="${CTG}:${BS}-${BE}"
T2=$(( THREADS > 6 ? THREADS - 3 : 2 ))   # threads for each side of the capture pipe

# sra-tools must run offline: without this it tries to re-resolve the accession
# through the name-resolver service even when handed a local file path.
export NCBI_SETTINGS="$ROOT/$WORK/ncbi.mkfg"
if [[ ! -s "$NCBI_SETTINGS" ]]; then
  cat > "$NCBI_SETTINGS" <<EOF
/LIBS/GUID = "$(python3 -c 'import uuid;print(uuid.uuid4())')"
/repository/remote/disabled = "true"
/repository/user/cache-disabled = "true"
EOF
fi

# ---- pass 1: promiscuous capture against the bait -------------------------
if [[ ! -s "$TMP/capture.bam" ]]; then
  echo "[03] $ACC pass 1: streaming reads against the bait"
  fasterq-dump -e "$T2" -m 1G --split-spot --skip-technical --stdout \
      -t "$TMP/fqtmp" "$WORK/sra/${ACC}.sra" 2> "$OUT/fasterq.log" \
    | minimap2 -ax sr -t "$T2" "$REF/bait.fa" - 2> "$OUT/minimap2_capture.log" \
    | samtools view -b -F 4 - > "$TMP/capture.bam"
fi
samtools fastq -@ 2 -0 /dev/null -s /dev/null "$TMP/capture.bam" > "$TMP/candidates.fq" 2>/dev/null
N_CAND=$(( $(wc -l < "$TMP/candidates.fq") / 4 ))
echo "[03] pass 1 kept $N_CAND candidate reads"

# ---- pass 2: competitive re-mapping against the whole genome --------------
if [[ ! -s "$OUT/locus.bam" ]]; then
  echo "[03] $ACC pass 2: competitive re-mapping against $ASM_ACC"
  minimap2 -ax sr -t "$THREADS" -I 8G "$FNA" "$TMP/candidates.fq" 2> "$OUT/minimap2_genome.log" \
    | samtools sort -@ 4 -m 1G -o "$TMP/cand_vs_genome.bam"
  samtools index "$TMP/cand_vs_genome.bam"
  samtools view -b -q "$MIN_MAPQ" "$TMP/cand_vs_genome.bam" "$BAIT_REGION" > "$OUT/locus.bam"
  samtools index "$OUT/locus.bam"
fi
N_LOCUS=$(samtools view -c "$OUT/locus.bam")

# ---- depth QC over the reference CDS exons -------------------------------
: > "$OUT/exon_depth.tsv"
printf 'exon\tstart\tend\tbp\tmean_depth\tmin_depth\tbp_below_min_depth\n' >> "$OUT/exon_depth.tsv"
i=0
while read -r c s e _rest; do
  i=$(( i + 1 ))
  samtools depth -a -r "${c}:$(( s + 1 ))-${e}" "$OUT/locus.bam" \
    | awk -v i="$i" -v s="$(( s + 1 ))" -v e="$e" -v d="$MIN_DEPTH" \
      '{t+=$3;n++; if(min==""||$3<min)min=$3; if($3<d)lo++}
       END{printf "%d\t%d\t%d\t%d\t%.1f\t%d\t%d\n", i,s,e,n,(n?t/n:0),(min==""?0:min),lo+0}' \
    >> "$OUT/exon_depth.tsv"
done < "$REF/cds_exons.bed"
GENE_DEPTH=$(samtools depth -a -r "${CTG}:${GS}-${GE}" "$OUT/locus.bam" | awk '{t+=$3;n++} END{printf "%.2f", (n?t/n:0)}')

# ---- consensus + miniprot on the reconstruction --------------------------
# Two flavours of the same call, same coordinates, differing only in how
# heterozygous sites are written:
#   locus_consensus.fa        majority base  -> the sequence we translate
#   locus_consensus_ambig.fa  IUPAC code (-A) -> used only to flag het sites
# Translating the IUPAC version would put X in the protein at every het site,
# which is why the protein comes from the majority call.
# The quality floors matter: with no --min-BQ, samtools consensus -A flagged six
# CDS positions of SRR2917427 as heterozygous for which `samtools mpileup -Q 15`
# shows zero mismatching bases -- low-quality base noise promoted to IUPAC codes.
# Using the same floors as the variant caller removes all six and leaves the one
# site that is genuinely 19/41 reads.
echo "[03] calling consensus (min-depth $MIN_DEPTH, call-fract $CALL_FRACT, min-BQ 15, min-MQ $MIN_MAPQ)"
samtools consensus -@ 4 -a --mode simple --min-depth "$MIN_DEPTH" --call-fract "$CALL_FRACT" \
  --min-BQ 15 --min-MQ "$MIN_MAPQ" \
  -r "$BAIT_REGION" -f fasta -o "$TMP/consensus_major.fa" "$OUT/locus.bam"
{ echo ">${TARGET_LABEL}_${GENE_LABEL:-PROCR}_locus_consensus_${ACC}"; grep -v '^>' "$TMP/consensus_major.fa"; } \
  > "$OUT/locus_consensus.fa"
samtools consensus -@ 4 -a -A --mode simple --min-depth "$MIN_DEPTH" --call-fract "$CALL_FRACT" \
  --min-BQ 15 --min-MQ "$MIN_MAPQ" \
  -r "$BAIT_REGION" -f fasta -o "$TMP/consensus_ambig.fa" "$OUT/locus.bam"
{ echo ">${TARGET_LABEL}_${GENE_LABEL:-PROCR}_locus_consensus_ambig_${ACC}"; grep -v '^>' "$TMP/consensus_ambig.fa"; } \
  > "$OUT/locus_consensus_ambig.fa"

miniprot -t "$THREADS" -j "$SPLICE_MODEL" --gff --trans --aln --outn=3 \
  "$OUT/locus_consensus.fa" "$RESULTS/query_${QUERY_ACC}.faa" > "$OUT/miniprot.gff" 2> "$OUT/miniprot.log"
# The -A flavour is NOT coordinate-compatible with the majority call (ambiguous
# insertions are retained, so the two differ in length by ~150-200 bp over 53 kb).
# Give it its own miniprot pass rather than lifting coordinates across flavours.
miniprot -t "$THREADS" -j "$SPLICE_MODEL" --gff --trans --aln --outn=3 \
  "$OUT/locus_consensus_ambig.fa" "$RESULTS/query_${QUERY_ACC}.faa" \
  > "$OUT/miniprot_ambig.gff" 2> "$OUT/miniprot_ambig.log"
# miniprot's own translation of the majority-base consensus. At a 50/50
# heterozygous site the majority call is N (no base reaches --call-fract), so
# this file can carry an X; step 04 writes protein.faa, where such codons are
# resolved IUPAC-aware and only genuinely coding-changing het sites stay X.
awk '/^##STA/{print ">'"${TARGET_LABEL}"'_EPCR_miniprot_majority '"${ACC}"'_reconstruction"; print $2; exit}' \
  "$OUT/miniprot.gff" | sed 's/\*$//' > "$OUT/protein_miniprot_majority.faa"

# ---- independent variant call over the locus -----------------------------
bcftools mpileup -f "$FNA" -r "$BAIT_REGION" -q "$MIN_MAPQ" -Q 15 -d 200 -a AD -Ou "$OUT/locus.bam" 2>/dev/null \
  | bcftools call -mv -Oz -o "$OUT/locus.vcf.gz" 2>/dev/null
bcftools index -f -t "$OUT/locus.vcf.gz"
N_VAR=$(bcftools view -H "$OUT/locus.vcf.gz" | wc -l | tr -d ' ')
N_SNV=$(bcftools view -H -v snps "$OUT/locus.vcf.gz" | wc -l | tr -d ' ')
N_IND=$(bcftools view -H -v indels "$OUT/locus.vcf.gz" | wc -l | tr -d ' ')
N_GENE=$(bcftools view -H -r "${CTG}:${GS}-${GE}" "$OUT/locus.vcf.gz" | wc -l | tr -d ' ')
: > "$OUT/exonic_variants.tsv"
while read -r c s e _rest; do
  bcftools view -H -r "${c}:$(( s + 1 ))-${e}" "$OUT/locus.vcf.gz" | cut -f1,2,4,5,6 >> "$OUT/exonic_variants.tsv"
done < "$REF/cds_exons.bed"
N_EX=$(wc -l < "$OUT/exonic_variants.tsv" | tr -d ' ')

printf 'accession\t%s\ncandidate_reads_pass1\t%s\nlocus_reads_pass2_q%s\t%s\ngene_mean_depth\t%s\nlocus_variants\t%s\nlocus_snvs\t%s\nlocus_indels\t%s\ngene_span_variants\t%s\ncds_exonic_variants\t%s\n' \
  "$ACC" "$N_CAND" "$MIN_MAPQ" "$N_LOCUS" "$GENE_DEPTH" "$N_VAR" "$N_SNV" "$N_IND" "$N_GENE" "$N_EX" \
  > "$OUT/qc.tsv"
cat "$OUT/qc.tsv"
echo "[03] $ACC done: $OUT/protein.faa"
