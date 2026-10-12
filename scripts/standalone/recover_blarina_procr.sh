#!/usr/bin/env bash
#
# Recover the Blarina brevicauda (northern short-tailed shrew) PROCR/EPCR
# protein sequence de novo from the unannotated reference genome assembly.
#
# There is NO annotated PROCR record for this species in UniProt or NCBI.
# This script locates the gene by splice-aware protein-to-genome alignment
# of a close relative's PROCR (Sorex araneus), then translates the exon
# structure miniprot infers. The result is a gene PREDICTION.
#
# Expected output: a 242-aa protein beginning MLTPPLL... ending ...TGGRRC
#
# Prerequisites: conda with the bioconda channel, curl.
# Runtime: genome download ~1-2 min (668 MB), miniprot ~1.5 min.
set -euo pipefail

WORKDIR="blarina_procr"
mkdir -p "$WORKDIR" && cd "$WORKDIR"

# --- 1. Tooling -------------------------------------------------------------
# miniprot = splice-aware aligner of a protein onto a genome (Heng Li).
if ! command -v miniprot >/dev/null 2>&1; then
  conda install -y -c bioconda -c conda-forge miniprot
fi

# --- 2. Reference genome ----------------------------------------------------
# Blarina brevicauda assembly GCA_051295465.2 (mBlaBre1_p1.2).
GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/295/465/GCA_051295465.2_mBlaBre1_p1.2/GCA_051295465.2_mBlaBre1_p1.2_genomic.fna.gz"
if [ ! -f blarina.fna ]; then
  echo "Downloading genome (~668 MB compressed)..."
  curl -sSL -o blarina.fna.gz "$GENOME_URL"
  gunzip -f blarina.fna.gz
fi

# --- 3. Query: PROCR protein of a close relative (Sorex araneus) ------------
cat > query.fa <<'FASTA'
>Sorex_araneus_PROCR
MLTPPLLLSLLPLLLLPGGALCGQEASEGPQNLHMLQISYFRNPNQVWHHGNASLGGRPTHKLEGWGTNV
TFVQLQPLQDPESWAIMEESLRLYIEQFHSLVRLVHQERHLDFPMTVRCSLGCELPTEGSKASAFFDVTVN
GSAFVSFQPKTALWVAGSGKTSNVASYTLQQLNAYNRTRFELQEFLQDTCVQYVQEHSSVKDDGTGSQKGR
SYTSLVLGVLMGCFIIAGVAVGIFLCTGGRRRC
FASTA

# --- 4. Splice-aware alignment + translation --------------------------------
# --trans prints the translated CDS of each alignment on a "##STA" line.
# The top-ranked hit (Rank=1) is PROCR on scaffold JBKFBT020000153.1 (+ strand).
echo "Aligning PROCR onto the genome with miniprot..."
miniprot -t4 --trans blarina.fna query.fa > miniprot_out.txt 2>/dev/null

# Extract the top translated protein, drop the trailing stop (*).
awk '/^##STA/{print $2; exit}' miniprot_out.txt | sed 's/\*$//' > blarina_procr_aa.txt

# --- 5. Write the FASTA and report ------------------------------------------
PROT=$(cat blarina_procr_aa.txt)
LEN=${#PROT}
{
  echo ">Blarina_brevicauda_PROCR de novo prediction from GCA_051295465.2 (miniprot); ${LEN} aa"
  echo "$PROT" | fold -w 60
} > blarina_procr_protein.fasta

echo ""
echo "=== Recovered Blarina brevicauda PROCR ==="
echo "length: ${LEN} aa"
echo "starts with M : $([ "${PROT:0:1}" = "M" ] && echo yes || echo NO)"
cat blarina_procr_protein.fasta
