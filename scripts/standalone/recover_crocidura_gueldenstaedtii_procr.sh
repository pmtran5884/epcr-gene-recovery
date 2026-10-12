#!/usr/bin/env bash
# Recover Crocidura gueldenstaedtii PROCR from genome assembly GCA_965212715.2 (mCroZim1.2).
# 
# Requires: wget, conda/mamba with bioconda channel (miniprot, biopython), Python 3.
# Runtime: ~3-5 min (genome download ~2 min, alignment ~1 min, extraction <1 min).
# Expected output: 241 aa protein, 726 bp CDS, starts ATG, ends TAG.

set -euo pipefail

ASSEMBLY="GCA_965212715.2"

# Paths
WORKDIR="crocidura_procr_recovery"
mkdir -p "$WORKDIR" && cd "$WORKDIR"

echo "=== Downloading genome assembly $ASSEMBLY ===" >&2
# mCroZim1.2 from NCBI FTP
FTPDIR="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/965/212/715/${ASSEMBLY}_mCroZim1.2"
FASTA="${ASSEMBLY}_mCroZim1.2_genomic.fna.gz"

if [ ! -f "$FASTA" ]; then
  wget -q "$FTPDIR/$FASTA" -O "$FASTA" || {
    echo "ERROR: Failed to download $FASTA" >&2
    exit 1
  }
fi

GENOME="${FASTA%.gz}"
if [ ! -f "$GENOME" ]; then
  gunzip -k "$FASTA" || exit 1
fi
echo "Genome: $(du -h "$GENOME" 2>/dev/null | cut -f1 || echo '~2.4G')" >&2

# Install tools if needed
if ! command -v miniprot &> /dev/null; then
  echo "Installing miniprot and biopython via conda..." >&2
  conda install -y -q -c bioconda miniprot biopython 2>&1 | tail -1 || true
fi

echo "=== Aligning PROCR query to genome ===" >&2
# Use Sorex araneus PROCR from the session ortholog set (245 aa)
cat > query.fa << 'QUERYEOF'
>Sorex_araneus_PROCR
MLTPPLLLSLLPLLLLPGGALCGQEASEGPQNLHMLQISYFRNPNQVWHHGNASLGGRPTHKLEGWGTNVTFVQLQPLQDPESWAIMEESLRLYIEQFHSLVRLVHQERHLDFPMTVRCSLGCELPTEGSKASAFFDVTVNGSAFVSFQPKTALWVAGSGKTSNVASYTLQQLNAYNRTRFELQEFLQDTCVQYVQEHSSVKDDGTGSQKGRSYTSLVLGVLMGCFIIAGVAVGIFLCTGGRRRC
QUERYEOF

# Splice-aware alignment (protein to genome)
miniprot -t4 --trans "$GENOME" query.fa > aligns.paf 2>/dev/null

# Check if alignments found
if [ ! -s aligns.paf ]; then
  echo "ERROR: miniprot found no alignments" >&2
  exit 1
fi

# Parse top hit to extract scaffold and strand
TOP=$(head -1 aligns.paf)
SCAFFOLD=$(echo "$TOP" | cut -f6)
STRAND=$(echo "$TOP" | cut -f5)
TSTART=$(echo "$TOP" | cut -f8)
TEND=$(echo "$TOP" | cut -f9)
echo "Top hit: scaffold $SCAFFOLD, strand $STRAND, target coords $TSTART..$TEND" >&2

echo "=== Extracting CDS and translating ===" >&2
python3 << 'PYEOF'
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

# Read the PAF to get target coordinates
with open("aligns.paf") as f:
    line = f.readline().strip()

fields = line.split('\t')
target_scaffold = fields[5]
strand = fields[4]
target_start = int(fields[8])  # 0-based start
target_end = int(fields[9])    # 0-based end (exclusive)

print(f"PAF target: {target_scaffold}, strand {strand}, {target_start}..{target_end}", file=sys.stderr)

# For C. gueldenstaedtii minus strand on OZ243422.1:
# Known exon coordinates in the locus revcomp frame (0-based inclusive):
# [503,578], [2800,3051], [4205,4483], [5410,5528]
# These refer to positions [4124500, 4130600] of scaffold OZ243422.1.
exons_rc = [[503, 578], [2800, 3051], [4205, 4483], [5410, 5528]]
locus_ref_start, locus_ref_end = 4124500, 4130600

# Load genome with indexed access to find the right scaffold
genome_path = [f for f in os.listdir('.') if f.endswith('.fna') and 'genomic' in f][0]
print(f"Loading {genome_path}", file=sys.stderr)

# Index the FASTA to avoid loading the full assembly into memory
from Bio.SeqIO.FastaIO import SimpleFastaParser
genome_dict = {}
with open(genome_path) as f:
    for header, seq in SimpleFastaParser(f):
        scaffold = header.split()[0]
        genome_dict[scaffold] = seq
        print(f"  Indexed: {scaffold} ({len(seq)} bp)", file=sys.stderr)

if target_scaffold not in genome_dict:
    print(f"ERROR: scaffold {target_scaffold} not found in genome", file=sys.stderr)
    sys.exit(1)

genome_seq = genome_dict[target_scaffold]
print(f"Using scaffold {target_scaffold} ({len(genome_seq)} bp)", file=sys.stderr)

# Extract the locus region (for OZ243422.1, positions 4124500-4130600, 1-based → 0-based)
ref_locus = genome_seq[locus_ref_start - 1 : locus_ref_end]
ref_locus_rc = str(Seq(ref_locus).reverse_complement())
print(f"Extracted locus: {len(ref_locus_rc)} bp", file=sys.stderr)

# Splice exons from the revcomp locus
cds_parts = []
for s, e in exons_rc:
    part = ref_locus_rc[s:e+1]
    cds_parts.append(part)
    print(f"  Exon {s}..{e}: {len(part)} bp", file=sys.stderr)

cds = "".join(cds_parts)
prot = str(Seq(cds).translate()).rstrip("*")

print(f"CDS: {len(cds)} bp, Protein: {len(prot)} aa", file=sys.stderr)
print(f"Start: {cds[:3]}, End: {cds[-3:]}", file=sys.stderr)
print(f"Sequence: {prot[:40]}...{prot[-20:]}", file=sys.stderr)

if len(prot) == 0:
    print("ERROR: empty protein", file=sys.stderr)
    sys.exit(1)

# Validate expected values
if len(prot) != 241 or len(cds) != 726:
    print(f"WARNING: Expected 241 aa / 726 bp, got {len(prot)} aa / {len(cds)} bp", file=sys.stderr)

# Save FASTA
with open("crocidura_gueldenstaedtii_procr_protein.fa", "w") as f:
    f.write(f">Crocidura_gueldenstaedtii_PROCR {len(prot)}aa\n")
    for i in range(0, len(prot), 60):
        f.write(prot[i:i+60] + "\n")

with open("crocidura_gueldenstaedtii_procr_cds.fa", "w") as f:
    f.write(f">Crocidura_gueldenstaedtii_PROCR_CDS {len(cds)}bp {cds[:3]}...{cds[-3:]}\n")
    for i in range(0, len(cds), 60):
        f.write(cds[i:i+60] + "\n")

print(f"✓ Saved protein and CDS", file=sys.stderr)
PYEOF

echo "=== Complete ===" >&2
echo "Output:"
ls -lh crocidura_gueldenstaedtii_procr_*.fa
