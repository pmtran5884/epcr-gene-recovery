# Configuration for the PROCR/EPCR miniprot extraction + locus reconstruction pipeline.
# Every value a user might reasonably want to change lives here.

# --- Query protein (protein-to-genome query for miniprot) -------------------
QUERY_ACC="NP_001020904.1"            # Rattus norvegicus EPCR precursor, 241 aa

# --- Reference assembly (the taxon that HAS a genome) ----------------------
ASM_ACC="GCF_011064425.1"
ASM_NAME="Rrattus_CSIRO_v1"
ASM_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/064/425/GCF_011064425.1_Rrattus_CSIRO_v1/GCF_011064425.1_Rrattus_CSIRO_v1_genomic.fna.gz"
REF_LABEL="rattus_rattus"

# --- Reconstruction target (the taxon that has only reads) ----------------
# Rattus flavipectus is an NCBI synonym of Rattus tanezumi (taxid 35732):
# no assembly, no WGS contigs, no PROCR mRNA. These are two of its WGS runs.
RUNS="SRR3948161 SRR2917427"
TARGET_LABEL="rattus_flavipectus"

# --- Parameters ------------------------------------------------------------
FLANK=25000          # bp of flanking sequence each side of the gene in the bait
THREADS=8            # miniprot / minimap2 / samtools threads  (miniprot -t)
DL_STREAMS=8         # parallel ranged HTTP streams for the .sra download
MIN_MAPQ=20          # MAPQ floor after competitive whole-genome re-mapping
MIN_DEPTH=3          # samtools consensus --min-depth
CALL_FRACT=0.6       # samtools consensus --call-fract
SPLICE_MODEL=2       # miniprot -j : 2 = vertebrate/insect splice model

# --- Layout ----------------------------------------------------------------
WORK="${WORK:-work}"          # large intermediates (genome, .sra, BAMs)
RESULTS="${RESULTS:-results}" # small, publishable outputs
GENE_LABEL="PROCR"
