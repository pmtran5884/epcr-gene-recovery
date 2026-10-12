# Configuration for the EPCR recovery pipeline.
#
# Species, assemblies and query proteins live in config/targets.tsv.
# This file holds the parameters.

# --- Targets (see config/targets.tsv for accessions and URLs) --------------
SCAFFOLD_LABEL="rattus_rattus"        # mode 'scaffold': assembly + annotation
GENOME_TARGETS="blarina_brevicauda crocidura_gueldenstaedtii"   # mode 'genome'
READ_TARGET="rattus_flavipectus"      # mode 'reads': no assembly at all
RUNS="SRR3948161 SRR2917427"          # WGS runs for the reads target

# Backwards-compatible aliases used by steps 03 and 04.
REF_LABEL="$SCAFFOLD_LABEL"
TARGET_LABEL="$READ_TARGET"
QUERY_ACC="NP_001020904.1"            # R. norvegicus EPCR, query for the Rattus clade
ASM_ACC="GCF_011064425.1"             # scaffold assembly (Rrattus_CSIRO_v1)

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
