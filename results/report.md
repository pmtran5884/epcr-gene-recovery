# PROCR / EPCR extraction and alignment

Query protein: `NP_001020904.1`  |  reference assembly: `GCF_011064425.1` (`NC_046158.1:135103022-135106686`, identity to query 0.9793)

## Sequences

| label | source | protein length |
|---|---|---|
| rattus_rattus | GCF_011064425.1 assembly, miniprot | 240 aa |
| SRR3948161 | WGS reads, reconstruction (12.46x over gene) | 240 aa |
| SRR2917427 | WGS reads, reconstruction (32.21x over gene) | 240 aa |

## Pairwise comparisons

| A | B | CDS bp | CDS diffs | CDS identity | syn | nonsyn | protein diffs | protein identity |
|---|---|---|---|---|---|---|---|---|
| rattus_rattus | SRR3948161 | 723 | 4 | 99.447% | 4 | 0 | 0 | 100.000% |
| rattus_rattus | SRR2917427 | 723 | 5 | 99.308% | 5 | 0 | 0 | 100.000% |
| SRR3948161 | SRR2917427 | 723 | 1 | 99.862% | 1 | 0 | 0 | 100.000% |

## Nucleotide differences

| A | B | CDS pos | change | codon | residue | effect |
|---|---|---|---|---|---|---|
| rattus_rattus | SRR3948161 | c.48 | T>C | TGT>TGC | C16C | synonymous |
| rattus_rattus | SRR3948161 | c.66 | T>C | GAT>GAC | D22D | synonymous |
| rattus_rattus | SRR3948161 | c.294 | A>G | TTA>TTG | L98L | synonymous |
| rattus_rattus | SRR3948161 | c.369 | T>G | CCT>CCG | P123P | synonymous |
| rattus_rattus | SRR2917427 | c.48 | T>C | TGT>TGC | C16C | synonymous |
| rattus_rattus | SRR2917427 | c.66 | T>C | GAT>GAC | D22D | synonymous |
| rattus_rattus | SRR2917427 | c.294 | A>G | TTA>TTG | L98L | synonymous |
| rattus_rattus | SRR2917427 | c.369 | T>G | CCT>CCG | P123P | synonymous |
| rattus_rattus | SRR2917427 | c.426 | T>Y | TTT>TTY | F142F | synonymous_het |
| SRR3948161 | SRR2917427 | c.426 | T>Y | TTT>TTY | F142F | synonymous_het |

## Heterozygous CDS sites in the reconstructions

CDS positions where that individual's reads carry two alleles. The protein is
translated ambiguity-aware: a site like this contributes its single residue when
every resolution of the codon agrees, and an X only when they disagree.

| run | CDS pos | IUPAC | codon | residue(s) | effect |
|---|---|---|---|---|---|
| SRR2917427 | c.426 | Y | TTY | F | silent |

## Read-reconstruction QC

**SRR3948161** — pass-1 candidates 2443968, locus reads after competitive re-mapping 4252, mean gene depth 12.46x, locus variants 961 (819 SNVs / 142 indels), gene-span 35, CDS-exonic 4.

| exon | start | end | bp | mean depth | min depth | bp below min-depth |
|---|---|---|---|---|---|---|
| 1 | 135103022 | 135103088 | 67 | 7.6 | 6 | 0 |
| 2 | 135104974 | 135105231 | 258 | 11.0 | 7 | 0 |
| 3 | 135106059 | 135106340 | 282 | 16.3 | 12 | 0 |
| 4 | 135106571 | 135106686 | 116 | 10.6 | 9 | 0 |

**SRR2917427** — pass-1 candidates 5214020, locus reads after competitive re-mapping 10985, mean gene depth 32.21x, locus variants 1359 (1168 SNVs / 191 indels), gene-span 55, CDS-exonic 5.

| exon | start | end | bp | mean depth | min depth | bp below min-depth |
|---|---|---|---|---|---|---|
| 1 | 135103022 | 135103088 | 67 | 32.0 | 28 | 0 |
| 2 | 135104974 | 135105231 | 258 | 32.4 | 28 | 0 |
| 3 | 135106059 | 135106340 | 282 | 45.8 | 39 | 0 |
| 4 | 135106571 | 135106686 | 116 | 29.6 | 25 | 0 |

