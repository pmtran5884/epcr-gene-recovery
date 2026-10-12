# EPCR recovery in species with no annotated PROCR

## Sequences recovered

| taxon | mode | source | query | protein |
|---|---|---|---|---|
| rattus_rattus | scaffold | GCF_011064425.1 `NC_046158.1:135103022-135106686(+)` | NP_001020904.1 | 240 aa |
| blarina_brevicauda | genome | GCA_051295465.2 `JBKFBT020000153.1:87596-91813(+)` | XP_004612598.1 | 242 aa |
| crocidura_gueldenstaedtii | genome | GCA_965212715.2 `OZ243422.1:4125072-4130097(-)` | XP_004612598.1 | 241 aa |
| rattus_flavipectus | reads | SRR3948161, 12.46x over gene | NP_001020904.1 | 240 aa |
| rattus_flavipectus | reads | SRR2917427, 32.21x over gene | NP_001020904.1 | 240 aa |

## Extraction QC (assembly route)

| taxon | identity to query | exons | CDS bp | protein | starts M | internal stops | frameshifts | N-terminal trim |
|---|---|---|---|---|---|---|---|---|
| rattus_rattus | 0.9793 | 4 | 723 | 240 aa | yes | 0 | 0 | none (already starts with M) |
| blarina_brevicauda | 0.8694 | 4 | 729 | 242 aa | yes | 0 | 0 | none (already starts with M) |
| crocidura_gueldenstaedtii | 0.8178 | 4 | 726 | 241 aa | yes | 0 | 0 | trimmed 15 bp to first in-frame ATG (codon 6) |

## Protein identity (all taxa, from the MAFFT alignment)

| A | B | len A | len B | aligned columns | diffs | identity |
|---|---|---|---|---|---|---|
| rattus_rattus | blarina_brevicauda | 240 | 242 | 237 | 83 | 64.98% |
| rattus_rattus | crocidura_gueldenstaedtii | 240 | 241 | 237 | 88 | 62.87% |
| rattus_rattus | SRR3948161 | 240 | 240 | 240 | 0 | 100.00% |
| rattus_rattus | SRR2917427 | 240 | 240 | 240 | 0 | 100.00% |
| blarina_brevicauda | crocidura_gueldenstaedtii | 242 | 241 | 240 | 43 | 82.08% |
| blarina_brevicauda | SRR3948161 | 242 | 240 | 237 | 83 | 64.98% |
| blarina_brevicauda | SRR2917427 | 242 | 240 | 237 | 83 | 64.98% |
| crocidura_gueldenstaedtii | SRR3948161 | 241 | 240 | 237 | 88 | 62.87% |
| crocidura_gueldenstaedtii | SRR2917427 | 241 | 240 | 237 | 88 | 62.87% |
| SRR3948161 | SRR2917427 | 240 | 240 | 240 | 0 | 100.00% |

## CDS comparison within the *Rattus* group

Codon-level comparison is restricted to taxa whose CDS are the same length; across species the lengths differ and the protein table above applies instead.

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

