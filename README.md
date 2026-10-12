# PROCR / EPCR in *Rattus rattus* and *Rattus flavipectus*

miniprot extraction of the endothelial protein C receptor (EPCR, gene *PROCR*) coding sequence from
the *Rattus rattus* genome, and reconstruction of the same locus for *Rattus flavipectus* — a taxon
that has **no deposited genome** — from two independent WGS runs.

## Result

All three EPCR proteins are **identical: 240 aa, 100% identity, no gaps.** The differences are
entirely synonymous and sit at third codon positions:

| A | B | CDS bp | CDS diffs | CDS identity | synonymous | nonsynonymous | protein identity |
|---|---|---|---|---|---|---|---|
| *R. rattus* | SRR3948161 | 723 | 4 | 99.447% | 4 | 0 | 100% |
| *R. rattus* | SRR2917427 | 723 | 5 | 99.308% | 5 | 0 | 100% |
| SRR3948161 | SRR2917427 | 723 | 1 | 99.862% | 1 | 0 | 100% |

| CDS position | change | codon | residue | present in |
|---|---|---|---|---|
| c.48 | T→C | TGT→TGC | Cys16 | both runs |
| c.66 | T→C | GAT→GAC | Asp22 | both runs |
| c.294 | A→G | TTA→TTG | Leu98 | both runs |
| c.369 | T→G | CCT→CCG | Pro123 | both runs |
| c.426 | T→Y (T/C het) | TTT→TTY | Phe142 | SRR2917427 only |

The four fixed differences are recovered **identically in two different individuals sequenced at
12× and 32×**, which is the evidence that they are properties of the taxon rather than mapping
artefacts. The fifth is a heterozygous silent site private to the deeper individual.

*R. rattus* locus: `NC_046158.1:135,103,022–135,106,686 (+)`, 4 CDS exons, 97.93% identity /
99.17% positives to the *R. norvegicus* query. The translation is **byte-identical to RefSeq
`XP_032760446.1`**, which validates the extraction independently of this pipeline.

## Read this before citing the *flavipectus* sequence

*Rattus flavipectus* is an NCBI **synonym of *Rattus tanezumi*** (taxid 35732), and that taxon has
**no genome assembly, no WGS contigs and no PROCR mRNA** on record — as queried on 2026-10-12, all
2,688 of its nucleotide records are shorter than 20 kb (mitochondrial genomes, barcode amplicons).
There is nothing deposited for miniprot to align a protein against. The sequences here are
**read-based reconstructions** of the locus from raw WGS reads and must be cited as such, with the
run accession — not as deposited sequences.

Coverage support for the CDS (no base below the `--min-depth 3` floor in either run):

| run | gene depth | exon 1 | exon 2 | exon 3 | exon 4 | locus divergence vs *R. rattus* |
|---|---|---|---|---|---|---|
| SRR3948161 | 12.46× | 7.6× | 11.0× | 16.3× | 10.6× | 961 variants / 53.7 kb (819 SNV, 142 indel) |
| SRR2917427 | 32.21× | 32.0× | 32.4× | 45.8× | 29.6× | 1,359 variants / 53.7 kb (1,168 SNV, 191 indel) |

The locus-wide divergence matters: a CDS that comes out nearly identical to the reference is only
believable if the surrounding introns are visibly divergent, i.e. if the consensus is demonstrably
tracking the reads rather than reverting to the reference.

## Quickstart

```bash
conda env create -f environment.yml
conda activate procr-miniprot
./run_all.sh            # ~2 h per WGS run, dominated by a 15-39 GB download
```

Everything configurable is in [`config/params.sh`](config/params.sh) — query protein, assembly,
run accessions, flank size, thread count, MAPQ/depth thresholds. Point `RUNS` at other SRA
accessions, or `QUERY_ACC`/`ASM_URL` at another gene or assembly, and nothing else needs editing.

Disk: ~45 GB for the assembly plus the size of each `.sra` (15 GB and 39 GB for the two runs here).
Peak RAM ~10 GB (miniprot indexing the 2.4 Gb genome).

## Pipeline

| step | script | what it does |
|---|---|---|
| 00 | [`scripts/00_fetch_inputs.sh`](scripts/00_fetch_inputs.sh) | query protein from NCBI, reference assembly |
| 01 | [`scripts/01_extract_reference.sh`](scripts/01_extract_reference.sh) | `miniprot -t 8 -j 2` on the assembly; bait = gene ± 25 kb |
| 02 | [`scripts/02_fetch_sra.sh`](scripts/02_fetch_sra.sh) | `.sra` via the NCBI SDL manifest, md5-verified |
| 03 | [`scripts/03_reconstruct_locus.sh`](scripts/03_reconstruct_locus.sh) | two-pass read capture → consensus → miniprot → variants |
| 04 | [`scripts/04_align_and_report.py`](scripts/04_align_and_report.py) | MAFFT alignments, pairwise tables, `results/report.md` |

### The one thing to understand before reusing this

Capturing reads directly onto a small extracted locus **is not specific, and MAPQ does not fix it.**
A ~50 kb bait contains no competing copy of the repeat families in its introns, so reads from
homologous repeats elsewhere in the genome have nowhere better to go and are handed high MAPQ.
Measured here: **825× over PROCR at MAPQ ≥ 20 against an expected ~14×** — roughly 98% of retained
reads were off-target. A consensus called at that stage is contaminated by paralogous reads.

So pass 1 (capture against the bait) is used only as a sensitivity filter, and its survivors are
re-mapped **competitively against the whole 2.4 Gb assembly**, keeping MAPQ ≥ 20 alignments that land
back in the locus. Depth then drops to the run's genome-wide coverage, which is the signal that the
capture is now specific.

Two further traps, both documented in [`docs/METHODS.md`](docs/METHODS.md): `samtools consensus`
applies no base-quality floor by default (with `-A` it promoted six positions of noise to
heterozygous calls that `mpileup -Q 15` shows zero support for), and the majority-base and `-A`
consensus flavours are not coordinate-compatible, so each needs its own miniprot pass.

## Outputs

```
results/
  summary.tsv                  one row per pairwise comparison
  codon_differences.tsv        one row per nucleotide difference, with effect
  heterozygous_cds_sites.tsv   het CDS sites per run
  report.md                    the tables above, rendered
  proteins_all.faa  cds_all.fna
  alignments/                  MAFFT L-INS-i protein and CDS alignments
  rattus_rattus/               miniprot.gff, protein.faa, bait.fa, cds_exons.bed, locus.txt
  <SRR accession>/             protein.faa, locus_consensus.fa, miniprot.gff, locus.bam,
                               locus.vcf.gz, qc.tsv, exon_depth.tsv, sra_manifest.txt, logs
```

## Data and software

Data: `GCF_011064425.1` (*R. rattus* assembly) · `NP_001020904.1` (*R. norvegicus* query) ·
`XP_032760446.1` (RefSeq check) · `SRR3948161`, `SRR2917427` (*R. tanezumi* WGS).

Software: miniprot 0.18-r281 · minimap2 · samtools/bcftools 1.x · sra-tools 3.4.1 · MAFFT 7.526.

Please cite [miniprot](https://github.com/lh3/miniprot) (Li 2023, *Bioinformatics*) and
[minimap2](https://github.com/lh3/minimap2) (Li 2018) if you use this.

## License

MIT — see [LICENSE](LICENSE).
