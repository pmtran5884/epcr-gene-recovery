# EPCR recovery in species with no annotated *PROCR*

Recovering the endothelial protein C receptor (EPCR, gene *PROCR*) coding sequence for three
mammals that have **no annotated *PROCR* record** in NCBI or UniProt — two from unannotated
genome assemblies, one from raw WGS reads because it has no assembly at all.

| target species | what is deposited | route | result |
|---|---|---|---|
| ***Blarina brevicauda*** (northern short-tailed shrew) | assembly `GCA_051295465.2`, **no annotation** | miniprot on the assembly | **242 aa**, `JBKFBT020000153.1:87,596-91,813 (+)`, 4 exons, 729 bp |
| ***Crocidura gueldenstaedtii*** (Güldenstädt's shrew) | assembly `GCA_965212715.2`, **no annotation** | miniprot + start-codon trim | **241 aa**, `OZ243422.1:4,125,072-4,130,097 (−)`, 4 exons, 726 bp |
| ***Rattus flavipectus*** (= *R. tanezumi*) | **no assembly at all** — only mito/barcode records <20 kb, plus WGS reads | two-pass read capture → consensus → miniprot | **240 aa**, reconstructed independently from two WGS runs |

*Rattus rattus* (`GCF_011064425.1`) is in the pipeline but is **not** a target: it has both an
assembly and a RefSeq annotation, and serves two purposes — the mapping scaffold for the
*flavipectus* reads, and an end-to-end control, because miniprot's translation of its locus must
come out byte-identical to RefSeq `XP_032760446.1`. It does.

Every protein here is a **gene prediction**, not a deposited sequence. Cite them as such, with the
assembly accession (shrews) or the run accession (*flavipectus*).

## Results

Pairwise protein identity, from the MAFFT alignment of all recovered sequences:

| | *R. rattus* | *B. brevicauda* | *C. gueldenstaedtii* |
|---|---|---|---|
| *R. rattus* | — | 64.98% | 62.87% |
| *B. brevicauda* | 64.98% | — | 82.08% |
| *C. gueldenstaedtii* | 62.87% | 82.08% | — |
| *R. flavipectus* (both runs) | **100%** | 64.98% | 62.87% |

The two shrews sit at 82% to each other and ~63-65% to the rats, consistent with EPCR's fast-
evolving MHC-like fold. Within *Rattus*, the *flavipectus* reconstructions are **identical** to
*R. rattus* at the protein level; their 723 bp CDS differ at 4-5 positions, every one synonymous:

| CDS position | change | codon | residue | seen in |
|---|---|---|---|---|
| c.48 | T→C | TGT→TGC | Cys16 | both runs |
| c.66 | T→C | GAT→GAC | Asp22 | both runs |
| c.294 | A→G | TTA→TTG | Leu98 | both runs |
| c.369 | T→G | CCT→CCG | Pro123 | both runs |
| c.426 | T→Y (T/C het) | TTT→TTY | Phe142 | SRR2917427 only |

The four fixed differences were recovered identically from two individuals in two different
studies (PRJNA330658 at ~12×, PRJNA301523 at ~32×), which is what distinguishes a property of the
taxon from a mapping artefact.

## Two routes, and when each applies

**Assembly route** (`scripts/01_extract_from_assembly.sh`) — for a genome with no annotation.
miniprot aligns a relative's EPCR protein onto the assembly and infers the exon structure itself;
no GFF, GTF or gene model is required. Query proteins: *R. norvegicus* `NP_001020904.1` for the
rats, *Sorex araneus* `XP_004612598.1` for the shrews.

**Reads route** (`scripts/02_fetch_sra.sh`, `scripts/03_reconstruct_locus.sh`) — for a taxon with
no assembly. Reads are captured onto the scaffold species' locus, re-mapped competitively against
its whole genome, and turned into a consensus that miniprot then treats as a genome. This is
reference-guided reconstruction, so it can only report positions where reads mapped; unsupported
positions are written as N rather than inheriting the scaffold base.

## Four traps this pipeline exists to handle

1. **Bait-only read capture is not specific, and MAPQ does not fix it.** A ~50 kb bait holds no
   competing copy of the repeat families in its introns, so repeat reads from elsewhere in the
   genome have nowhere better to go and are handed high MAPQ: 825× over PROCR against an expected
   ~14×, i.e. ~98% off-target. Hence the competitive second pass against the whole 2.4 Gb assembly.
2. **miniprot can over-extend past the start codon.** With a diverged query (82% for *Crocidura*)
   the end bonus makes it align the query's N-terminus into upstream genomic sequence, so the raw
   translation began `NVSPR`-`MLTALLP…`. The fix is a deterministic trim to the first in-frame ATG —
   15 bp here — not hand-written coordinates; see [`scripts/standalone/README.md`](scripts/standalone/README.md).
3. **`samtools consensus` has no base-quality floor by default.** With `-A` it promoted low-quality
   base noise to six "heterozygous" CDS calls in SRR2917427 that `mpileup -Q 15` shows zero support
   for. `--min-BQ 15 --min-MQ 20` removes all six and leaves the one real 19/41 site.
4. **The majority-base and `-A` consensus flavours are not coordinate-compatible** (`-A` keeps
   ambiguous insertions; they differ by ~150-200 bp per 53 kb), so each gets its own miniprot pass
   and the report step refuses to compare them if the lifted CDS lengths disagree.

## Quickstart

```bash
conda env create -f environment.yml
conda activate procr-miniprot
./run_all.sh
```

Runtime: the assembly route is ~2 min per species (download + miniprot ~20 s, peak ~20 GB RAM for
the 2.4 Gb genomes). The reads route is ~1-2 h per run, dominated by the 15-39 GB `.sra` download
and one streaming pass over all reads.

To add a species, append a line to [`config/targets.tsv`](config/targets.tsv) — label, mode,
assembly accession, query accession, URL — and nothing else needs editing. Thresholds live in
[`config/params.sh`](config/params.sh).

## Layout

```
config/targets.tsv      species, mode (genome|reads|scaffold), assembly, query protein
config/params.sh        thresholds: flank, threads, MAPQ, min-depth, call-fract
scripts/01_extract_from_assembly.sh   assembly route (any unannotated genome)
scripts/extract_cds_from_gff.py       GFF -> CDS/protein + start-codon trim + QC checks
scripts/02_fetch_sra.sh               .sra via the NCBI SDL manifest, md5-verified
scripts/03_reconstruct_locus.sh       two-pass capture -> consensus -> miniprot -> variants
scripts/04_align_and_report.py        MAFFT alignments, identity tables, results/report.md
scripts/standalone/                   the original single-file shrew scripts, kept verbatim
results/<taxon>/        protein.faa, cds.fna, miniprot.gff, cds_exons.bed, locus.txt, qc.tsv
results/<SRR>/          protein.faa, locus_consensus.fa, locus.bam, locus.vcf.gz, qc.tsv, depth
results/                proteins_all.faa, cds_all.fna, protein_identity.tsv, summary.tsv,
                        codon_differences.tsv, heterozygous_cds_sites.tsv, alignments/, report.md
```

## Data and software

Assemblies `GCA_051295465.2`, `GCA_965212715.2`, `GCF_011064425.1` · query proteins
`XP_004612598.1`, `NP_001020904.1` · control `XP_032760446.1` · WGS runs `SRR3948161`,
`SRR2917427`.

miniprot 0.18-r281 · minimap2 · samtools/bcftools 1.x · sra-tools 3.4.1 · MAFFT 7.526 · seqkit.
Please cite [miniprot](https://github.com/lh3/miniprot) (Li 2023, *Bioinformatics*) and
[minimap2](https://github.com/lh3/minimap2) (Li 2018).

## License

MIT — see [LICENSE](LICENSE).
