# Methods

## Problem

Three species need an EPCR coding sequence and none of them has an annotated *PROCR* record.
They fail in two different ways, which is why the pipeline has two routes.

**Unannotated assembly** — *Blarina brevicauda* (`GCA_051295465.2`, 419 scaffolds, 2.25 Gb) and
*Crocidura gueldenstaedtii* (`GCA_965212715.2`, 145 scaffolds, 2.57 Gb) both have genomes but
carry **no NCBI annotation of any kind**. Annotation is not required: miniprot aligns a relative's
protein onto raw assembly sequence and infers the exon structure itself. This distinction matters
because "no annotated PROCR" is often read as "nothing can be extracted" — the gene is in the
assembly, it simply has no gene model pointing at it.

**No assembly at all** — *Rattus flavipectus* is an NCBI **synonym of *Rattus tanezumi***
(taxid 35732). Queried on 2026-10-12, that taxon has no genome assembly, no WGS contigs and no
PROCR mRNA; all 2,688 of its nucleotide records are shorter than 20 kb (mitochondrial genomes and
barcode amplicons), and the only genome-scale data are 9 WGS runs in the SRA. Nothing deposited can
be handed to miniprot, so the locus is reconstructed from reads first.

Every sequence produced here is a gene prediction and should be cited with its assembly or run
accession, never as a deposited sequence.

## Assembly route

```
miniprot -t 8 -j 2 --gff --trans --aln --outn=5 <assembly.fna> <query.faa>
```

The `Rank=1` hit defines the contig, strand and CDS exons; `scripts/extract_cds_from_gff.py`
splices them (reverse-complementing per exon for a minus-strand gene), translates, and runs the
checks below. Queries are chosen by clade: *R. norvegicus* `NP_001020904.1` for the rats,
*Sorex araneus* `XP_004612598.1` for the shrews.

### The start-codon correction

miniprot reports the alignment of the *query*, which is not always the gene's own reading frame.
With the end bonus (`-B`, default 5) it prefers to align the full query, so a diverged query's
N-terminus can run past the target's start codon into upstream genomic sequence. On *Crocidura*
(82% identity to the *Sorex* query) it extended 15 bp / 5 codons, and the raw translation began
`NVSPR`-`MLTALLP…`.

The correction is deterministic: if the translation does not begin with M, trim to the first
in-frame ATG within `--max-trim-codons` (default 20) and adjust the exon coordinates to match. On
*Crocidura* this trims 15 bp (codon 6) and gives a 726 bp CDS, ATG…TAG, 241 aa; on *Blarina*, whose
translation already starts with M, it trims nothing. Earlier work on *Crocidura* achieved the same
result with hardcoded exon coordinates — `scripts/standalone/README.md` shows that those
coordinates are miniprot's own exons with exon 1 trimmed by exactly those 15 bp, and why encoding
the rule beats freezing the coordinates.

### Checks applied to every assembly-route extraction

Recorded in `results/<taxon>/qc.tsv`, and the script exits non-zero if any fails: CDS is a multiple
of three, translation starts with M, CDS ends in a stop codon, no internal stop, no ambiguous
residue, plus the reported frameshift count (`fs:i:` in miniprot's PAF) and the exon count. These
are exactly the checks that catch an N-terminal over-extension, a missed exon or a frameshifted
assembly region.

## The scaffold species

*Rattus rattus* (`GCF_011064425.1`, Rrattus_CSIRO_v1, 2,382,787,281 bp in 2,173 sequences) is run
through the assembly route like any other species, but it is not a target. It has a RefSeq
annotation, so its extraction is an end-to-end control: miniprot's translation must equal
`XP_032760446.1`, and it does. It also supplies the bait locus (gene ± 25 kb) for the reads route.

`-j 2` selects the vertebrate/insect splice model; `--trans` emits the translated CDS as a `##STA`
record.

## Reads route: reconstruction for a taxon with no assembly

### Getting the reads

`prefetch`/`fasterq-dump` resolve accessions through an NCBI name-resolver service. Where that
service is unreachable, both fail with `Failed to call external services` — even when handed a local
file path. The pipeline therefore queries the SDL endpoint itself for the delivery URL, size and
md5, pulls the object over ranged parallel HTTP, verifies the md5 against NCBI's manifest, and runs
sra-tools offline via an `NCBI_SETTINGS` file containing:

```
/repository/remote/disabled = "true"
/repository/user/cache-disabled = "true"
```

### Two-pass read selection

This is the step that determines whether the result is real.

**Pass 1 — capture.** Every read of the run is streamed
(`fasterq-dump --split-spot --skip-technical --stdout`) into `minimap2 -ax sr` against the bait
(gene ±25 kb, ~53.7 kb) and mapped reads are kept. Reads are mapped single-end; pairing is not
needed because the consensus is called per base.

**Pass 1 is not specific, and MAPQ does not fix it.** A ~50 kb reference contains no competing copy
of the repeat families inside the introns, so reads from homologous repeats elsewhere in the genome
have nowhere better to go and are assigned high MAPQ. Measured on SRR3948161: mean depth over PROCR
was **825× at MAPQ ≥ 20 against an expected ~14×**, i.e. roughly 98% of retained reads were
off-target. A consensus called at this stage would be contaminated by paralogous repeat reads.

**Pass 2 — competitive re-mapping.** The pass-1 survivors are re-mapped against the **whole 2.4 Gb
assembly** (`minimap2 -ax sr -t 8 -I 8G`), and only alignments with MAPQ ≥ 20 that land back inside
the locus are kept. Depth then falls to the run's genome-wide coverage, which is the check that the
capture is specific. Pass 1 is a sensitivity filter; pass 2 supplies the specificity.

### Consensus and extraction

```
samtools consensus -a [-A] --mode simple --min-depth 3 --call-fract 0.6 \
                   --min-BQ 15 --min-MQ 20 -r <locus>
```

`--min-depth 3` writes N where fewer than 3 reads cover a base, so unsupported positions cannot
silently inherit the *R. rattus* reference base. miniprot is then run on the consensus with the same
parameters as for the assembly, so the gene is re-found de novo and no coordinate lift-over is
involved (the consensus differs in length from the reference window because of indels).

**The base-quality floor is load-bearing.** `samtools consensus` applies no `--min-BQ` by default,
and with `-A` it then promotes low-quality base noise to IUPAC codes: on SRR2917427 it flagged six
CDS positions as heterozygous for which `samtools mpileup -q 20 -Q 15` reports *zero* mismatching
bases, all six clustered in the deepest exon. Running the consensus with the same floors as the
variant caller (`--min-BQ 15 --min-MQ 20`) removes all six and leaves the single site that is
genuinely split 19/41 reads.

**Two consensus flavours.** The pipeline calls the locus twice from the same BAM:
`locus_consensus.fa` (majority base) and `locus_consensus_ambig.fa` (`-A`, IUPAC at heterozygous
sites). The protein is *not* translated from the majority call alone: at a 50/50 heterozygous site no
base reaches `--call-fract`, so the majority call is N and miniprot emits X. Instead the CDS is
translated ambiguity-aware — a codon contributes its single residue when every resolution of the
IUPAC code agrees (a silent heterozygous site) and X only when they disagree.

The two flavours are **not coordinate-compatible** (`-A` retains ambiguous insertions, so they differ
by ~150-200 bp over 53 kb). Each therefore gets its own miniprot pass, and step 04 refuses to compare
them if the lifted CDS lengths disagree.

### Cross-checks

1. The *R. rattus* miniprot translation is compared with RefSeq `XP_032760446.1`.
2. `bcftools mpileup | bcftools call -mv` is run over the locus independently of the consensus; the
   exonic variants it reports must match, position for position, the CDS differences found by direct
   sequence comparison.
3. Locus-wide divergence (variants per kb across the full bait) is reported as evidence that the
   consensus tracks the reads rather than reverting to the reference — a CDS that comes out identical
   to the reference is only believable if the surrounding introns are visibly divergent.
4. Per-exon depth and the number of bases below `--min-depth` are tabulated, so coverage support for
   the CDS is explicit.
5. Two independent WGS runs of the same taxon (SRR3948161 ~12×, SRR2917427 ~32×, different
   individuals) are reconstructed separately and compared. This is the check that matters most for a
   read-based reconstruction: a reference-bias or mapping artefact would not be expected to
   reproduce, whereas the four synonymous differences from *R. rattus* recovered identically in both
   runs, at two and a half times different depth, are a property of the taxon.

## Reference bias — what this design does and does not control

Reference bias is the main threat to a mapping-based reconstruction. The design controls it at
positions where reads mapped: with ≥3 reads and a 0.6 calling fraction, a true substitution carried
by the reads is called, not the reference base. It does **not** control positions where no read of
the target taxon maps at all — e.g. a *flavipectus*-specific insertion, or a region too divergent
for `minimap2 -ax sr`. Those positions appear as N (below `--min-depth`) rather than as reference
bases, which makes the gap visible rather than silent, and the per-exon depth table reports how many
such bases fall in CDS. A large structural difference in the gene would not be recoverable this way
and would need a de novo assembly of the run.

## Alignment

MAFFT L-INS-i (`--localpair --maxiterate 1000`) on proteins and on CDS. Pairwise identity, the
synonymous/nonsynonymous split and per-codon annotation are computed from the aligned sequences with
the standard genetic code.

## Software

miniprot 0.18-r281; minimap2; samtools/bcftools 1.x; sra-tools 3.4.1; MAFFT 7.526; seqkit.

## Data

- `GCA_051295465.2` — *Blarina brevicauda* assembly (mBlaBre1_p1.2), unannotated
- `GCA_965212715.2` — *Crocidura gueldenstaedtii* assembly (mCroZim1.2), unannotated
- `GCF_011064425.1` — *R. rattus* assembly (Rrattus_CSIRO_v1), scaffold and control
- `XP_004612598.1` — *Sorex araneus* EPCR, query for the shrews (245 aa)
- `NP_001020904.1` — *R. norvegicus* EPCR, query for the rats (241 aa)
- `XP_032760446.1` — *R. rattus* RefSeq EPCR, used only as the control
- `SRR3948161` (PRJNA330658, ~12×), `SRR2917427` (PRJNA301523, ~32×) — *R. tanezumi* WGS runs
