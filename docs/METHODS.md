# Methods

## Problem

*Rattus rattus* has a chromosome-level assembly, so extracting its PROCR/EPCR coding sequence is a
protein-to-genome alignment and nothing more. *Rattus flavipectus* does not: the name is an NCBI
**synonym of *Rattus tanezumi*** (taxid 35732), and that taxon has **no genome assembly, no WGS
contigs and no PROCR mRNA** deposited. Queried on 2026-10-12, its 2,688 nucleotide records are all
shorter than 20 kb (mitochondrial genomes and barcode amplicons); the only genome-scale data are
9 WGS runs in the SRA. There is therefore no deposited sequence for miniprot to align a protein
against, and the *flavipectus* sequence reported here is a **read-based reconstruction** of the
locus — it should be cited as such, with the run accession, not as a deposited sequence.

## Query and reference

The query protein is *R. norvegicus* EPCR precursor `NP_001020904.1` (241 aa). The reference
assembly is `GCF_011064425.1` (Rrattus_CSIRO_v1, 2,382,787,281 bp in 2,173 sequences).

miniprot 0.18-r281 is run with the vertebrate splice model and 8 threads:

```
miniprot -t 8 -j 2 --gff --trans --aln --outn=5 <genome.fna> <query.faa>
```

`-j 2` selects the vertebrate/insect splice model; `--trans` emits the translated CDS as a `##STA`
record, which is what we take as the protein. The best hit (`Rank=1`) defines the gene span and the
four CDS exons used downstream.

## Reconstruction from reads

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

- `GCF_011064425.1` — *R. rattus* assembly (Rrattus_CSIRO_v1)
- `NP_001020904.1` — *R. norvegicus* EPCR query protein
- `XP_032760446.1` — *R. rattus* RefSeq EPCR, used only as an independent check
- `SRR3948161`, `SRR2917427` — *R. tanezumi* (= *flavipectus*) WGS runs
