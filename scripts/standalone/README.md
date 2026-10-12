# Standalone recovery scripts (originals, kept verbatim)

These are the two single-file scripts that first recovered EPCR from the shrew assemblies.
They are self-contained — each downloads its own genome, installs miniprot if missing, and
writes a FASTA — and are kept here unmodified as the record of how the sequences were first
obtained. The integrated route (`scripts/01_extract_from_assembly.sh` +
`scripts/extract_cds_from_gff.py`) reproduces both results; see below.

| script | assembly | result | reproduced by the pipeline? |
|---|---|---|---|
| `recover_blarina_procr.sh` | GCA_051295465.2 (mBlaBre1_p1.2) | 242 aa, `JBKFBT020000153.1:87596-91813 (+)` | yes, byte-identical |
| `recover_crocidura_gueldenstaedtii_procr.sh` | GCA_965212715.2 (mCroZim1.2) | 241 aa / 726 bp, `OZ243422.1:4125072-4130097 (−)` | yes, byte-identical (protein and CDS) |

Both use the same query, *Sorex araneus* EPCR — pasted inline in the scripts, and identical to
RefSeq **XP_004612598.1** (245 aa), which is what the pipeline fetches by accession.

## What the Crocidura script is actually doing

A reader can easily conclude from this script that miniprot failed for *Crocidura* and was
replaced by hand-built coordinates. That is not what happened. miniprot found the gene correctly;
three of its four CDS exons are **identical** to the script's hardcoded values:

| exon | hardcoded (via the revcomp window) | miniprot (`Rank=1`) | agree? |
|---|---|---|---|
| 1 | OZ243422.1:4130022-4130097 | OZ243422.1:4130022-**4130112** | start yes, end 15 bp apart |
| 2 | OZ243422.1:4127549-4127800 | OZ243422.1:4127549-4127800 | identical |
| 3 | OZ243422.1:4126117-4126395 | OZ243422.1:4126117-4126395 | identical |
| 4 | OZ243422.1:4125072-4125190 | OZ243422.1:4125072-4125190 | identical |

The single disagreement is the 5′ end of exon 1 (the gene is on the minus strand), and it is
exactly 15 bp = 5 codons. miniprot's raw translation therefore begins `NVSPR`-`MLTALLP…` instead
of `MLTALLP…`: with the end bonus (`-B`, default 5) it prefers to align the full query, so the
N-terminus of the 82%-identical *Sorex* query runs past the *Crocidura* start codon into upstream
genomic sequence. The hardcoded exon 1 is miniprot's exon 1 with those 15 bp trimmed off — a
hand-curated N-terminal trim, not a workaround for a failed alignment.

So the fix is real and correct; only the *way* it is encoded is fragile:

1. **The coordinates are frozen to one assembly version.** `GCA_965212715.2` is already a `.2`
   release; any re-release shifts the window and the script would splice the wrong bases.
2. **The scaffold is assumed, not checked.** The script reads the top hit's scaffold out of the PAF
   and then applies a window hardcoded for `OZ243422.1` to *whatever* that scaffold is. If the best
   hit ever landed elsewhere, it would silently splice unrelated sequence.
3. **The validation warns instead of failing** (`WARNING: Expected 241 aa…`) and never checks the
   things that would have caught the N-terminal problem in the first place: starts with ATG,
   no internal stop, CDS a multiple of three.
4. **The comment and the code disagree about memory.** "Index the FASTA to avoid loading the full
   assembly into memory" sits above a loop that loads every scaffold of a 2.4 Gb genome into a dict
   and prints a line per scaffold.
5. `ftp://` with `wget` is fragile on networks that block FTP; `https://` works.
6. Biopython is pulled in for a reverse-complement and a translation.

`extract_cds_from_gff.py` replaces the frozen coordinates with the rule they encode: take
miniprot's exon structure, and if the translation does not start with M, trim to the first in-frame
ATG and adjust the exon coordinates accordingly. On *Crocidura* that trims 15 bp (codon 6) and
yields the same 726 bp CDS (ATG…TAG) and the same 241 aa protein; on *Blarina*, whose translation
already starts with M, it trims nothing. The checks in point 3 run on every species and fail loudly.
