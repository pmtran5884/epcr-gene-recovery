#!/usr/bin/env python3
"""Turn a miniprot GFF3 hit into a CDS, a protein and a QC record.

miniprot reports the alignment of the *query* protein, which is not always the
gene's own reading frame: with the end bonus (-B) it will happily align the
query's N-terminus into genomic sequence upstream of the real start codon when
the query is diverged. On Crocidura gueldenstaedtii (82% identity to the Sorex
araneus query) it extends 15 bp / 5 codons past the ATG, so its raw translation
begins NVSPR-MLTALLP... rather than MLTALLP...

This script therefore applies one deterministic correction: if the translation
does not begin with M, trim to the first in-frame ATG within --max-trim-codons
and adjust the exon coordinates to match. That reproduces a hand-curated
N-terminal trim without hard-coding any coordinate.

Outputs, in --out-dir:
    protein.faa     trimmed protein, no trailing stop
    cds.fna         the CDS that encodes it
    cds_exons.bed   exons after trimming (0-based half-open, BED convention)
    locus.txt       contig / gene span / strand / identity
    qc.tsv          every check, pass or fail

Standard library only.
"""
import argparse
import os
import sys

BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {b1 + b2 + b3: AAS[i * 16 + j * 4 + k]
         for i, b1 in enumerate(BASES)
         for j, b2 in enumerate(BASES)
         for k, b3 in enumerate(BASES)}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(seq):
    return "".join(CODON.get(seq[i:i + 3].upper(), "X") for i in range(0, len(seq) - 2, 3))


def load_contig(path, want):
    """Stream a single FASTA record out of a multi-GB assembly."""
    chunks, on = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if on:
                    break
                on = line[1:].split()[0] == want
                continue
            if on:
                chunks.append(line.strip())
    if not chunks:
        sys.exit(f"[extract] contig {want} not found in {path}")
    return "".join(chunks).upper()


def parse_gff(path):
    """Rank=1 mRNA and CDS records of the best hit, plus its ##PAF stats."""
    exons, mrna, paf = [], None, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##PAF") and not paf:
                f = line.rstrip("\n").split("\t")
                paf = {"query": f[1], "query_len": f[2], "strand": f[5],
                       "fs": next((t.split(":")[-1] for t in f if t.startswith("fs:i:")), "NA")}
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or "Rank=1" not in f[8]:
                continue
            if f[2] == "mRNA" and mrna is None:
                ident = [a for a in f[8].split(";") if a.startswith("Identity=")]
                mrna = {"contig": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                        "identity": ident[0].split("=")[1] if ident else "NA"}
            elif f[2] == "CDS" and (mrna is None or f[0] == mrna["contig"]):
                exons.append([int(f[3]), int(f[4])])
    if mrna is None or not exons:
        sys.exit(f"[extract] no Rank=1 hit in {path}")
    exons.sort()
    return mrna, exons, paf


def splice(contig_seq, exons, strand):
    """Concatenate exons in translation order."""
    if strand == "+":
        return "".join(contig_seq[s - 1:e] for s, e in exons)
    return "".join(revcomp(contig_seq[s - 1:e]) for s, e in reversed(exons))


def trim_exons(exons, strand, nbases):
    """Drop nbases from the 5' end of the CDS, in genomic coordinates."""
    out = [list(x) for x in exons]
    order = range(len(out)) if strand == "+" else range(len(out) - 1, -1, -1)
    left = nbases
    for i in order:
        s, e = out[i]
        span = e - s + 1
        if left == 0:
            break
        if left >= span:
            out[i] = None
            left -= span
        else:
            if strand == "+":
                out[i][0] = s + left
            else:
                out[i][1] = e - left
            left = 0
    return [x for x in out if x is not None]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--query-acc", default="NA")
    ap.add_argument("--assembly-acc", default="NA")
    ap.add_argument("--max-trim-codons", type=int, default=20,
                    help="how far into the alignment to look for the start codon")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    mrna, exons, paf = parse_gff(a.gff)
    contig = load_contig(a.genome, mrna["contig"])
    strand = mrna["strand"]
    cds = splice(contig, exons, strand)
    prot = translate(cds)

    # --- deterministic N-terminal correction
    trim_bp, trim_reason = 0, "none (already starts with M)"
    if not prot.startswith("M"):
        starts = [i for i in range(min(a.max_trim_codons, len(prot)))
                  if cds[i * 3:i * 3 + 3] == "ATG"]
        if starts:
            trim_bp = starts[0] * 3
            trim_reason = f"trimmed {trim_bp} bp to first in-frame ATG (codon {starts[0] + 1})"
            exons = trim_exons(exons, strand, trim_bp)
            cds = splice(contig, exons, strand)
            prot = translate(cds)
        else:
            trim_reason = f"NO in-frame ATG within {a.max_trim_codons} codons -- left untrimmed"

    prot_nostop = prot.rstrip("*")
    checks = {
        "label": a.label,
        "assembly": a.assembly_acc,
        "query": a.query_acc or paf.get("query", "NA"),
        "contig": mrna["contig"],
        "strand": strand,
        "gene_start": mrna["start"],      # miniprot's alignment span (pre-trim)
        "gene_end": mrna["end"],
        "cds_start": min(s for s, _ in exons),   # after the N-terminal trim
        "cds_end": max(e for _, e in exons),
        "identity_to_query": mrna["identity"],
        "frameshifts_reported": paf.get("fs", "NA"),
        "cds_exons": len(exons),
        "cds_bp": len(cds),
        "cds_multiple_of_3": "yes" if len(cds) % 3 == 0 else "NO",
        "protein_aa": len(prot_nostop),
        "starts_with_M": "yes" if prot_nostop.startswith("M") else "NO",
        "ends_with_stop": "yes" if prot.endswith("*") else "NO",
        "internal_stops": prot_nostop.count("*"),
        "ambiguous_residues": prot_nostop.count("X"),
        "n_terminal_trim": trim_reason,
    }

    with open(os.path.join(a.out_dir, "protein.faa"), "w") as fh:
        fh.write(f">{a.label}_EPCR miniprot/{a.assembly_acc} query={checks['query']} "
                 f"{mrna['contig']}:{checks['cds_start']}-{checks['cds_end']}({strand}) "
                 f"{len(prot_nostop)}aa CDS={len(cds)}bp exons={len(exons)}\n")
        for i in range(0, len(prot_nostop), 60):
            fh.write(prot_nostop[i:i + 60] + "\n")
    with open(os.path.join(a.out_dir, "cds.fna"), "w") as fh:
        fh.write(f">{a.label}_PROCR_CDS {len(cds)}bp\n")
        for i in range(0, len(cds), 60):
            fh.write(cds[i:i + 60] + "\n")
    with open(os.path.join(a.out_dir, "cds_exons.bed"), "w") as fh:
        for s, e in exons:
            fh.write(f"{mrna['contig']}\t{s - 1}\t{e}\t{a.label}_CDS\t0\t{strand}\n")
    with open(os.path.join(a.out_dir, "locus.txt"), "w") as fh:
        for k in ("contig", "gene_start", "gene_end", "cds_start", "cds_end",
                  "strand", "identity_to_query"):
            fh.write(f"{k}\t{checks[k]}\n")
    with open(os.path.join(a.out_dir, "qc.tsv"), "w") as fh:
        for k, v in checks.items():
            fh.write(f"{k}\t{v}\n")

    bad = [k for k in ("cds_multiple_of_3", "starts_with_M", "ends_with_stop") if checks[k] == "NO"]
    if checks["internal_stops"]:
        bad.append("internal_stops")
    print(f"[extract] {a.label}: {checks['protein_aa']} aa, {checks['cds_bp']} bp, "
          f"{checks['cds_exons']} exons, identity {checks['identity_to_query']}, "
          f"{trim_reason}" + (f"  ** FAILED: {bad}" if bad else "  (all checks pass)"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
