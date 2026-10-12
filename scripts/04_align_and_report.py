#!/usr/bin/env python3
"""Collect the reference extraction and every run reconstruction, align them,
and write the summary tables.

Inputs (produced by steps 01 and 03), discovered from config/params.sh:
    results/<REF_LABEL>/{protein.faa,miniprot.gff,bait.fa,locus.txt,cds_exons.bed}
    results/<ACC>/{protein.faa,miniprot.gff,locus_consensus.fa,qc.tsv}

Outputs:
    results/proteins_all.faa            every protein, one FASTA
    results/cds_all.fna                 every CDS, one FASTA
    results/alignments/protein_alignment.{fasta,clustal}
    results/alignments/cds_alignment.{fasta,clustal}
    results/summary.tsv                 one row per pairwise comparison
    results/codon_differences.tsv       one row per nucleotide difference
    results/report.md                   human-readable summary

Only the standard library plus MAFFT on PATH.
"""
import os
import re
import subprocess
import sys
from itertools import combinations

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BASES = "TCAG"
AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
CODON = {b1 + b2 + b3: AAS[i * 16 + j * 4 + k]
         for i, b1 in enumerate(BASES)
         for j, b2 in enumerate(BASES)
         for k, b3 in enumerate(BASES)}


def params():
    """Read the handful of shell variables we need out of config/params.sh."""
    out = {}
    with open(os.path.join(ROOT, "config", "params.sh")) as fh:
        for line in fh:
            m = re.match(r'^(\w+)="?([^"#\n]*?)"?\s*(?:#.*)?$', line.strip())
            if not m:
                continue
            val = m.group(2)
            # resolve ${VAR:-default} the way the shell would, honouring the
            # environment first so WORK=/scratch ./scripts/04_...py works too
            def expand(mm):
                return os.environ.get(mm.group(1), mm.group(2))
            val = re.sub(r'\$\{(\w+):-([^}]*)\}', expand, val)
            out[m.group(1)] = val
    return out


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def cds_exons(gff):
    """Rank=1 CDS intervals, 1-based inclusive, in the coordinates of the GFF."""
    ex = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "CDS" and "Rank=1" in f[8]:
                ex.append((int(f[3]), int(f[4])))
    return sorted(ex)


IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
         "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT",
         "V": "ACG", "N": "ACGT"}


def translate(seq):
    return "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def aa_set(codon):
    """Every amino acid a (possibly ambiguous) codon can encode."""
    out = set()
    for b1 in IUPAC.get(codon[0], "ACGT"):
        for b2 in IUPAC.get(codon[1], "ACGT"):
            for b3 in IUPAC.get(codon[2], "ACGT"):
                out.add(CODON[b1 + b2 + b3])
    return out


def classify(codon_a, codon_b):
    """Classify a codon difference, tolerating IUPAC ambiguity codes.

    synonymous      -> both codons can only encode the same single residue
    synonymous_het  -> at least one codon is ambiguous, but every resolution of
                       both encodes that same residue (a heterozygous silent
                       site, which is what a diploid WGS consensus produces)
    nonsynonymous   -> some resolution changes the residue
    """
    sa, sb = aa_set(codon_a), aa_set(codon_b)
    a_lab, b_lab = "/".join(sorted(sa)), "/".join(sorted(sb))
    if sa != sb:
        return "nonsynonymous", a_lab, b_lab
    if len(sa) == 1:
        amb = any(c not in "ACGT" for c in codon_a + codon_b)
        return ("synonymous_het" if amb else "synonymous"), a_lab, b_lab
    return "nonsynonymous", a_lab, b_lab


def mafft(infile, outfile, clustal=False):
    cmd = ["mafft", "--quiet", "--maxiterate", "1000", "--localpair"]
    if clustal:
        cmd.append("--clustalout")
    cmd.append(infile)
    with open(outfile, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)


def pair_stats(a, b):
    diffs = [(i + 1, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y]
    ident = 100.0 * (1 - len(diffs) / len(a)) if a else 0.0
    return diffs, ident


def main():
    P = params()
    res = os.path.join(ROOT, P.get("RESULTS", "results"))
    ref_label = P["REF_LABEL"]
    runs = [r for r in P["RUNS"].split() if os.path.isdir(os.path.join(res, r))]
    os.makedirs(os.path.join(res, "alignments"), exist_ok=True)

    # --- reference CDS, lifted out of the bait using genomic exon coordinates
    refdir = os.path.join(res, ref_label)
    locus = dict(l.split("\t") for l in open(os.path.join(refdir, "locus.txt")).read().splitlines())
    bait = list(read_fasta(os.path.join(refdir, "bait.fa")).values())[0]
    off = int(locus["bait_start"])
    cds = {ref_label: "".join(bait[s - off:e - off + 1] for s, e in cds_exons(os.path.join(refdir, "miniprot.gff")))}
    prot = {ref_label: list(read_fasta(os.path.join(refdir, "protein.faa")).values())[0]}

    # --- one CDS per run, lifted out of that run's own consensus
    het = {}   # acc -> [(cds_pos, iupac_code, codon, residue(s), effect), ...]
    for acc in runs:
        d = os.path.join(res, acc)
        ex = cds_exons(os.path.join(d, "miniprot.gff"))
        con = list(read_fasta(os.path.join(d, "locus_consensus.fa")).values())[0]
        cds[acc] = "".join(con[s - 1:e] for s, e in ex)

        # Heterozygous CDS sites, from the IUPAC flavour of the same call
        # (identical coordinates: same BAM, same depth/indel handling).
        amb_path = os.path.join(d, "locus_consensus_ambig.fa")
        amb_gff = os.path.join(d, "miniprot_ambig.gff")
        het[acc] = []
        cds_amb = cds[acc]
        if os.path.exists(amb_path) and os.path.exists(amb_gff):
            amb = list(read_fasta(amb_path).values())[0]
            ex_amb = cds_exons(amb_gff)
            cds_amb = "".join(amb[s - 1:e] for s, e in ex_amb)
            if len(cds_amb) != len(cds[acc]):
                raise SystemExit(f"[04] {acc}: CDS length mismatch between consensus flavours "
                                 f"({len(cds[acc])} vs {len(cds_amb)}) -- refusing to compare")
            for i, b in enumerate(cds_amb):
                if b not in "ACGT":
                    c = i // 3
                    codon = cds_amb[c * 3:c * 3 + 3]
                    aas = aa_set(codon)
                    het[acc].append((i + 1, b, codon, "/".join(sorted(aas)),
                                     "silent" if len(aas) == 1 else "coding-changing"))

        # Compare on the ambiguity-aware CDS: at a 50/50 het site the majority
        # call is N, which would read as a residue change; the IUPAC code carries
        # both alleles and classify() can see that both encode the same residue.
        cds[acc] = cds_amb

        # The protein of record: translate the ambiguity-aware CDS, so a silent
        # heterozygous site yields its single residue instead of an X, and only a
        # genuinely coding-changing het site stays X.
        prot[acc] = "".join(
            (lambda s: next(iter(s)) if len(s) == 1 else "X")(aa_set(cds_amb[i:i + 3]))
            for i in range(0, len(cds_amb) - 2, 3)).rstrip("*")
        with open(os.path.join(d, "protein.faa"), "w") as fh:
            fh.write(f">{P['TARGET_LABEL']}_EPCR {acc}_reconstruction iupac_aware_translation\n{prot[acc]}\n")

    labels = [ref_label] + runs
    with open(os.path.join(res, "proteins_all.faa"), "w") as fh:
        for k in labels:
            fh.write(f">{k}\n{prot[k]}\n")
    with open(os.path.join(res, "cds_all.fna"), "w") as fh:
        for k in labels:
            fh.write(f">{k}\n{cds[k]}\n")

    if len(labels) > 1:
        mafft(os.path.join(res, "proteins_all.faa"), os.path.join(res, "alignments", "protein_alignment.fasta"))
        mafft(os.path.join(res, "proteins_all.faa"), os.path.join(res, "alignments", "protein_alignment.clustal"), clustal=True)
        mafft(os.path.join(res, "cds_all.fna"), os.path.join(res, "alignments", "cds_alignment.fasta"))
        mafft(os.path.join(res, "cds_all.fna"), os.path.join(res, "alignments", "cds_alignment.clustal"), clustal=True)

    # --- pairwise tables
    rows, codon_rows = [], []
    for a, b in combinations(labels, 2):
        if len(cds[a]) != len(cds[b]):
            # fall back to the aligned CDS if lengths differ (indel in the CDS)
            aln = read_fasta(os.path.join(res, "alignments", "cds_alignment.fasta"))
            sa, sb = aln[a], aln[b]
        else:
            sa, sb = cds[a], cds[b]
        ndiff, nid = pair_stats(sa, sb)
        pdiff, pid = pair_stats(prot[a], prot[b]) if len(prot[a]) == len(prot[b]) else ([], float("nan"))
        syn = 0
        for pos, x, y in ndiff:
            c = (pos - 1) // 3
            ca, cb = sa[c * 3:c * 3 + 3], sb[c * 3:c * 3 + 3]
            kind, aa_a, aa_b = classify(ca, cb)
            syn += kind.startswith("synonymous")
            codon_rows.append((a, b, pos, f"{x}>{y}", c + 1, f"{ca}>{cb}", f"{aa_a}{c+1}{aa_b}", kind))
        rows.append((a, b, len(sa), len(ndiff), f"{nid:.3f}", syn, len(ndiff) - syn,
                     len(prot[a]), len(prot[b]), len(pdiff), f"{pid:.3f}"))

    with open(os.path.join(res, "summary.tsv"), "w") as fh:
        fh.write("seq_a\tseq_b\tcds_len\tcds_diffs\tcds_identity_pct\tsynonymous\tnonsynonymous\t"
                 "prot_len_a\tprot_len_b\tprot_diffs\tprot_identity_pct\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(res, "codon_differences.tsv"), "w") as fh:
        fh.write("seq_a\tseq_b\tcds_pos\tnt_change\tcodon_no\tcodon_change\taa_change\teffect\n")
        for r in codon_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    # --- report
    out = ["# PROCR / EPCR extraction and alignment", ""]
    out.append(f"Query protein: `{P['QUERY_ACC']}`  |  reference assembly: `{P['ASM_ACC']}` "
               f"(`{locus['contig']}:{locus['gene_start']}-{locus['gene_end']}`, "
               f"identity to query {locus['identity_to_query']})")
    out += ["", "## Sequences", "", "| label | source | protein length |", "|---|---|---|",
            f"| {ref_label} | {P['ASM_ACC']} assembly, miniprot | {len(prot[ref_label])} aa |"]
    for acc in runs:
        qc = dict(l.split("\t") for l in open(os.path.join(res, acc, "qc.tsv")).read().splitlines())
        out.append(f"| {acc} | WGS reads, reconstruction ({qc['gene_mean_depth']}x over gene) | {len(prot[acc])} aa |")
    out += ["", "## Pairwise comparisons", "",
            "| A | B | CDS bp | CDS diffs | CDS identity | syn | nonsyn | protein diffs | protein identity |",
            "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        out.append(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]}% | {r[5]} | {r[6]} | {r[9]} | {r[10]}% |")
    if codon_rows:
        out += ["", "## Nucleotide differences", "",
                "| A | B | CDS pos | change | codon | residue | effect |", "|---|---|---|---|---|---|---|"]
        for r in codon_rows:
            out.append(f"| {r[0]} | {r[1]} | c.{r[2]} | {r[3]} | {r[5]} | {r[6]} | {r[7]} |")
    if any(het.values()):
        out += ["", "## Heterozygous CDS sites in the reconstructions", "",
                "CDS positions where that individual's reads carry two alleles. The protein is",
                "translated ambiguity-aware: a site like this contributes its single residue when",
                "every resolution of the codon agrees, and an X only when they disagree.", "",
                "| run | CDS pos | IUPAC | codon | residue(s) | effect |", "|---|---|---|---|---|---|"]
        for acc in runs:
            for pos, b, codon, aas, eff in het[acc]:
                out.append(f"| {acc} | c.{pos} | {b} | {codon} | {aas} | {eff} |")
    with open(os.path.join(res, "heterozygous_cds_sites.tsv"), "w") as fh:
        fh.write("run\tcds_pos\tiupac\tcodon\tresidues\teffect\n")
        for acc in runs:
            for r in het[acc]:
                fh.write(acc + "\t" + "\t".join(str(x) for x in r) + "\n")

    out += ["", "## Read-reconstruction QC", ""]
    for acc in runs:
        qc = dict(l.split("\t") for l in open(os.path.join(res, acc, "qc.tsv")).read().splitlines())
        out.append(f"**{acc}** — pass-1 candidates {qc['candidate_reads_pass1']}, "
                   f"locus reads after competitive re-mapping {qc[[k for k in qc if k.startswith('locus_reads_pass2')][0]]}, "
                   f"mean gene depth {qc['gene_mean_depth']}x, locus variants {qc['locus_variants']} "
                   f"({qc['locus_snvs']} SNVs / {qc['locus_indels']} indels), "
                   f"gene-span {qc['gene_span_variants']}, CDS-exonic {qc['cds_exonic_variants']}.")
        out.append("")
        out.append("| exon | start | end | bp | mean depth | min depth | bp below min-depth |")
        out.append("|---|---|---|---|---|---|---|")
        for line in open(os.path.join(res, acc, "exon_depth.tsv")).read().splitlines()[1:]:
            out.append("| " + " | ".join(line.split("\t")) + " |")
        out.append("")
    with open(os.path.join(res, "report.md"), "w") as fh:
        fh.write("\n".join(out) + "\n")

    print("\n".join(out[:40]))
    print(f"\n[04] wrote {res}/summary.tsv, codon_differences.tsv, report.md and alignments/")


if __name__ == "__main__":
    sys.exit(main())
