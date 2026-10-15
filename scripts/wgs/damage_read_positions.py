#!/usr/bin/env python3
"""Read-level evidence at candidate SNV sites: for every read covering a site, which allele it
carries and where that base sits in the read and the fragment.

Sites: a TSV with chrom, pos (1-based), ref, alt, group (written by analysis/07_roq_damage.R).
Per read (primary, mapped, not duplicate / QC-fail): allele (REF / ALT / other), cycle (position
in the read from its 5' end, i.e. sequencing cycle, 1-based), dist_frag_end (bases to the nearer
fragment end, proper pairs only - from the template length), base quality, mapping quality,
mate (R1 / R2) and strand. Deamination during library prep (heat, end repair of single-stranded
overhangs) puts the ALT near fragment ends; in fresh-frozen tissue a negative result is
uninformative.

Usage: damage_read_positions.py --bam X.cram --fasta ref.fa --sites sites.tsv --out reads.tsv.gz
"""
import argparse
import csv
import gzip
import sys

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--sample", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    bam = pysam.AlignmentFile(a.bam, reference_filename=a.fasta)
    n_sites = n_rows = 0
    with open(a.sites) as f, gzip.open(a.out, "wt") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "chrom", "pos", "ref", "alt", "group", "allele", "cycle", "read_len",
                    "dist_frag_end", "bq", "mq", "mate", "strand"])
        for s in csv.DictReader(f, delimiter="\t"):
            n_sites += 1
            chrom, pos0 = s["chrom"], int(s["pos"]) - 1
            for r in bam.fetch(chrom, pos0, pos0 + 1):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail:
                    continue
                qpos = next((q for q, rp in r.get_aligned_pairs(matches_only=True) if rp == pos0), None)
                if qpos is None:                       # deletion / skip at the site
                    continue
                base = r.query_sequence[qpos].upper()
                allele = "REF" if base == s["ref"] else "ALT" if base == s["alt"] else "other"
                rlen = r.query_length
                cycle = rlen - qpos if r.is_reverse else qpos + 1
                dist = ""
                if r.is_proper_pair and r.template_length != 0:
                    fstart = min(r.reference_start, r.next_reference_start)
                    fend = fstart + abs(r.template_length) - 1
                    dist = min(pos0 - fstart, fend - pos0)
                w.writerow([a.sample, chrom, s["pos"], s["ref"], s["alt"], s["group"], allele, cycle, rlen,
                            dist, r.query_qualities[qpos], r.mapping_quality,
                            "R1" if r.is_read1 else "R2", "-" if r.is_reverse else "+"])
                n_rows += 1
    print(f"{a.sample}: {n_sites} sites, {n_rows} reads -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
