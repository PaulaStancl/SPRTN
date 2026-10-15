#!/usr/bin/env python3
"""GIV score (Chen et al. 2017, Science 355:752) from one BAM/CRAM: read-1 vs read-2 mismatch
rates per substitution type, in read orientation.

Single-strand damage before library prep is read on one strand of the original fragment, so it
shows up in read 1 and its complement in read 2: 8-oxoG -> G>T excess in R1 (C>A in R2),
deamination -> C>T excess in R1 (G>A in R2). Per type X>Y:
    rate_Rn = (X>Y mismatches in Rn) / (X reference bases in Rn), bases in read orientation
    GIV     = rate_R1 / rate_R2          (Chen et al.: GIV > 1.5 = damaged; 1 = undamaged)
Germline and somatic variants count equally in R1 and R2, so they pull GIV towards 1, never
away from it. Reads: primary, mapped, properly paired, not duplicate / QC-fail, MAPQ >= --min-mq;
bases: BQ >= --min-bq, aligned (M/=/X), reference not N.

Usage: damage_giv.py --bam X.cram --fasta ref.fa --regions chr20 [--max-reads 5000000] --out X.giv.tsv
"""
import argparse
import sys

import numpy as np
import pysam

LUT = np.full(256, 4, dtype=np.uint8)          # A C G T -> 0..3, anything else (N) -> 4
for i, b in enumerate(b"ACGT"):
    LUT[b] = i
    LUT[ord(chr(b).lower())] = i
BASES = "ACGT"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True, help="reference the BAM/CRAM was aligned to")
    ap.add_argument("--regions", default="chr20", help="comma-separated contigs or chr:start-end")
    ap.add_argument("--max-reads", type=int, default=5_000_000, help="reads counted (after filters), all regions together")
    ap.add_argument("--min-bq", type=int, default=30)
    ap.add_argument("--min-mq", type=int, default=20)
    ap.add_argument("--sample", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    bam = pysam.AlignmentFile(a.bam, reference_filename=a.fasta)
    fa = pysam.FastaFile(a.fasta)
    counts = np.zeros(2 * 5 * 5, dtype=np.int64)  # mate (0 = R1, 1 = R2) x ref x read base
    n_reads = 0
    for region in a.regions.split(","):
        contig = region.split(":")[0]
        ref = LUT[np.frombuffer(fa.fetch(contig).encode(), dtype=np.uint8)]
        for r in bam.fetch(region=region):
            if (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail
                    or not r.is_proper_pair or r.mapping_quality < a.min_mq or r.query_qualities is None):
                continue
            seq = LUT[np.frombuffer(r.query_sequence.encode(), dtype=np.uint8)]
            bq = np.asarray(r.query_qualities)
            q_idx, r_idx = [], []
            qp, rp = 0, r.reference_start
            for op, ln in r.cigartuples:
                if op in (0, 7, 8):                   # M = X: aligned
                    q_idx.append(np.arange(qp, qp + ln)); r_idx.append(np.arange(rp, rp + ln))
                    qp += ln; rp += ln
                elif op in (1, 4):                    # I S: query only
                    qp += ln
                elif op in (2, 3):                    # D N: reference only
                    rp += ln
            if not q_idx:
                continue
            qi = np.concatenate(q_idx); ri = np.concatenate(r_idx)
            rb, qb = ref[ri], seq[qi]
            keep = (bq[qi] >= a.min_bq) & (rb < 4) & (qb < 4)
            rb, qb = rb[keep], qb[keep]
            if r.is_reverse:                          # read orientation: complement (A<->T, C<->G = 3 - x)
                rb, qb = 3 - rb, 3 - qb
            mate = 0 if r.is_read1 else 1
            counts += np.bincount(mate * 25 + rb.astype(np.int64) * 5 + qb, minlength=50)
            n_reads += 1
            if n_reads >= a.max_reads:
                break
        if n_reads >= a.max_reads:
            break

    c = counts.reshape(2, 5, 5)
    with open(a.out, "w") as out:
        out.write("sample\ttype\tref_bases_R1\tmismatch_R1\tref_bases_R2\tmismatch_R2\trate_R1\trate_R2\tGIV\treads\n")
        for x in range(4):
            for y in range(4):
                if x == y:
                    continue
                tot1, tot2 = c[0, x, :4].sum(), c[1, x, :4].sum()
                m1, m2 = c[0, x, y], c[1, x, y]
                r1 = m1 / tot1 if tot1 else float("nan")
                r2 = m2 / tot2 if tot2 else float("nan")
                giv = r1 / r2 if r2 else float("nan")
                out.write(f"{a.sample}\t{BASES[x]}_{BASES[y]}\t{tot1}\t{m1}\t{tot2}\t{m2}\t{r1:.4g}\t{r2:.4g}\t{giv:.4g}\t{n_reads}\n")
    print(f"{a.sample}: {n_reads} reads -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
