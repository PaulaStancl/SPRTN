#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# h2m_model.py - the only Python in cross_species/h2m: runs H2M (bioh2m 1.0.3, direction h2m)
# on a table prepared by 01_snv_coding_h2m.R. H2M is a Python package, so 01 calls this
# helper; everything before (variant selection, transcripts, orthologs) and after
# (interpretation, liftOver to GRCm38) is R.
#
#   python h2m_model.py --input in.tsv --out raw.tsv --failed failed.tsv \
#       --human-fasta GRCh38.fa --mouse-fasta GRCm39.fa --human-gtf h.gtf.gz --mouse-gtf m.gtf.gz \
#       --workdir tmp/ [--max-alternative 5] [--no-alternative]
#
# in.tsv: index tx_id_h tx_id_m start_h end_h ref_seq_h alt_seq_h type_h (MAF-style, H2M's input).
# H2M is given (a) gffutils DBs holding only the requested transcripts, built from the
# Ensembl GTFs, and (b) a lazy FASTA reader with the same shape as h2m.genome_loader()
# (records[index_list[i]].seq, human chr1-22,X,Y / mouse 1-19,X,Y), so neither GENCODE DBs
# nor whole genomes in memory are needed. Human genome version passed to H2M: 38.
# ---------------------------------------------------------------------------
import argparse
import gzip
import re
from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace

import pandas as pd


class LazyGenome:
    """records[index_list[i]].seq like h2m.genome_loader, one chromosome in memory at a time."""

    def __init__(self, fasta, n_autosomes, keep=2):
        import pysam
        self.fa = pysam.FastaFile(str(fasta))
        refs = set(self.fa.references)
        self.names = []
        for n in [str(i) for i in range(1, n_autosomes + 1)] + ["X", "Y"]:
            hit = next((x for x in (n, f"chr{n}") if x in refs), None)
            if hit is None:
                raise SystemExit(f"{fasta}: chromosome {n} not found")
            self.names.append(hit)
        self.cache, self.keep = OrderedDict(), keep
        self.index_list = list(range(len(self.names)))

    def __len__(self):
        return len(self.names)

    def __getitem__(self, i):
        from Bio.Seq import Seq
        if i not in self.cache:
            self.cache[i] = SimpleNamespace(seq=Seq(self.fa.fetch(self.names[i])), description=self.names[i])
            while len(self.cache) > self.keep:
                self.cache.popitem(last=False)
        return self.cache[i]


def norm_chrom(c):
    c = c[3:] if c.startswith("chr") else c
    return "chrM" if c in ("MT", "M") else f"chr{c}"


def mini_db(gtf, tx_ids, dest):
    """gffutils DB with only these transcripts (transcript / exon / CDS / start / stop lines)."""
    import gffutils
    keep_feat = {"transcript", "exon", "CDS", "start_codon", "stop_codon"}
    tx_ids = set(tx_ids)
    sub = Path(str(dest) + ".gtf")
    n = 0
    with gzip.open(gtf, "rt") as fi, open(sub, "w") as fo:
        for line in fi:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] not in keep_feat:
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if m and m.group(1) in tx_ids:
                f[0] = norm_chrom(f[0])
                fo.write("\t".join(f))
                n += 1
    print(f"{sub.name}: {n} GTF lines for {len(tx_ids)} transcripts")
    return gffutils.create_db(str(sub), dbfn=str(dest), force=True, keep_order=True,
                              merge_strategy="create_unique", disable_infer_genes=True,
                              disable_infer_transcripts=True,
                              id_spec={"transcript": "transcript_id", "gene": "gene_id"})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for a in ("--input", "--out", "--failed", "--human-fasta", "--mouse-fasta", "--human-gtf", "--mouse-gtf", "--workdir"):
        ap.add_argument(a, required=True)
    ap.add_argument("--max-alternative", type=int, default=5)
    ap.add_argument("--no-alternative", action="store_true")
    a = ap.parse_args()
    import bioh2m as h2m

    inp = pd.read_csv(a.input, sep="\t", dtype={"ref_seq_h": str, "alt_seq_h": str, "tx_id_h": str, "tx_id_m": str})
    wd = Path(a.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    db_h = mini_db(a.human_gtf, set(inp.tx_id_h), wd / "h2m_human.db")
    db_m = mini_db(a.mouse_gtf, set(inp.tx_id_m), wd / "h2m_mouse.db")
    gh, gm = LazyGenome(a.human_fasta, 22), LazyGenome(a.mouse_fasta, 19)
    res, fail = h2m.model_batch(inp, gh, gh.index_list, gm, gm.index_list, db_h, db_m, 38,
                                search_alternative=not a.no_alternative, max_alternative=a.max_alternative)
    res.to_csv(a.out, sep="\t", index=False)
    fail.to_csv(a.failed, sep="\t", index=False)
    print(f"H2M: {len(inp)} inputs, {len(res)} result rows, {len(fail)} failed")


if __name__ == "__main__":
    main()
