#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 01b_sarek_signatures.py - COSMIC mutational signatures in the sarek PASS sets.
# Run after 01_sarek.R (section 6 writes the inputs), in the sigprofiler env:
#
#   micromamba activate /common/WORK/pstancl/envs/sigprofiler
#   python 01b_sarek_signatures.py                  # SBS96 + DBS78 + ID83
#   python 01b_sarek_signatures.py --from-r-matrix  # SBS96 only, from 01_sarek.R's matrix
#
# Input  <OUT>/sarek/signatures/input/vcf/<set>.vcf - mutect2 (all PASS), strelka (all PASS),
#        mutect2_strelka (PASS in both); each file is one "sample".
# Output <OUT>/sarek/signatures/sigprofiler/
#   matrix_generator/output/{SBS,DBS,ID}/   SigProfilerMatrixGenerator matrices
#   SBS96/<set>/ DBS78/<set>/ ID83/<set>/   SigProfilerAssignment cosmic_fit results - one
#                                           cosmic_fit call per set, so each is fitted on its own
#   activities_<context>.csv                signature, mutations assigned, fraction - all sets
#   fit_stats_<context>.csv                 cosine similarity etc. of each set's reconstruction
#   sbs96_crosscheck.csv                    matrix generator vs 01_sarek.R SBS96 counts
#
# cosmic_fit is SigProfilerAssignment's single-sample refitting (the SigProfilerSingleSample
# package is deprecated and was folded into it): each set's profile is explained as a
# non-negative mix of COSMIC reference signatures - no de novo extraction, which needs a cohort.
#
# The matrix generator needs its own copy of GRCh38 (a one-off download of a few GB):
#   python -c "from SigProfilerMatrixGenerator import install as g; g.install('GRCh38')"
# Without it, --from-r-matrix fits SBS96 from the matrix 01_sarek.R built from the GATK FASTA.
# ---------------------------------------------------------------------------
import argparse
import glob
import os
import re
import shutil
import sys

import pandas as pd

PATIENT = "RJALS"
RESULTS = os.environ.get("SPRTN_RESULTS", "/common/WORK/pstancl/projects/MariaBoskovic/SPRTN/results/wgs")
OUT = os.environ.get("SPRTN_OUT", os.path.join(RESULTS, "analysis"))
if re.search(r"OneDrive|CloudStorage|Dropbox|iCloud|Google Drive", os.path.abspath(OUT), re.I):
    sys.exit(f"Refusing to use a cloud-synced folder for patient data: {OUT}")

SIG = os.path.join(OUT, "sarek", "signatures")
IN_DIR = os.path.join(SIG, "input")
FIT_DIR = os.path.join(SIG, "sigprofiler")
CONTEXTS = {"SBS96": ("SBS", "SBS96"), "DBS78": ("DBS", "DBS78"), "ID83": ("ID", "ID83")}

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--contexts", default="SBS96,DBS78,ID83", help="comma-separated: SBS96, DBS78, ID83")
ap.add_argument("--from-r-matrix", action="store_true",
                help="fit SBS96 from 01_sarek.R's matrix; no SigProfiler genome needed")
ap.add_argument("--cosmic-version", type=float, default=None, help="COSMIC version (default: the package's)")
ap.add_argument("--exclude-subgroups", default="",
                help="comma-separated signature subgroups to leave out, e.g. Chemotherapy_signatures")
ap.add_argument("--cpu", type=int, default=int(os.environ.get("NCPUS", "-1")))
args = ap.parse_args()

from SigProfilerAssignment import Analyzer as Analyze  # noqa: E402 - after argparse so --help is fast

vcf_dir = os.path.join(IN_DIR, "vcf")
if not glob.glob(os.path.join(vcf_dir, "*.vcf")):
    sys.exit(f"No PASS-set VCFs in {vcf_dir} - run 01_sarek.R (section 6) first")
os.makedirs(FIT_DIR, exist_ok=True)
r_matrix = os.path.join(IN_DIR, f"{PATIENT}.SBS96.from_R.txt")

# ---- matrices -----------------------------------------------------------------
matrices = {}
if args.from_r_matrix:
    matrices["SBS96"] = r_matrix
else:
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
    # the generator writes output/ next to its input, so work on a copy of the VCFs
    mg_dir = os.path.join(FIT_DIR, "matrix_generator")
    shutil.rmtree(mg_dir, ignore_errors=True)
    shutil.copytree(vcf_dir, mg_dir)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PATIENT, "GRCh38", mg_dir + "/", exome=False, bed_file=None,
                                              chrom_based=False, plot=True, tsb_stat=False, seqInfo=False)
    except Exception as e:  # most often: the GRCh38 reference is not installed
        sys.exit(f"SigProfilerMatrixGenerator failed: {e}\nInstall its genome once with\n"
                 "  python -c \"from SigProfilerMatrixGenerator import install as g; g.install('GRCh38')\"\n"
                 "or run with --from-r-matrix (SBS96 only).")
    for ctx in args.contexts.split(","):
        sub, name = CONTEXTS[ctx.strip()]
        matrices[ctx.strip()] = os.path.join(mg_dir, "output", sub, f"{PATIENT}.{name}.all")

    # cross-check: the generator's SBS96 vs 01_sarek.R's (same mutations, other reference copy)
    mg96 = matrices.get("SBS96")
    if mg96 and os.path.exists(mg96) and os.path.exists(r_matrix):
        a = pd.read_csv(mg96, sep="\t", index_col=0)
        b = pd.read_csv(r_matrix, sep="\t", index_col=0).reindex(index=a.index, columns=a.columns).fillna(0)
        cc = pd.DataFrame({"set": a.columns, "matrix_generator": a.sum().values, "from_R": b.sum().values,
                           "rows_differing": (a != b).sum().values})
        cc.to_csv(os.path.join(FIT_DIR, "sbs96_crosscheck.csv"), index=False)
        print(cc.to_string(index=False))
        if cc["rows_differing"].any():
            print("WARNING: SBS96 counts differ between SigProfiler and 01_sarek.R - see sbs96_crosscheck.csv")

# ---- COSMIC refit per context -----------------------------------------------------
for ctx, path in matrices.items():
    if not os.path.exists(path):
        print(f"{ctx}: no matrix ({path}) - skipped")
        continue
    m = pd.read_csv(path, sep="\t", index_col=0)
    print(f"{ctx}: mutations per set\n{m.sum().to_string()}")
    acts, stats = [], []
    # one cosmic_fit call per set, on a one-column matrix: each set is fitted on its own and
    # gets its own folder, <context>/<set>/
    for st in m.columns:
        if m[st].sum() == 0:                   # e.g. no doublets in a set - cosmic_fit would stop
            print(f"{ctx} {st}: no mutations - skipped")
            continue
        out = os.path.join(FIT_DIR, ctx, st)
        os.makedirs(out, exist_ok=True)
        m_path = os.path.join(out, f"{PATIENT}.{ctx}.{st}.input.txt")
        m[[st]].to_csv(m_path, sep="\t")
        kw = dict(samples=m_path, output=out, input_type="matrix", genome_build="GRCh38", exome=False,
                  make_plots=True, cpu=args.cpu)
        if args.cosmic_version:
            kw["cosmic_version"] = args.cosmic_version
        if args.exclude_subgroups:
            kw["exclude_signature_subgroups"] = [s.strip() for s in args.exclude_subgroups.split(",")]
        Analyze.cosmic_fit(**kw)

        # activities: mutations assigned to each signature
        act = [f for f in glob.glob(os.path.join(out, "**", "*Activities.txt"), recursive=True)
               if "Assignment_Solution" in f]
        if act:
            a = pd.read_csv(act[0], sep="\t", index_col=0)
            long = a.reset_index().melt(id_vars=a.index.name or "index", var_name="signature", value_name="mutations")
            long.columns = ["set", "signature", "mutations"]
            acts.append(long[long["mutations"] > 0])
        sf = glob.glob(os.path.join(out, "**", "*Samples_Stats.txt"), recursive=True)
        if sf:
            stats.append(pd.read_csv(sf[0], sep="\t"))

    # all sets of this context in one table each
    if acts:
        long = pd.concat(acts)
        long["fraction"] = long["mutations"] / long.groupby("set")["mutations"].transform("sum")
        long = long.sort_values(["set", "mutations"], ascending=[True, False])
        long.to_csv(os.path.join(FIT_DIR, f"activities_{ctx}.csv"), index=False)
        print(long.to_string(index=False))
    if stats:
        pd.concat(stats).to_csv(os.path.join(FIT_DIR, f"fit_stats_{ctx}.csv"), index=False)

print(f"done: {FIT_DIR}")
