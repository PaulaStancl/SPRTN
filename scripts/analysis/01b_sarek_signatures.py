#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 01b_sarek_signatures.py - COSMIC mutational signatures in the sarek PASS sets.
# Run after 01_sarek.R (section 7 writes the inputs), in the sigprofiler env:
#
#   micromamba activate /common/WORK/pstancl/envs/sigprofiler
#   python 01b_sarek_signatures.py                     # SBS96 + DBS78 + ID83
#   python 01b_sarek_signatures.py --contexts SBS96    # one context only
#   python 01b_sarek_signatures.py --pipeline oncoanalyser   # SAGE PASS set from 02_oncoanalyser.R
#
# Input  <OUT>/<pipeline>/signatures/input/vcf/<set>.vcf - mutect2 (all PASS), strelka (all PASS),
#        mutect2_strelka (PASS in both) - or, for oncoanalyser, sage; each file is one
#        "sample", fitted on its own.
# Output <OUT>/<pipeline>/signatures/sigprofiler/
#   matrix_generator/output/{SBS,DBS,ID}/   SigProfilerMatrixGenerator matrices
#   SBS96/<set>/ DBS78/<set>/ ID83/<set>/   SigProfilerAssignment cosmic_fit results - one
#                                           cosmic_fit call per set, so each is fitted on its own
#   activities_<context>.csv                signature, mutations assigned, fraction - all sets
#   fit_stats_<context>.csv                 cosine similarity etc. of each set's reconstruction
#   sbs96_compare.csv / sbs96_crosscheck.csv  SigProfiler vs 01_sarek.R SBS96, per class / per set
#                                           (comparison only - fits use SigProfiler's matrices)
#
# cosmic_fit is SigProfilerAssignment's single-sample refitting (the SigProfilerSingleSample
# package is deprecated and was folded into it): each set's profile is explained as a
# non-negative mix of COSMIC reference signatures - no de novo extraction, which needs a cohort.
#
# All matrices (SBS96, DBS78, ID83) are built by SigProfilerMatrixGenerator from the VCFs,
# on its own GRCh38 - installed once (a few GB):
#   python -c "from SigProfilerMatrixGenerator import install as g; g.install('GRCh38')"
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

CONTEXTS = {"SBS96": ("SBS", "SBS96"), "DBS78": ("DBS", "DBS78"), "ID83": ("ID", "ID83")}

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--pipeline", default="sarek", choices=["sarek", "oncoanalyser"],
                help="whose PASS sets: sarek (01_sarek.R) or oncoanalyser (02_oncoanalyser.R)")
ap.add_argument("--contexts", default="SBS96,DBS78,ID83", help="comma-separated: SBS96, DBS78, ID83")
ap.add_argument("--cosmic-version", type=float, default=None, help="COSMIC version (default: the package's)")
ap.add_argument("--exclude-subgroups", default="",
                help="comma-separated signature subgroups to leave out, e.g. Chemotherapy_signatures")
ap.add_argument("--cpu", type=int, default=int(os.environ.get("NCPUS", "-1")))
args = ap.parse_args()
SIG = os.path.join(OUT, args.pipeline, "signatures")
IN_DIR = os.path.join(SIG, "input")
FIT_DIR = os.path.join(SIG, "sigprofiler")

from SigProfilerAssignment import Analyzer as Analyze  # noqa: E402 - after argparse so --help is fast

vcf_dir = os.path.join(IN_DIR, "vcf")
if not glob.glob(os.path.join(vcf_dir, "*.vcf")):
    sys.exit(f"No PASS-set VCFs in {vcf_dir} - run 01_sarek.R (section 7) first")
os.makedirs(FIT_DIR, exist_ok=True)

# ---- matrices: SigProfilerMatrixGenerator --------------------------------------
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen  # noqa: E402
# the generator writes output/ next to its input, so work on a copy of the VCFs
mg_dir = os.path.join(FIT_DIR, "matrix_generator")
shutil.rmtree(mg_dir, ignore_errors=True)
shutil.copytree(vcf_dir, mg_dir)
try:
    matGen.SigProfilerMatrixGeneratorFunc(PATIENT, "GRCh38", mg_dir + "/", exome=False, bed_file=None,
                                          chrom_based=False, plot=True, tsb_stat=False, seqInfo=False)
except Exception as e:  # most often: the GRCh38 reference is not installed
    sys.exit(f"SigProfilerMatrixGenerator failed: {e}\nInstall its genome once with\n"
             "  python -c \"from SigProfilerMatrixGenerator import install as g; g.install('GRCh38')\"")
matrices = {}
for ctx in args.contexts.split(","):
    sub, name = CONTEXTS[ctx.strip()]
    matrices[ctx.strip()] = os.path.join(mg_dir, "output", sub, f"{PATIENT}.{name}.all")

# Comparison only - the fits below use the matrix generator's matrices, never 01_sarek.R's:
# SBS96 per class and set, SigProfiler vs 01_sarek.R (sbs96_compare.csv, every row) and
# totals per set (sbs96_crosscheck.csv).
r_matrix = os.path.join(IN_DIR, f"{PATIENT}.SBS96.from_R.txt")
mg96 = matrices.get("SBS96")
if mg96 and os.path.exists(mg96) and os.path.exists(r_matrix):
    a = pd.read_csv(mg96, sep="\t", index_col=0)
    b = pd.read_csv(r_matrix, sep="\t", index_col=0).reindex(index=a.index, columns=a.columns).fillna(0).astype(int)
    cmp = (a.rename_axis("MutationType").reset_index().melt(id_vars="MutationType", var_name="set", value_name="sigprofiler")
           .merge(b.rename_axis("MutationType").reset_index().melt(id_vars="MutationType", var_name="set", value_name="from_R"),
                  on=["MutationType", "set"]))
    cmp["difference"] = cmp["sigprofiler"] - cmp["from_R"]
    cmp.to_csv(os.path.join(FIT_DIR, "sbs96_compare.csv"), index=False)
    cc = cmp.groupby("set").agg(sigprofiler=("sigprofiler", "sum"), from_R=("from_R", "sum"),
                                rows_differing=("difference", lambda d: int((d != 0).sum()))).reset_index()
    cc.to_csv(os.path.join(FIT_DIR, "sbs96_crosscheck.csv"), index=False)
    print(cc.to_string(index=False))
    if cc["rows_differing"].any():
        print("SBS96 differs between SigProfiler and 01_sarek.R - rows with difference != 0 in sbs96_compare.csv")

# ---- COSMIC refit per context -----------------------------------------------------
for ctx, path in matrices.items():
    if not os.path.exists(path):
        print(f"{ctx}: no matrix ({path}) - skipped")
        continue
    m = pd.read_csv(path, sep="\t", index_col=0)
    print(f"{ctx}: mutations per set\n{m.sum().to_string()}")
    acts, stats = [], []
    shutil.rmtree(os.path.join(FIT_DIR, ctx), ignore_errors=True)   # no fits left from sets of an earlier run
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
