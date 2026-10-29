#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# h2m_visualize.py - H2M's own result plot (h2m.model(show_sequence=True) + h2m.visualization, as in
# the H2M tutorial) for selected variants; called by check_class4.R --h2m-plots.
#
#   python h2m_visualize.py --input in.tsv --outdir plots/ --status status.tsv \
#       --human-fasta ... --mouse-fasta ... --human-gtf ... --mouse-gtf ... --workdir tmp/ [--flank 4 --print-size 2]
#
# in.tsv: label tx_id_h tx_id_m start_h end_h ref_seq_h alt_seq_h type_h (MAF-style, as for h2m_model.py).
# Writes <outdir>/<label>.png per variant and status.tsv (label, h2m_class, plotted, message).
# visualization() is tried for every variant, also classes H2M did not model (3-6); status.tsv records
# the outcome (a plot, no figure, or the error H2M raised).
# ---------------------------------------------------------------------------
import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h2m_model import LazyGenome, mini_db  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for a in ("--input", "--outdir", "--status", "--human-fasta", "--mouse-fasta", "--human-gtf", "--mouse-gtf", "--workdir"):
        ap.add_argument(a, required=True)
    ap.add_argument("--flank", type=int, default=4)
    ap.add_argument("--print-size", type=int, default=2)
    a = ap.parse_args()
    import bioh2m as h2m

    inp = pd.read_csv(a.input, sep="\t", dtype=str)
    out, wd = Path(a.outdir), Path(a.workdir)
    out.mkdir(parents=True, exist_ok=True)
    wd.mkdir(parents=True, exist_ok=True)
    db_h = mini_db(a.human_gtf, set(inp.tx_id_h), wd / "viz_human.db")
    db_m = mini_db(a.mouse_gtf, set(inp.tx_id_m), wd / "viz_mouse.db")
    gh, gm = LazyGenome(a.human_fasta, 22), LazyGenome(a.mouse_fasta, 19)

    rows = []
    for r in inp.itertuples():
        png = out / f"{r.label}.png"
        cls, plotted, msg = None, False, ""
        try:
            res = h2m.model(gh, gh.index_list, gm, gm.index_list, db_h, db_m, r.tx_id_h, r.tx_id_m,
                            int(r.start_h), int(r.end_h), r.ref_seq_h, r.alt_seq_h, ty_h=r.type_h, ver=38,
                            show_sequence=True)
            cls = res[0].get("class")
            # always try, also for classes H2M did not model, and record what happens
            plt.close("all")
            plt.show = lambda *x, **k: plt.savefig(png, dpi=130, bbox_inches="tight")
            h2m.visualization(res, flank_size=a.flank, print_size=a.print_size)
            plotted = png.exists()
            msg = "ok" if plotted else f"visualization returned no figure (class {cls})"
        except Exception as e:  # one failing variant must not stop the others
            msg = f"error: {type(e).__name__}: {e}"
        plt.close("all")
        rows.append({"label": r.label, "h2m_class": cls, "plotted": plotted, "message": msg})
        print(f"{r.label}: {msg}")
    pd.DataFrame(rows).to_csv(a.status, sep="\t", index=False)


if __name__ == "__main__":
    main()
