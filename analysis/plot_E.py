"""Experiment E: manual packing vs MPI derived datatypes, by face direction.

Stacked bars: compute, Waitall, pack + unpack, each averaged over ranks (the stack
is the average time per step); the dash marks t_step, the slowest rank. With datatypes the packing happens
inside MPI and shows up in the Waitall time.

usage: PYTHONPATH=analysis python3 analysis/plot_E.py   (writes figs/E_packing.png)
"""
import collections
import csv
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import INK2, SERIES

GRIDS = [("48x1x1", "48x1x1\nx-faces"), ("1x48x1", "1x48x1\ny-faces"),
         ("1x1x48", "1x1x48\nz-faces"), ("4x4x3", "4x4x3\nall")]
PARTS = [("comp_int_max", "compute"), ("wait_max", "Waitall"), ("pack", "pack + unpack")]


def main():
    runs = collections.defaultdict(list)
    for r in csv.DictReader(open("results/E.csv")):
        grid = "x".join(r[c] for c in ("Px", "Py", "Pz"))
        runs[(int(r["N"]), grid, r["pack"])].append(r)
    Ns = sorted({k[0] for k in runs})
    fig, axes = plt.subplots(1, len(Ns), figsize=(5.6 * len(Ns), 4.4))
    w = 0.36
    for ax, N in zip(axes, Ns):
        for gi, (grid, _) in enumerate(GRIDS):
            for pi, pack in enumerate(["manual", "datatype"]):
                rs = runs[(N, grid, pack)]
                med = lambda c: statistics.median(float(r[c]) for r in rs) * 1e3
                vals = [med("comp_int_avg"), med("wait_avg"), med("pack_avg") + med("unpack_avg")]
                x = gi + (pi - 0.5) * (w + 0.04)
                bottom = 0
                for k, v in enumerate(vals):
                    ax.bar(x, v, width=w, bottom=bottom, color=SERIES[k], edgecolor="#fcfcfb",
                           linewidth=1, label=PARTS[k][1] if (gi == 0 and pi == 0) else None)
                    bottom += v
                t = med("t_step")
                ax.plot([x - w / 2, x + w / 2], [t, t], color=INK2, lw=1.5)
                ax.text(x, max(t, bottom), "M" if pack == "manual" else "D", ha="center",
                        va="bottom", fontsize=9, color=INK2)
        ax.set_xticks(range(len(GRIDS)), [g[1] for g in GRIDS])
        ax.set(title=f"N = {N}")
        ax.grid(True, axis="y")
    axes[0].set_ylabel("time per step (ms), average over ranks")
    axes[0].legend(frameon=False, loc="upper right", fontsize=9)
    fig.suptitle("Experiment E: packing of non-contiguous faces (P = 48).  M = manual, D = datatype, dash = slowest rank", fontsize=11)
    fig.tight_layout()
    fig.savefig("figs/E_packing.png", dpi=150)


if __name__ == "__main__":
    main()
