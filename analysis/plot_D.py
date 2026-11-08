"""Experiment D: exchange modes (blocking / nonblocking / overlap / overlap_test) x transport.

usage: PYTHONPATH=analysis python3 analysis/plot_D.py   (writes figs/D_overlap.png)
"""
import collections
import csv
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import SERIES

MODES = ["blocking", "nonblocking", "overlap", "overlap_test"]
COLORS = SERIES + ["#eda100"]  # categorical slots 1-4
GROUPS = [("4x4x3", "sm+self+rc", "4x4x3\nshm"), ("4x4x3", "rc", "4x4x3\nrc"),
          ("1x1x48", "sm+self+rc", "1x1x48\nshm"), ("1x1x48", "rc", "1x1x48\nrc")]
METRICS = [("t_step", "time per step (ms)"), ("wait_max", "Waitall / Sendrecv time (ms)")]


def main():
    runs = collections.defaultdict(list)
    for r in csv.DictReader(open("results/D.csv")):
        grid = "x".join(r[c] for c in ("Px", "Py", "Pz"))
        runs[(int(r["N"]), grid, r["tls"], r["mode"])].append(r)
    Ns = sorted({k[0] for k in runs})
    fig, axes = plt.subplots(len(Ns), 2, figsize=(11, 3.8 * len(Ns)))
    w = 0.2
    for i, N in enumerate(Ns):
        for j, (col, ylabel) in enumerate(METRICS):
            ax = axes[i][j]
            for k, mode in enumerate(MODES):
                xs = [g + (k - 1.5) * (w + 0.01) for g in range(len(GROUPS))]
                ys = [statistics.median(float(r[col]) for r in runs[(N, grid, tls, mode)]) * 1e3
                      for grid, tls, _ in GROUPS]
                ax.bar(xs, ys, width=w, color=COLORS[k], label=mode, edgecolor="#fcfcfb",
                       linewidth=1)
            ax.set_xticks(range(len(GROUPS)), [g[2] for g in GROUPS])
            ax.set_title(f"N = {N}: {ylabel.split(' (')[0]}")
            ax.set_ylabel(ylabel)
            ax.grid(True, axis="y")
    axes[0][0].legend(frameon=False, loc="upper left", fontsize=9)
    fig.suptitle("Experiment D: overlap of halo exchange and interior compute (P = 48)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/D_overlap.png", dpi=150)


if __name__ == "__main__":
    main()
