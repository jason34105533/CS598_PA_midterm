"""Experiment C: placement (block / cyclic / packed) x transport (shm / rc).

usage: PYTHONPATH=analysis python3 analysis/plot_C.py   (writes figs/C_placement.png)
"""
import collections
import csv
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import SERIES

GRIDS = ["4x4x3", "1x1x48", "48x1x1"]
PLACEMENTS = [("block_spread", "block, 1 rank/CCD"), ("cyclic_spread", "cyclic, 1 rank/CCD"),
              ("block_packed", "block, 16 ranks on 2 CCDs")]
TLS = [("sm+self+rc", "shm"), ("rc", "rc (site default)")]


def main():
    runs = collections.defaultdict(list)
    for r in csv.DictReader(open("results/C.csv")):
        grid = "x".join(r[c] for c in ("Px", "Py", "Pz"))
        runs[(r["tls"], int(r["N"]), grid, r["tag"])].append(float(r["t_step"]) * 1e3)
    Ns = sorted({k[1] for k in runs})
    fig, axes = plt.subplots(len(Ns), 2, figsize=(10, 3.6 * len(Ns)), sharey="row")
    w = 0.26
    for i, N in enumerate(Ns):
        for j, (tls, tname) in enumerate(TLS):
            ax = axes[i][j]
            for k, (tag, label) in enumerate(PLACEMENTS):
                xs = [g + (k - 1) * (w + 0.02) for g in range(len(GRIDS))]
                ys = [statistics.median(runs[(tls, N, g, tag)]) for g in GRIDS]
                ax.bar(xs, ys, width=w, color=SERIES[k], label=label, edgecolor="#fcfcfb",
                       linewidth=1)
            ax.set_xticks(range(len(GRIDS)), GRIDS)
            ax.set_title(f"N = {N}, {tname}")
            ax.grid(True, axis="y")
            if j == 0:
                ax.set_ylabel("time per step (ms)")
    axes[0][0].legend(frameon=False, loc="upper left", fontsize=9)
    fig.suptitle("Experiment C: rank placement and MPI transport (P = 48)", fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/C_placement.png", dpi=150)


if __name__ == "__main__":
    main()
