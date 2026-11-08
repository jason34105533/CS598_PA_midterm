"""Experiment A: a two-term model of the time per step,
    t_step = a + b * (halo bytes per rank) + c * log2(Px),
fitted per N over the 45 grids. b is an effective inverse bandwidth under the real
concurrency (all 48 ranks exchanging at once); c is the penalty per halving of
the unit-stride (x) extent. Prints the fit and writes figs/A_model.png.

usage: PYTHONPATH=analysis python3 analysis/model_A.py
"""
import collections
import csv
import math
import statistics

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import INK2, SERIES

MARKERS = ["o", "s", "^"]


def main():
    runs = collections.defaultdict(list)
    for r in csv.DictReader(open("results/A.csv")):
        runs[(int(r["N"]), int(r["Px"]), int(r["Py"]), int(r["Pz"]))].append(r)
    Ns = sorted({k[0] for k in runs})
    fig, axes = plt.subplots(1, len(Ns), figsize=(4.2 * len(Ns), 4.2))
    for ax, (i, N) in zip(axes, enumerate(Ns)):
        keys = [k for k in runs if k[0] == N]
        y = np.array([statistics.median(float(r["t_step"]) for r in runs[k]) * 1e3 for k in keys])
        halo = np.array([float(runs[k][0]["halo_bytes_max"]) / 1e6 for k in keys])
        lpx = np.array([math.log2(k[1]) for k in keys])
        off = np.array([float(runs[k][0]["bytes_offnode"]) / 3 / 1e6 for k in keys])
        A = np.c_[np.ones(len(y)), halo, lpx]
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        pred = A @ coef
        r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        r2_halo = 1 - ((y - np.polyval(np.polyfit(halo, y, 1), halo)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        print(f"N = {N}: t = {coef[0]:.3f} + {coef[1]:.3f} ms/MB * halo + {coef[2]:.3f} ms * log2(Px)"
              f"   R2 = {r2:.2f} (halo only {r2_halo:.2f});"
              f" effective bandwidth {1 / coef[1]:.2f} GB/s per rank;"
              f" off-node MB per node {off.min():.2f}-{off.max():.2f}")
        lim = (0, max(y.max(), pred.max()) * 1.08)
        ax.plot(lim, lim, color=INK2, lw=1, ls="--")
        ax.scatter(pred, y, s=34, color=SERIES[i], marker=MARKERS[i], edgecolor="#fcfcfb",
                   linewidth=0.8, zorder=3)
        ax.set(xlim=lim, ylim=lim, xlabel="model (ms)", title=f"N = {N}   R² = {r2:.2f}")
        ax.grid(True)
    axes[0].set_ylabel("measured time per step (ms)")
    fig.suptitle("Experiment A: t = a + b·halo + c·log2(Px), fitted per N", fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/A_model.png", dpi=150)


if __name__ == "__main__":
    main()
