"""Experiment D2: overlap gain vs N, and the peel mode.

Gain = t_step(nonblocking) / t_step(best of overlap_test and peel), per grid and
transport, at N = 240, 480, 960. Reads results/D.csv and results/D2.csv.

usage: PYTHONPATH=analysis python3 analysis/plot_D2.py   (writes figs/D2_gain.png)
"""
import collections
import csv
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import INK2, MARKERS, SERIES

COLORS = SERIES + ["#eda100"]
LINES = [("4x4x3", "rc", "4x4x3, rc"), ("1x1x48", "rc", "1x1x48, rc"),
         ("4x4x3", "sm+self+rc", "4x4x3, shm"), ("1x1x48", "sm+self+rc", "1x1x48, shm")]
NS = [240, 480, 960]


def main():
    runs = collections.defaultdict(list)
    for path in ("results/D.csv", "results/D2.csv"):
        for r in csv.DictReader(open(path)):
            if r["pack"] != "manual":
                continue
            grid = "x".join(r[c] for c in ("Px", "Py", "Pz"))
            runs[(int(r["N"]), grid, r["tls"], r["mode"])].append(float(r["t_step"]))
    t = {k: statistics.median(v) for k, v in runs.items()}

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for i, (grid, tls, label) in enumerate(LINES):
        g_test = [t[(N, grid, tls, "nonblocking")] / t[(N, grid, tls, "overlap_test")] for N in NS]
        g_peel = [t[(N, grid, tls, "nonblocking")] / t[(N, grid, tls, "peel")] for N in NS]
        axes[0].plot(NS, g_test, marker=(MARKERS + ["D"])[i], color=COLORS[i], label=label)
        axes[1].plot(NS, g_peel, marker=(MARKERS + ["D"])[i], color=COLORS[i], label=label)
    for ax, title in zip(axes, ["overlap_test (send u, then interior)",
                                "peel (boundary first, then interior)"]):
        ax.axhline(1, color=INK2, lw=1, ls="--")
        ax.set(xscale="log", xlabel="N", title=title, ylim=(0.5, 2.8))
        ax.set_xticks(NS, [str(n) for n in NS])
        ax.set_xticks([], minor=True)
        ax.grid(True)
    axes[0].set_ylabel("speedup over nonblocking")
    axes[0].legend(frameon=False, loc="upper right", fontsize=9)
    fig.suptitle("Overlap gain vs problem size (P = 48)", fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/D2_gain.png", dpi=150)


if __name__ == "__main__":
    main()
