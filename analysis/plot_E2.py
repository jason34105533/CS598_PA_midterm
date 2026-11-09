"""Experiment E2: packing through index lists (indirect addressing) vs structured
loops and subarray datatypes. manual / indexed / indexed_dt come from
results/E2.csv; datatype (subarray) from results/E.csv.

usage: PYTHONPATH=analysis python3 analysis/plot_E2.py   (writes figs/E2_indexed.png)
"""
import collections
import csv
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import SERIES

COLORS = SERIES + ["#eda100"]
GRIDS = [("48x1x1", "48x1x1\nx-faces"), ("1x48x1", "1x48x1\ny-faces"),
         ("1x1x48", "1x1x48\nz-faces"), ("4x4x3", "4x4x3\nall")]
PACKS = [("E2", "manual", "manual loops"), ("E2", "indexed", "index list, by hand"),
         ("E", "datatype", "subarray datatype"), ("E2", "indexed_dt", "index list as datatype")]


def main():
    runs = collections.defaultdict(list)
    for src in ("E", "E2"):
        for r in csv.DictReader(open(f"results/{src}.csv")):
            grid = "x".join(r[c] for c in ("Px", "Py", "Pz"))
            runs[(src, int(r["N"]), grid, r["pack"])].append(float(r["t_step"]) * 1e3)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    w = 0.2
    for ax, N in zip(axes, [240, 960]):
        for k, (src, pack, label) in enumerate(PACKS):
            xs = [g + (k - 1.5) * (w + 0.01) for g in range(len(GRIDS))]
            ys = [statistics.median(runs[(src, N, g, pack)]) for g, _ in GRIDS]
            ax.bar(xs, ys, width=w, color=COLORS[k], label=label, edgecolor="#fcfcfb", linewidth=1)
        ax.set_xticks(range(len(GRIDS)), [g[1] for g in GRIDS])
        ax.set(title=f"N = {N}")
        ax.grid(True, axis="y")
    axes[0].set_ylabel("time per step (ms)")
    axes[0].legend(frameon=False, loc="upper right", fontsize=9)
    fig.suptitle("Packing with index lists (P = 48, nonblocking)", fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/E2_indexed.png", dpi=150)


if __name__ == "__main__":
    main()
