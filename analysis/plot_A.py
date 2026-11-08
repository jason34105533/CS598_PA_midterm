"""Experiment A: time per step vs halo volume for all 45 grids of P = 48.

usage: python3 analysis/plot_A.py   (from the repo root; writes figs/A_shape.png)
"""
import collections
import csv
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import INK2, SERIES  # shared style (sets rcParams on import)

MARKERS = ["o", "s", "^"]
GROUPS = [("Px = 1", lambda px: px == 1), ("Px = 2-4", lambda px: 2 <= px <= 4),
          ("Px >= 6", lambda px: px >= 6)]
LABEL = {"1x1x48", "48x1x1", "1x48x1"}


def main():
    runs = collections.defaultdict(list)
    for r in csv.DictReader(open("results/A.csv")):
        runs[(int(r["N"]), int(r["Px"]), int(r["Py"]), int(r["Pz"]))].append(r)
    Ns = sorted({k[0] for k in runs})
    fig, axes = plt.subplots(1, len(Ns), figsize=(4.2 * len(Ns), 4.2))
    for ax, N in zip(axes, Ns):
        pts = []
        for (n, px, py, pz), rs in runs.items():
            if n != N:
                continue
            t = statistics.median(float(r["t_step"]) for r in rs) * 1e3
            pts.append((float(rs[0]["halo_bytes_max"]) / 1e6, t, px, f"{px}x{py}x{pz}"))
        for i, (name, sel) in enumerate(GROUPS):
            sub = [p for p in pts if sel(p[2])]
            ax.scatter([p[0] for p in sub], [p[1] for p in sub], s=36, color=SERIES[i],
                       marker=MARKERS[i], label=name, edgecolor="#fcfcfb", linewidth=0.8, zorder=3)
        for x, y, _, g in pts:
            if g in LABEL:
                ax.annotate(g, (x, y), xytext=(5, 3), textcoords="offset points",
                            fontsize=8, color=INK2)
        ax.set(xlabel="halo per rank per step (MB)", title=f"N = {N}")
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=0)
        ax.grid(True)
    axes[0].set_ylabel("time per step (ms)")
    axes[0].legend(frameon=False, loc="lower right", title="x-direction split", title_fontsize=9)
    fig.suptitle("Experiment A: all 45 process grids of P = 48", fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/A_shape.png", dpi=150)


if __name__ == "__main__":
    main()
