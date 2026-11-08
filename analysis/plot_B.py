"""Experiment B: time per step vs P = 36..48 (auto grid); primes are forced to 1x1xP.

usage: PYTHONPATH=analysis python3 analysis/plot_B.py   (writes figs/B_badP.png)
"""
import collections
import csv
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import INK2, SERIES

PRIMES = {37, 41, 43, 47}


def main():
    runs = collections.defaultdict(list)
    for r in csv.DictReader(open("results/B.csv")):
        runs[(int(r["N"]), int(r["P"]))].append(r)
    Ns = sorted({k[0] for k in runs})
    fig, axes = plt.subplots(1, len(Ns), figsize=(5.2 * len(Ns), 4.2))
    for ax, N in zip(axes, Ns):
        Ps = sorted(p for n, p in runs if n == N)
        t = {P: statistics.median(float(r["t_step"]) for r in runs[(N, P)]) * 1e3 for P in Ps}
        grid = {P: "x".join(runs[(N, P)][0][c] for c in ("Px", "Py", "Pz")) for P in Ps}
        ax.plot(Ps, [t[P] for P in Ps], color=INK2, lw=1, zorder=1)
        other = [P for P in Ps if P not in PRIMES]
        ax.scatter(other, [t[P] for P in other], s=40, color=SERIES[0], marker="o",
                   label="non-prime (balanced grid)", zorder=3)
        ax.scatter(sorted(PRIMES), [t[P] for P in sorted(PRIMES)], s=46, color=SERIES[1],
                   marker="s", label="prime (1x1xP)", zorder=3)
        for P in (40, 45, 47, 48):
            ax.annotate(grid[P], (P, t[P]), xytext=(0, 8), textcoords="offset points",
                        ha="center", fontsize=8, color=INK2)
        ax.set(xlabel="ranks P", title=f"N = {N}")
        ax.set_ylim(0, max(t.values()) * 1.15)
        ax.set_xticks(Ps)
        ax.grid(True, axis="y")
    axes[0].set_ylabel("time per step (ms)")
    axes[0].legend(frameon=False, loc="lower left")
    fig.suptitle("Experiment B: prime P vs fewer ranks", fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/B_badP.png", dpi=150)


if __name__ == "__main__":
    main()
