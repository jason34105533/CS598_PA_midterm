"""Experiment F: for prime P, 1x1xP vs the layered irregular decomposition vs a
regular grid with fewer ranks (the next non-prime below P).

usage: PYTHONPATH=analysis python3 analysis/plot_F.py   (writes figs/F_irregular.png)
"""
import collections
import csv
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import SERIES

PRIMES = [37, 41, 43, 47]
FEWER = {37: 36, 41: 40, 43: 42, 47: 45}


def load(path):
    runs = collections.defaultdict(list)
    for r in csv.DictReader(open(path)):
        runs[(int(r["N"]), int(r["P"]))].append(float(r["t_step"]) * 1e3)
    return {k: statistics.median(v) for k, v in runs.items()}


def main():
    irr, ref = load("results/F.csv"), load("results/F_ref.csv")
    Ns = sorted({k[0] for k in irr})
    fig, axes = plt.subplots(1, len(Ns), figsize=(5.6 * len(Ns), 4.2))
    w = 0.26
    for ax, N in zip(axes, Ns):
        bars = [("1x1xP", lambda P: ref[(N, P)]),
                ("layered, P ranks", lambda P: irr[(N, P)]),
                ("regular grid, fewer ranks", lambda P: ref[(N, FEWER[P])])]
        for k, (label, f) in enumerate(bars):
            xs = [i + (k - 1) * (w + 0.02) for i in range(len(PRIMES))]
            ax.bar(xs, [f(P) for P in PRIMES], width=w, color=SERIES[k], label=label,
                   edgecolor="#fcfcfb", linewidth=1)
        ax.set_xticks(range(len(PRIMES)), [f"P = {P}\n(fewer: {FEWER[P]})" for P in PRIMES])
        ax.set(title=f"N = {N}")
        ax.grid(True, axis="y")
    axes[0].set_ylabel("time per step (ms)")
    axes[0].legend(frameon=False, loc="upper left", fontsize=9)
    fig.suptitle("Experiment F: prime P with a layered irregular decomposition", fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/F_irregular.png", dpi=150)


if __name__ == "__main__":
    main()
