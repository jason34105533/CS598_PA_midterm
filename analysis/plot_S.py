"""Experiment S: strong and weak scaling plots, and n_0.8.

usage: python3 analysis/plot_S.py   (from the repo root; writes figs/S_*.png)
"""
import collections
import csv
import math
import os
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]  # categorical slots 1-3, fixed order
MARKERS = ["o", "s", "^"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "grid.color": GRID, "grid.linewidth": 0.8,
    "font.size": 11, "lines.linewidth": 2, "lines.markersize": 6,
})


def load(path="results/S.csv"):
    t = collections.defaultdict(list)
    for r in csv.DictReader(open(path)):
        t[(r["tag"], int(r["N"]), int(r["P"]))].append(float(r["t_step"]))
    return {k: statistics.median(v) for k, v in t.items()}


def n08(points):
    """points: (n/P, efficiency) sorted by decreasing n/P. Returns the n/P where the
    efficiency first falls below 0.8, interpolated in log(n/P)."""
    for (x0, e0), (x1, e1) in zip(points, points[1:]):
        if e0 >= 0.8 > e1:
            f = (e0 - 0.8) / (e0 - e1)
            return math.exp(math.log(x0) + f * (math.log(x1) - math.log(x0)))
    return None


def main():
    med = load()
    os.makedirs("figs", exist_ok=True)
    Ns = sorted({k[1] for k in med if k[0] == "strong"})

    # strong scaling: time vs P
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for i, N in enumerate(Ns):
        Ps = sorted(k[2] for k in med if k[0] == "strong" and k[1] == N)
        ts = [med[("strong", N, P)] * 1e3 for P in Ps]
        ax.plot(Ps, ts, marker=MARKERS[i], color=SERIES[i], label=f"N = {N}")
        ax.plot(Ps, [ts[0] / P for P in Ps], ls=":", lw=1, color=SERIES[i])
        ax.annotate(f"N = {N}", (Ps[-1], ts[-1]), xytext=(6, 0), textcoords="offset points",
                    va="center", color=INK2, fontsize=10)
    ax.set(xscale="log", yscale="log", xlabel="ranks P", ylabel="time per step (ms)",
           title="Strong scaling (dotted: ideal)")
    ax.set_xticks([1, 2, 4, 8, 16, 48], labels=["1", "2", "4", "8", "16", "48"])
    ax.grid(True, which="major")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig("figs/S_strong_time.png", dpi=150)

    # strong scaling: efficiency vs n/P
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    print("n_0.8 (cells per rank where strong-scaling efficiency drops below 0.8):")
    for i, N in enumerate(Ns):
        Ps = sorted(k[2] for k in med if k[0] == "strong" and k[1] == N)
        t1 = med[("strong", N, 1)]
        pts = [(N**3 / P, t1 / (P * med[("strong", N, P)])) for P in Ps]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker=MARKERS[i], color=SERIES[i],
                label=f"N = {N}")
        x = n08(pts)
        print(f"  N = {N}: {x:.3g}" if x else f"  N = {N}: not reached")
    ax.axhline(0.8, color=INK2, lw=1, ls="--")
    ax.text(6e8, 0.82, "80%", color=INK2, fontsize=10)
    ax.set(xscale="log", xlabel="cells per rank  n/P", ylabel="parallel efficiency",
           title="Strong scaling efficiency", ylim=(0, 1.1))
    ax.invert_xaxis()
    ax.grid(True, which="major")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig("figs/S_strong_eff.png", dpi=150)

    # weak scaling: time per step vs P, normalized to P = 1
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for i, L in enumerate([60, 240]):
        tag = f"weak{L}"
        ks = sorted((k for k in med if k[0] == tag), key=lambda k: k[2])
        Ps = [k[2] for k in ks]
        eff = [med[ks[0]] / med[k] for k in ks]
        ax.plot(Ps, eff, marker=MARKERS[i], color=SERIES[i], label=f"local about {L}^3")
        ax.annotate(f"{L}^3", (Ps[-1], eff[-1]), xytext=(6, 0), textcoords="offset points",
                    va="center", color=INK2, fontsize=10)
    for x, label in [(8.5, "socket 0 full"), (16.5, "node 0 full")]:
        ax.axvline(x, color=INK2, lw=1, ls=":")
        ax.text(x * 0.94, 0.08, label, color=INK2, fontsize=9, rotation=90, ha="right", va="bottom")
    ax.set(xscale="log", xlabel="ranks P", ylabel="weak efficiency  t(1) / t(P)",
           title="Weak scaling", ylim=(0, 1.1))
    ax.set_xticks([1, 2, 4, 8, 16, 48], labels=["1", "2", "4", "8", "16", "48"])
    ax.grid(True, which="major")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig("figs/S_weak.png", dpi=150)


if __name__ == "__main__":
    main()
