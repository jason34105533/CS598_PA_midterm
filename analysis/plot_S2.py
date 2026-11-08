"""Experiment S2 vs S: the same P <= 16 runs with ranks filling socket 0 first
(SPREAD) or alternating between the two sockets (ALT).

usage: PYTHONPATH=analysis python3 analysis/plot_S2.py   (writes figs/S2_socket.png)
"""
import collections
import csv
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import INK2, SERIES

PS = [1, 2, 3, 4, 6, 8, 12, 16]
PANELS = [("strong", 480, "strong, N = 480"), ("strong", 960, "strong, N = 960"),
          ("weak240", None, "weak, local about 240^3")]


def main():
    t = collections.defaultdict(list)
    for r in csv.DictReader(open("results/S.csv")):
        t[(r["tag"], int(r["N"]), int(r["P"]))].append(float(r["t_step"]))
    med = {k: statistics.median(v) for k, v in t.items()}

    def get(tag, N, P):
        if N is not None:
            return med.get((tag, N, P))
        hits = [v for k, v in med.items() if k[0] == tag and k[2] == P]
        return hits[0] if hits else None

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for ax, (tag, N, title) in zip(axes, PANELS):
        t1 = get(tag, N, 1)
        for i, (suffix, label) in enumerate([("", "socket 0 first (SPREAD)"),
                                             ("_alt", "alternate sockets (ALT)")]):
            Ps, eff = [], []
            for P in PS:
                v = get(tag, N, 1) if P == 1 else get(tag + suffix, N, P)
                if v is None:
                    continue
                Ps.append(P)
                eff.append(t1 / (P * v) if tag == "strong" else t1 / v)
            ax.plot(Ps, eff, marker="os"[i], color=SERIES[i], label=label)
        ax.axvline(8.5, color=INK2, lw=1, ls=":")
        ax.set(xscale="log", xlabel="ranks P (one node)", title=title, ylim=(0, 1.2))
        ax.set_xticks(PS, [str(p) for p in PS])
        ax.grid(True)
    axes[0].set_ylabel("efficiency")
    axes[0].legend(frameon=False, loc="lower left")
    fig.suptitle("Same ranks, one CCD each: placing them on one socket or two", fontsize=12)
    fig.tight_layout()
    fig.savefig("figs/S2_socket.png", dpi=150)


if __name__ == "__main__":
    main()
