"""Bandwidth per pair with many pairs active at once (exchange pattern, 1 MB
messages), against the effective bandwidth fitted in the A model.

usage: PYTHONPATH=analysis python3 analysis/plot_ppmulti.py   (writes figs/pingpong_multi.png)
"""
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_S import GRID, INK2, MARKERS, SERIES

SIZE = 1 << 20
LEVELS = [("other_node", "other node"), ("other_socket", "same node, other socket"),
          ("same_socket", "same socket, other CCD")]


def main():
    bw = {}
    for r in csv.DictReader(open("results/pingpong_multi.csv")):
        if r["pattern"] == "exchange" and int(r["bytes"]) == SIZE:
            bw[(r["tag"], int(r["pairs"]))] = float(r["gbps_med"])
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.axhspan(0.65, 1.3, color=GRID, zorder=0)
    ax.text(1.05, 0.97, "A model: 0.65-1.3 GB/s per rank", color=INK2, fontsize=9, va="center")
    for i, (tag, label) in enumerate(LEVELS):
        pts = sorted((p, v) for (t, p), v in bw.items() if t == tag)
        ax.plot([p for p, _ in pts], [v for _, v in pts], marker=MARKERS[i], color=SERIES[i],
                label=label)
    ax.set(xscale="log", xlabel="pairs exchanging at once", ylabel="GB/s per pair, each direction",
           ylim=(0, None), title="Bandwidth under concurrency (1 MB messages, shm)")
    ax.set_xticks([1, 2, 4, 8, 16], ["1", "2", "4", "8", "16"])
    ax.grid(True)
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig("figs/pingpong_multi.png", dpi=150)


if __name__ == "__main__":
    main()
