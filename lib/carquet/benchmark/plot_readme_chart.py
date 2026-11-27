#!/usr/bin/env python3
"""Render the README benchmark figure (res/img/bench_10m_m3_{light,dark}.png).

Usage: python3 benchmark/plot_readme_chart.py benchmark/bench_<ver>_<os>_<arch>_<date>.json [--version 0.7.3]
                                             [--name bench_10m_m3]

Draws every library present in the JSON (carquet, Arrow C++, PyArrow, fastparquet, polars).
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator, NullLocator

CODECS = ["none", "snappy", "zstd", "lz4"]
LIBS = [("carquet", None), ("arrow_cpp", "Arrow C++"), ("pyarrow", "PyArrow"), ("fastparquet", "fastparquet"),
        ("polars", "polars")]
THEMES = {
    "light": dict(bg="#fcfcfb", grid="#e6e5df", text="#1a1a19", muted="#8a897f", tick="#5d5c55",
                  bars=["#2a78d6", "#eb6834", "#1baf7a", "#8e5cd9", "#c48a12"]),
    "dark": dict(bg="#1a1a19", grid="#33332f", text="#f0efec", muted="#8f8e85", tick="#c3c2b7",
                 bars=["#3987e5", "#d95926", "#199e70", "#9d6fe6", "#dba22e"]),
}


def label(v):
    if v >= 100:
        return f"{v:.0f}"
    if v >= 1:
        return f"{v:.1f}"
    return f"{v:.2f}"


def render(data, version, theme, out):
    t = THEMES[theme]
    # dataset == "large" keeps the opt-in large_nullable configs (also 10M rows) out of the chart
    rows = {(r["library"], r["compression"]): r for r in data["results"]
            if r["rows"] == 10_000_000 and r.get("dataset", "large") == "large"}
    libs = [(i, lib, f"carquet {version}" if name is None else f"{name} {data[lib + '_version']}")
            for i, (lib, name) in enumerate(LIBS) if (lib, CODECS[0]) in rows]
    colors = [t["bars"][i] for i, _, _ in libs]
    step = (0.732 if len(libs) < 5 else 0.86) / len(libs)
    value_fs = 7.5 if len(libs) < 5 else 6.2
    iters = data["config"]["bench_iterations"]["large"]

    plt.rcParams["font.family"] = "DejaVu Sans"
    fig = plt.figure(figsize=(10.5, 4.7), dpi=200, facecolor=t["bg"])
    fig.text(0.06, 0.925, f"Parquet at 10 million rows (INT64 + DOUBLE + INT32), {data['system']['cpu']}, {data['timestamp']}",
             fontsize=12, fontweight="bold", color=t["text"], va="baseline")
    fig.text(0.06, 0.87, f"Trimmed median of {iters} runs, page cache purged between write and read. "
             "Uncompressed reads use carquet's mmap zero-copy path.", fontsize=8.5, color=t["muted"], va="baseline")
    fig.legend(handles=[Patch(color=c) for c in colors], labels=[n for _, _, n in libs], loc="upper left",
               bbox_to_anchor=(0.055, 0.855), ncol=len(libs), frameon=False, fontsize=9.5, labelcolor=t["text"],
               handlelength=1.3, handleheight=1.0, columnspacing=1.6, handletextpad=0.9)

    for i, (key, title, log) in enumerate([("write_ms", "Write time, ms (lower is better)", False),
                                           ("read_ms", "Read time, ms, log scale (lower is better)", True)]):
        ax = fig.add_axes([0.06 + i * 0.5105, 0.1207, 0.4185, 0.565], facecolor=t["bg"])
        top = 0
        for j, (_, lib, _) in enumerate(libs):
            vals = [rows[(lib, c)][key] for c in CODECS]
            top = max(top, *vals)
            xs = [k + (j - (len(libs) - 1) / 2) * step for k in range(len(CODECS))]
            ax.bar(xs, vals, width=step * 0.9, color=colors[j], zorder=3)
            for x, v in zip(xs, vals):
                ax.annotate(label(v), (x, v), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom",
                            fontsize=value_fs, color=t["tick"])
        if log:
            ax.set_yscale("log")
            ax.set_ylim(0.1, max(top * 1.6, 120))
            ax.yaxis.set_major_locator(FixedLocator([0.1, 1, 10, 100, 1000]))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            ax.yaxis.set_minor_locator(NullLocator())
        else:
            ax.set_ylim(0, top * 1.12)
            ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.set_xlim(-0.55, len(CODECS) - 0.45)
        ax.set_xticks(range(len(CODECS)), CODECS)
        ax.grid(axis="y", color=t["grid"], linewidth=1, zorder=0)
        ax.tick_params(length=0, colors=t["tick"], labelsize=10, pad=5)
        ax.tick_params(axis="x", labelsize=10.5)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(title, loc="left", fontsize=10, color=t["text"], pad=12)

    fig.savefig(out, dpi=200, facecolor=t["bg"])
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("json")
    ap.add_argument("--version", help="carquet version shown in the legend (default: from the JSON)")
    ap.add_argument("--name", default="bench_10m_m3", help="output basename, _{light,dark}.png is appended")
    ap.add_argument("--out-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "res", "img"))
    args = ap.parse_args()
    data = json.load(open(args.json))
    for theme in THEMES:
        out = os.path.join(args.out_dir, f"{args.name}_{theme}.png")
        render(data, args.version or data["carquet_version"], theme, out)
        print(out)


if __name__ == "__main__":
    main()
