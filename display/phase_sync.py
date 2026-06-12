#!/usr/bin/env uv run python

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

# Коды режимов: бит 0 = s01≈0, бит 1 = s02≈0, бит 2 = s12≈0
REGIME_LABELS = [
    "нет синхр.",     # 0: 000
    "0↔1",            # 1: 001
    "0↔2",            # 2: 010
    "0↔1, 0↔2",      # 3: 011
    "1↔2",            # 4: 100
    "0↔1, 1↔2",      # 5: 101
    "0↔2, 1↔2",      # 6: 110
    "полная синхр.",  # 7: 111
]

REGIME_COLORS = [
    "#d3d3d3",  # 0 — нет синхр.
    "#4477aa",  # 1
    "#66ccee",  # 2
    "#228833",  # 3
    "#ccbb44",  # 4
    "#ee6677",  # 5
    "#aa3377",  # 6
    "#222222",  # 7 — полная синхр.
]


def sync_code(s01, s02, s12, theta):
    return int(s01 < theta) | (int(s02 < theta) << 1) | (int(s12 < theta) << 2)


def to_grid(df, col):
    pivot = df.pivot_table(index="delta1", columns="delta2", values=col).sort_index()
    d1 = pivot.index.values
    d2 = pivot.columns.values
    return d1, d2, pivot.values


def cell_edges(v):
    step = v[1] - v[0]
    return np.append(v - step / 2, v[-1] + step / 2)


def plot_eps(ax_coh, ax_sync, sub_df, eps, theta, vmin, vmax):
    d1, d2, L_grid = to_grid(sub_df, "L")

    sub_df = sub_df.copy()
    sub_df["code"] = sub_df.apply(
        lambda r: sync_code(r["s01"], r["s02"], r["s12"], theta), axis=1
    )
    _, _, code_grid = to_grid(sub_df, "code")

    D2, D1 = np.meshgrid(cell_edges(d2), cell_edges(d1))

    # --- когерентность ---
    pcm = ax_coh.pcolormesh(D2, D1, L_grid, cmap="YlOrRd", vmin=vmin, vmax=vmax, shading="flat")
    plt.colorbar(pcm, ax=ax_coh, label="L")
    ax_coh.set_title(f"ε = {eps:.3f}  —  когерентность L")
    ax_coh.set_xlabel("$\\delta_2$")
    ax_coh.set_ylabel("$\\delta_1$")

    # --- режимы синхронизации ---
    cmap_disc = mcolors.ListedColormap(REGIME_COLORS)
    norm_disc = mcolors.BoundaryNorm(np.arange(-0.5, 8.5), ncolors=8)
    ax_sync.pcolormesh(D2, D1, code_grid, cmap=cmap_disc, norm=norm_disc, shading="flat")
    ax_sync.set_title(f"ε = {eps:.3f}  —  режимы синхронизации  (θ = {theta})")
    ax_sync.set_xlabel("$\\delta_2$")
    ax_sync.set_ylabel("$\\delta_1$")

    present = sorted(sub_df["code"].unique())
    ax_sync.legend(
        handles=[Patch(facecolor=REGIME_COLORS[c], label=REGIME_LABELS[c]) for c in present],
        loc="upper right",
        fontsize=7,
    )


def main():
    parser = argparse.ArgumentParser(description="Phase synchronization regimes visualization")
    parser.add_argument("file", nargs="?", default="results.csv")
    parser.add_argument(
        "--theta", type=float, default=0.05,
        help="Порог для 'наклон ≈ 0' (default: 0.05)"
    )
    parser.add_argument("--show", action="store_true", help="Открыть интерактивное окно")
    args = parser.parse_args()

    df = pd.read_csv(args.file, comment="#").replace("nan", np.nan).dropna()
    eps_values = sorted(df["eps"].unique())
    vmin, vmax = df["L"].min(), df["L"].max()

    n = len(eps_values)
    fig, axes = plt.subplots(n, 2, figsize=(12, 5 * n))
    if n == 1:
        axes = axes[np.newaxis, :]

    for i, eps in enumerate(eps_values):
        plot_eps(axes[i, 0], axes[i, 1], df[df["eps"] == eps], eps, args.theta, vmin, vmax)

    plt.tight_layout()
    out = args.file.replace(".csv", "_phase_sync.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
