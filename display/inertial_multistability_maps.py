#!/usr/bin/env uv run python

import argparse
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


T_OBS = 1600 - 240
THETA = 2 * np.pi / T_OBS
COUPLING_INERTIAL = 0

REGIME_LABELS = ["нет синхр.", "0↔1", "0↔2", "1↔2", "полная"]
REGIME_COLORS = ["#eeeeee", "#4477aa", "#66ccee", "#ccbb44", "#8b1a89"]


def read_results(path: Path) -> pd.DataFrame:
    cols = ["delta1", "delta2", "eps", "coupling_type", "s01", "s02", "s12"]
    if path.suffix == ".parquet":
        return pd.read_parquet(path, columns=cols)
    return pd.read_csv(path, comment="#", usecols=cols)


def prepare_points(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = df[df["coupling_type"] == COUPLING_INERTIAL].copy()

    b01 = df["s01"].to_numpy() < THETA
    b02 = df["s02"].to_numpy() < THETA
    b12 = df["s12"].to_numpy() < THETA
    raw_code = b01.astype(np.int8) | (b02.astype(np.int8) << 1) | (b12.astype(np.int8) << 2)
    popcount = b01.astype(np.int8) + b02.astype(np.int8) + b12.astype(np.int8)
    df["code"] = np.select(
        [popcount >= 2, raw_code == 1, raw_code == 2, raw_code == 4],
        [4, 1, 2, 3],
        default=0,
    ).astype(np.int8)

    keys = ["delta1", "delta2", "eps"]
    counts = df.groupby(keys + ["code"]).size().rename("n").reset_index()
    totals = counts.groupby(keys)["n"].sum().rename("total").reset_index()
    counts = counts.merge(totals, on=keys)
    counts["prob"] = counts["n"] / counts["total"]

    dom = counts.sort_values("n").drop_duplicates(keys, keep="last")
    points = dom[keys + ["code", "n", "total"]].rename(
        columns={"code": "dominant_code", "n": "dominant_n"}
    )
    points["agree_frac"] = points["dominant_n"] / points["total"]

    n_codes = counts.groupby(keys)["code"].nunique().rename("n_codes").reset_index()
    points = points.merge(n_codes, on=keys)
    points["minority_n"] = points["total"] - points["dominant_n"]
    points["multi_strength"] = 1 - points["agree_frac"]
    points["multi"] = points["n_codes"] > 1

    probs = counts.pivot_table(
        index=keys,
        columns="code",
        values="prob",
        fill_value=0.0,
    ).reset_index()
    probs = probs.rename(columns={code: f"p{code}" for code in range(len(REGIME_LABELS))})
    for code in range(len(REGIME_LABELS)):
        col = f"p{code}"
        if col not in probs:
            probs[col] = 0.0
    points = points.merge(probs[keys + [f"p{code}" for code in range(len(REGIME_LABELS))]], on=keys)

    return points, counts


def grid_from(points: pd.DataFrame, eps: float, col: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    piv = (
        points[points["eps"] == eps]
        .pivot(index="delta1", columns="delta2", values=col)
        .sort_index()
    )
    return piv.index.to_numpy(), piv.columns.to_numpy(), piv.to_numpy()


def plot_dominant_regime(
    points: pd.DataFrame,
    eps_values: list[float],
    out_dir: Path,
    suffix: str = "all-eps",
    title_suffix: str = "",
) -> None:
    fig, axes = plt.subplots(1, len(eps_values), figsize=(4.2 * len(eps_values), 4.2), constrained_layout=True)
    if len(eps_values) == 1:
        axes = np.array([axes])

    cmap = mcolors.ListedColormap(REGIME_COLORS)
    norm = mcolors.BoundaryNorm(np.arange(-0.5, len(REGIME_LABELS) + 0.5), cmap.N)
    for ax, eps in zip(axes, eps_values):
        d1, d2, code = grid_from(points, eps, "dominant_code")

        ax.imshow(
            code,
            origin="lower",
            extent=[d2.min(), d2.max(), d1.min(), d1.max()],
            interpolation="nearest",
            aspect="equal",
            cmap=cmap,
            norm=norm,
        )
        ax.plot([d2.min(), d2.max()], [d2.min(), d2.max()], color="0.45", lw=0.7, ls="--", alpha=0.6)
        ax.set_title(f"ε = {eps:.2f}")
        ax.set_xlabel("$\\delta_2$")
        ax.set_ylabel("$\\delta_1$")

    handles = [Patch(facecolor=REGIME_COLORS[i], label=REGIME_LABELS[i]) for i in range(len(REGIME_LABELS))]
    fig.legend(handles=handles, loc="lower center", ncol=len(REGIME_LABELS), frameon=False)
    fig.suptitle(f"Inertial: доминирующий режим синхронизации{title_suffix}")
    fig.savefig(out_dir / f"inertial-dominant-regime-{suffix}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_multistability(
    points: pd.DataFrame,
    eps_values: list[float],
    out_dir: Path,
    suffix: str = "all-eps",
    title_suffix: str = "",
) -> None:
    fig, axes = plt.subplots(1, len(eps_values), figsize=(4.2 * len(eps_values), 4.2), constrained_layout=True)
    if len(eps_values) == 1:
        axes = np.array([axes])

    pcm = None
    for ax, eps in zip(axes, eps_values):
        d1, d2, strength = grid_from(points, eps, "multi_strength")
        pcm = ax.imshow(
            strength,
            origin="lower",
            extent=[d2.min(), d2.max(), d1.min(), d1.max()],
            interpolation="nearest",
            aspect="equal",
            cmap="magma_r",
            vmin=0,
            vmax=2 / 3,
        )
        ax.plot([d2.min(), d2.max()], [d2.min(), d2.max()], color="white", lw=0.8, alpha=0.55)
        ax.set_title(f"ε = {eps:.2f}")
        ax.set_xlabel("$\\delta_2$")
        ax.set_ylabel("$\\delta_1$")

    fig.colorbar(pcm, ax=axes, shrink=0.75, label="$1 - p_{mode}$")
    fig.suptitle(f"Inertial: сила мультистабильности по 30 начальным условиям{title_suffix}")
    fig.savefig(out_dir / f"inertial-multistability-strength-{suffix}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def clip_center(points: pd.DataFrame, zoom_min: float, zoom_max: float) -> pd.DataFrame:
    return points[
        points["delta1"].between(zoom_min, zoom_max)
        & points["delta2"].between(zoom_min, zoom_max)
    ].copy()


def main() -> None:
    parser = argparse.ArgumentParser(description="Inertial multistability maps for all eps")
    parser.add_argument("file", nargs="?", default="../results_new.parquet")
    parser.add_argument("--out-dir", default="../plots")
    parser.add_argument("--zoom-min", type=float, default=-0.15)
    parser.add_argument("--zoom-max", type=float, default=0.25)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = read_results(Path(args.file))
    points, counts = prepare_points(df)
    eps_values = sorted(points["eps"].unique())

    plot_dominant_regime(points, eps_values, out_dir)
    plot_multistability(points, eps_values, out_dir)

    zoom_points = clip_center(points, args.zoom_min, args.zoom_max)
    center_title = f"; центр δ∈[{args.zoom_min:g}, {args.zoom_max:g}]"
    plot_dominant_regime(zoom_points, eps_values, out_dir, suffix="center", title_suffix=center_title)
    plot_multistability(zoom_points, eps_values, out_dir, suffix="center", title_suffix=center_title)

    summary = points.groupby("eps").agg(
        cells=("multi", "size"),
        multi_cells=("multi", "sum"),
        multi_frac=("multi", "mean"),
        mean_agree=("agree_frac", "mean"),
    )
    print(f"theta = {THETA:.6g}")
    print(summary.to_string(float_format=lambda x: f"{x:.4f}"))
    print(f"Saved plots to: {out_dir}")


if __name__ == "__main__":
    main()
