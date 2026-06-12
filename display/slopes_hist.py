#!/usr/bin/env uv run python

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_hist(csv_file: str):
    df = pd.read_csv(csv_file, comment='#')
    slopes = pd.concat([df["s01"], df["s02"], df["s12"]], ignore_index=True)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4))

    # Линейная шкала — видно общую форму распределения
    axes[0].hist(slopes, bins=100, color="steelblue", edgecolor="none")
    axes[0].set_xlabel("наклон s")
    axes[0].set_ylabel("количество")
    axes[0].set_title("Распределение наклонов (линейная шкала)")

    # Логарифмическая шкала — видно провал между горбами
    axes[1].hist(slopes, bins=100, color="steelblue", edgecolor="none", log=True)
    axes[1].set_xlabel("наклон s")
    axes[1].set_ylabel("количество (лог.)")
    axes[1].set_title("Распределение наклонов (лог. шкала)")

    plt.tight_layout()
    out_path = csv_file.replace(".csv", "_hist.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: slopes_hist.py <results.csv>")
        sys.exit(1)
    plot_hist(sys.argv[1])
