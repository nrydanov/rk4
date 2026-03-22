#!/usr/bin/env uv run python

import sys
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

def plot_heatmaps(csv_file: str):
    df = pd.read_csv(csv_file).replace('nan', np.nan)
    eps_values = sorted(df['eps'].dropna().unique())
    n = len(eps_values)

    cmap = sns.color_palette("RdBu_r", as_cmap=True)
    cmap.set_bad(color='lightgray')
    vmin, vmax = df['L'].min(), df['L'].max()

    cols = 2
    fig, axes = plt.subplots(math.ceil(n / cols), cols, figsize=(12, 5.5 * math.ceil(n / cols)))
    axes = axes.flatten()

    for i, eps in enumerate(eps_values):
        pivot = (df[df['eps'] == eps]
                 .pivot_table(index='delta1', columns='delta2', values='L')
                 .sort_index(ascending=False))
        sns.heatmap(pivot, ax=axes[i], cmap=cmap, vmin=vmin, vmax=vmax,
                    center=0, mask=pivot.isna(), square=True, cbar=False)
        axes[i].set_title(f'ε = {eps:.3f}')
        axes[i].set_xlabel('$\\delta_2$')
        axes[i].set_ylabel('$\\delta_1$', rotation=0, labelpad=20)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    fig.colorbar(ScalarMappable(Normalize(vmin, vmax), cmap), ax=axes[:n], shrink=0.6, label='L')
    plt.savefig(csv_file.replace('.csv', '.png'), dpi=150, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    plot_heatmaps(sys.argv[1])
