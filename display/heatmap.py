#!/usr/bin/env uv run python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sns.set_theme(style="whitegrid", palette="deep", font_scale=1.1)

def plot_heatmaps(csv_file):
    df = pd.read_csv(csv_file).replace('nan', np.nan)
    eps_values = sorted(df['eps'].dropna().unique())

    # Глобальные min/max ПО EPS
    vlims = {}
    for eps in eps_values:
        eps_data = df[(df['eps'] == eps) & df['L'].notna()]['L']
        vlims[eps] = (eps_data.min(), eps_data.max())

    # Seaborn FacetGrid
    g = sns.FacetGrid(df.dropna(subset=['L']), col='eps', col_wrap=2,
                     col_order=eps_values, height=5, aspect=1,
                     sharex=True, sharey=True)

    def plot_heatmap(data, **kws):
        eps = data['eps'].iloc[0]  # текущий eps
        vmin, vmax = vlims[eps]

        pivot = data.pivot(index='delta1', columns='delta2', values='L')
        sns.heatmap(pivot, cmap='RdBu_r', center=0,
                   vmin=vmin, vmax=vmax,  # ГЛОБАЛЬНЫЕ для eps!
                   cbar_kws={'shrink': 0.8, 'label': f'L (min={vmin:.3f}, max={vmax:.3f})'},
                   square=True)

    g.map_dataframe(plot_heatmap)
    g.set_titles('ε = {col_name:.3f}')
    g.set_axis_labels('δ₂', 'δ₁')
    g.fig.suptitle('VdP Synchronization Heatmaps', y=1.02, fontsize=16)
    g.tight_layout()

    out_png = csv_file.replace('.csv', '_heatmaps.png')
    g.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()

    print(f"✨ {out_png}")
    for eps in eps_values:
        vmin, vmax = vlims[eps]
        valid = len(df[(df['eps'] == eps) & df['L'].notna()])
        print(f"ε={eps}: {valid} pts | L ∈ [{vmin:.4f}, {vmax:.4f}]")

if __name__ == "__main__":
    plot_heatmaps(sys.argv[1])
