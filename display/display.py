import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

parser = argparse.ArgumentParser()
parser.add_argument("file", nargs="?", default="results.csv")
parser.add_argument("--show", action="store_true", help="Open interactive window")
args = parser.parse_args()

filename = args.file
try:
    data = pd.read_csv(filename)
except FileNotFoundError:
    print(f"Error: File {filename} not found.")
    sys.exit(1)

x_cols = [c for c in data.columns if c.startswith('x')]
n = len(x_cols)

sns.set_theme(style="whitegrid")
palette = sns.color_palette("husl", n)

plt.figure(num="Time Domain", figsize=(16, 8), dpi=100)
for i in range(n):
    plt.plot(data['t'], data[f'x{i}'],
             label=f'G{i}',
             color=palette[i],
             linewidth=1.2,
             alpha=0.8)

plt.title('Oscillations over time', fontsize=16)
plt.xlabel('Time (s)', fontsize=14)
plt.ylabel('Amplitude (X)', fontsize=14)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(filename.replace('.csv', '_time.png'), dpi=150, bbox_inches='tight')

plt.figure(num="Phase Portrait", figsize=(10, 10), dpi=100)
for i in range(n):
    plt.plot(data[f'x{i}'], data[f'y{i}'],
             color=palette[i],
             linewidth=0.8,
             alpha=0.4)

plt.title('Phase Space trajectories', fontsize=16)
plt.xlabel('Position (X)', fontsize=14)
plt.ylabel('Velocity (Y)', fontsize=14)
plt.tight_layout()
plt.savefig(filename.replace('.csv', '_phase.png'), dpi=150, bbox_inches='tight')

if args.show:
    plt.show()
