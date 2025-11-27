import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

filename = sys.argv[1] if len(sys.argv) > 1 else 'vanderpol.csv'
data = pd.read_csv(filename)

sns.set_theme()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(data['t'], data['x'])
ax1.set_xlabel('Time')
ax1.set_ylabel('Position')

ax2.plot(data['x'], data['v'])
ax2.set_xlabel('Position')
ax2.set_ylabel('Velocity')

plt.tight_layout()
plt.savefig('result.png')
plt.show()
