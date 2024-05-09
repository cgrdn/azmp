#!/usr/bin/python

from pathlib import Path
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', palette='colorblind')

f_down = Path('/Users/GordonC/Documents/cruises/AZMP_FALL_2021/D185A024_oce.csv')
df = pd.read_csv(f_down)

fig, axes = plt.subplots(1, 5, sharey=True)
df['oxygenError'] = df.oxygen - df.oxygen2
df['oxygenRawError'] = df.oxygenRaw - df.oxygenRaw2
df['oxygenPctError'] = df.oxygenError/df.oxygen*100

vars_to_plot = ['oxygen', 'oxygenRaw', 'oxygenError', 'oxygenRawError', 'oxygenPctError']

for k, ax in zip(vars_to_plot, axes):
    df.plot(k, 'pressure', ax=ax, legend=False)
axes[0].set_ylim((1800, 0))

w, h = fig.get_figwidth(), fig.get_figheight()
fig.set_size_inches(w*5/3, h)
fig.savefig(Path('figures/oxygen_profiles.png'), dpi=350, bbox_inches='tight')

fig, axes = plt.subplots(2, 2)
axes = axes.flatten()
df.plot('oxygen', 'oxygenError', c='pressure', kind='scatter', cmap='viridis', ax=axes[0])
df.plot('oxygenRaw', 'oxygenRawError', c='pressure', kind='scatter', cmap='viridis', ax=axes[1])
df.plot('oxygenRawError', 'oxygenError', c='pressure', kind='scatter', cmap='viridis', ax=axes[2])
df.plot('oxygenRaw', 'oxygen', c='oxygenError', kind='scatter', cmap='viridis', ax=axes[3])

fig.tight_layout()
fig.savefig(Path('../figures/oxygen_scatters.png'), dpi=350, bbox_inches='tight')

plt.show()