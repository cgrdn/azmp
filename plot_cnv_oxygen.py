#!/usr/bin/python

from pathlib import Path
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', palette='colorblind')

import bgcArgoDMQC as bgc

f_down = Path('/Users/GordonC/Documents/cruises/AZMP_FALL_2021/D185A024_oce.csv')
f_up   = f_down.parent / 'U185A024_oce.csv'

df_down = pd.read_csv(f_down)
df_up   = pd.read_csv(f_up)

fig, paxes = plt.subplots(1, 3, sharey=True)
fig, daxes = plt.subplots(1, 3, sharey=True)
for l, u, px, dx in zip(['temperature', 'salinity', 'oxygen'], ['{}C'.format(chr(176)), 'psu', 'mL L$^{-1}$'], paxes, daxes):
    px.plot(df_down[l], df_down['pressure'], label='down')
    dx.plot(df_down[l], df_down['sigma'], label='down')
    px.plot(df_up[l], df_up['pressure'], label='up')
    dx.plot(df_up[l], df_up['sigma'], label='up')
    px.set_xlabel(f'{l} ({u})')
    dx.set_xlabel(f'{l} ({u})')
    px.set_ylim((1800,0))
    dx.set_ylim((27.8,24))

paxes[0].legend(loc=4)
daxes[0].legend(loc=4)
paxes[0].set_ylabel('pressure (dbar)')
daxes[0].set_ylabel('density (kg L$^{-1}$)')