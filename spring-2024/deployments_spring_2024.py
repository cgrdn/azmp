
from pathlib import Path
from netCDF4 import Dataset

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmocean.cm as cmo

wmos = [4902622, 4902676]

def nc_to_dataframe(fn):
    nc = Dataset(fn)

    fdict = {}

    for v in nc.variables:
        if len(nc[v].dimensions) > 1 and nc[v].dimensions[1] == 'N_LEVELS':
            fdict[v] = nc[v][:][0,:].data
            fdict[v][fdict[v] == nc[v]._FillValue] = np.nan
            if v[-2:] == 'QC':
                fdict[v] = [f.decode() for f in fdict[v]]
    
    df = pd.DataFrame(fdict)
    return df

fig, axes = plt.subplots(2, 5, sharex=False, sharey=True)
varnames = ['TEMP', 'PSAL', 'DOXY', 'CHLA', 'BBP700']
varlabels = [f'Temperature ({chr(176)}C)', 'Salinity', 'Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)', 'Chlorophyll (mg m$^{-3}$)', '700nm Backscatter (m$^{-1}$)']
colors = ['red', cmo.haline(0.7), 'blue', 'green', 'brown']
bottom = False

for wmo, row in zip(wmos, axes):
    fn = Path(f'data/{wmo}/profiles/SR{wmo}_001.nc')
    df = nc_to_dataframe(fn)
    row[0].set_ylabel('Pressure (dbar)')
    for v, c, l, ax in zip(varnames, colors, varlabels, row):
        sl = df[v].notna()
        ax.plot(df[v].loc[sl], df.PRES.loc[sl], color=c)
        if bottom:
            ax.set_xlabel(l)
    bottom = True
    
ax.set_ylim((1525, -25))
plt.show()
