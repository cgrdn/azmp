
import numpy as np
import pandas as pd
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sns
sns.set_theme(style='ticks', context='paper', palette='colorblind')
import cmocean.cm as cmo

import gsw

def nc_to_dataframe(fn, wmo=None):
    nc = Dataset(fn)

    fdict = {}

    for v in nc.variables:
        if len(nc[v].dimensions) > 1 and nc[v].dimensions[1] == 'N_LEVELS':
            fdict[v] = nc[v][:][0,:].data
            if v[-2:] == 'QC':
                fdict[v] = [f.decode() for f in fdict[v]]
            fdict[v][fdict[v] == nc[v]._FillValue] = np.nan
    
    fdict['wmo'] = fdict['PRES'].shape[0]*[wmo]

    df = pd.DataFrame(fdict)
    return df


meta = pd.read_csv('../meta/ArgoFloat_DeploymentMetadata_EN728.csv')

fig, axes = plt.subplots(2, 2, sharey=True)

varcolors = {
    'TEMP':cmo.thermal(0.8),
    'PSAL':cmo.haline(0.3),
}

varlabels = {
    'TEMP':f'Temperature ({chr(176)}C)', 
    'PSAL':'Salinity', 
}

varlims = {
    'TEMP':(1.7, 13.1), 
    'PSAL':(32.05, 35.4), 
}

figlabels = [['a', 'b'], ['c', 'd']]

argo_vars = ['TEMP', 'PSAL']
bottom = False

for axrow, wmo, stn, figs in zip(axes, meta.WMO, meta.Station, figlabels):
    argo = nc_to_dataframe(f'data/{wmo}/profiles/R{wmo}_001.nc', wmo=wmo)
    argo['rho0'] = gsw.pot_rho_t_exact(
        gsw.SA_from_SP(argo.PSAL.values, argo.PRES.values, -58, 43),
        argo.TEMP.values, argo.PRES.values, 0
    )
    axrow[0].set_ylabel('Pressure (dbar)', fontsize=8)
    axrow[0].set_title(f'{stn} - {wmo}', loc='left', fontweight='bold', fontsize=10)
    for ax, av, f in zip(axrow, argo_vars, figs):
        sl = argo[av].notna()
        ax.plot(
            argo.loc[sl, av], 
            argo.loc[sl, 'PRES'], 
            color=varcolors[av], linewidth=1.25, zorder=1
        )
        at = AnchoredText(f'({f})', loc='lower left', prop=dict(size=8, fontweight='bold'), frameon=False)
        ax.add_artist(at)
        if bottom:
            ax.set_xlabel(varlabels[av], fontsize=8)
        else:
            ax.set_xticklabels([])
        ax.set_xlim(varlims[av])
    ax.set_ylim((2050,-50))
    bottom = True

fig.set_size_inches((fig.get_figwidth()*2/5, fig.get_figheight()))

fig.savefig('figures/2025_spring/azmp_2025_spring_first_profiles.png', dpi=350, bbox_inches='tight')
plt.close(fig)