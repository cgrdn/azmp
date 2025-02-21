
import numpy as np
import pandas as pd
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
sns.set_theme(style='ticks', context='paper', palette='colorblind')
import cmocean.cm as cmo

import bgcArgoDMQC as bgc
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


meta = pd.read_csv('../meta/ArgoFloat_DeploymentMetadata_DY18402.csv')
ship = pd.read_csv('data/DY18402_extracted20250212.csv')

fig, axes = plt.subplots(2, 5, sharey=True)

varcolors = {
    'TEMP':cmo.thermal(0.8),
    'PSAL':cmo.haline(0.3),
    'DOXY':sns.color_palette('Blues')[5],
    'CHLA':sns.color_palette('Greens')[4],
    'BBP700':'brown',
}

varlabels = {
    'TEMP':f'Temperature ({chr(176)}C)', 
    'PSAL':'Salinity', 
    'DOXY':'Dissolved\nOxygen (mL L$^{-1}$)', 
    'CHLA':'Chlorophyll\n(mg m$^{-3}$)', 
    'BBP700':'Backscatter\n(x1000 m$^{-1}$)'
}

varlims = {
    'TEMP':(1.7, 19.5), 
    'PSAL':(32.05, 35.4), 
    'DOXY':(3, 7.5), 
    'CHLA':(-0.1, 1.45), 
    'BBP700':(-0.2, 6.4)
}

figlabels = [['a', 'b', 'c', 'd', 'e'], ['f', 'g', 'h', 'i', 'j']]

gain = {wmo:pd.Series(bgc.sprof(wmo).calc_gains(ref='WOA')).mean() for wmo in meta.WMO}

argo_vars = ['TEMP', 'PSAL', 'DOXY', 'CHLA', 'BBP700']
bottle_vars = ['', 'Salinity_Sal_PSS | Salinity', 'O2_Winkler_Auto | O2', 'Chl_a_Holm-Hansen_F | Chlorophyll A', '']
ctd_vars = ['', 'Salinity_CTD | Salinity', 'O2_CTD_mLL | O2', 'Chl_Fluor_Voltage | CHL-SENSOR insitu', '']

bottom = False

for axrow, wmo, stn, figs in zip(axes, meta.WMO, meta.Station, figlabels):
    argo = nc_to_dataframe(f'data/{wmo}/profiles/SR{wmo}_001.nc', wmo=wmo)
    argo['rho0'] = gsw.pot_rho_t_exact(
        gsw.SA_from_SP(argo.PSAL.values, argo.PRES.values, -58, 43),
        argo.TEMP.values, argo.PRES.values, 0
    )
    argo['BBP700'] = argo['BBP700']*1000
    argo['DOXY'] = argo['DOXY']*gain[wmo]
    argo['DOXY'] = bgc.unit.umol_per_L_to_mL_per_L(argo.DOXY*argo.rho0/1000, argo.TEMP)
    axrow[0].set_ylabel('Pressure (dbar)', fontsize=8)
    axrow[0].set_title(f'{stn} - {wmo}', loc='left', fontweight='bold', fontsize=10)
    for ax, av, bv, cv, f in zip(axrow, argo_vars, bottle_vars, ctd_vars, figs):
        nearest_station = 'LL_09' if stn == 'M2300' else stn
        sl = argo[av].notna()
        ax.plot(
            argo.loc[sl, av], 
            argo.loc[sl, 'PRES'], 
            color=varcolors[av], linewidth=1.25, zorder=1
        )
        if bv:
            ax.scatter(
                ship.loc[ship['STATION'] == nearest_station, bv], 
                ship.loc[ship['STATION'] == nearest_station, 'Pressure | Pressure'], 
                marker='s', color=varcolors[av], linewidth=1, edgecolor='k', zorder=2
            )
        # if cv:
        #     rl = ship[cv].notna()
        #     sub = ship.loc[rl].sort_values('Pressure | Pressure')
        #     ax.plot(
        #         sub[cv], 
        #         sub['Pressure | Pressure'],
        #         color='k', linewidth=0.5, zorder=0
        #     )
        if av == 'CHLA':
            inset = inset_axes(ax, width='50%', height='60%', loc=7)
            inset.plot(argo.loc[(sl) & (argo.PRES < 200), av], argo.loc[(sl) & (argo.PRES < 200), 'PRES'], color=varcolors[av], zorder=1)
            inset.set_ylim((155, -5))
            inset.scatter(ship.loc[ship['STATION'] == nearest_station, bv], ship.loc[ship['STATION'] == nearest_station, 'Pressure | Pressure'], marker='s', color=varcolors[av], linewidth=1, edgecolor='k', zorder=2)
            inset.tick_params(axis='both', which='major', labelsize=6)
            inset.set_yticks((0, 75, 150))
        at = AnchoredText(f'({f})', loc='lower left', prop=dict(size=8, fontweight='bold'), frameon=False)
        ax.add_artist(at)
        if bottom:
            ax.set_xlabel(varlabels[av], fontsize=8)
        else:
            ax.set_xticklabels([])
        ax.set_xlim(varlims[av])
    ax.set_ylim((2050,-50))
    inset.tick_params(axis='both', which='major', labelsize=8)
    bottom = True

fig.savefig('figures/2024_fall/azmp_2024_fall_CTD_bottle_first_profiles.png', dpi=350, bbox_inches='tight')
plt.close(fig)