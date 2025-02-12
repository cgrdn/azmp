
from pathlib import Path
from netCDF4 import Dataset

import numpy as np
import pandas as pd
import seaborn as sns
sns.set_theme(style='ticks', palette='colorblind')

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid.inset_locator import inset_axes

import pycnv
import cmocean.cm as cmo
import gsw
import bgcArgoDMQC as bgc

wmos = [4902622, 4902676]

def nc_to_dataframe(fn, wmo=None):
    nc = Dataset(fn)

    fdict = {}

    for v in nc.variables:
        if len(nc[v].dimensions) > 1 and nc[v].dimensions[1] == 'N_LEVELS':
            fdict[v] = nc[v][:][0,:].data
            fdict[v][fdict[v] == nc[v]._FillValue] = np.nan
            if v[-2:] == 'QC':
                fdict[v] = [f.decode() for f in fdict[v]]
    
    fdict['wmo'] = fdict['PRES'].shape[0]*[wmo]

    df = pd.DataFrame(fdict)
    return df


## load azmp bottle data
salinity = pd.read_csv('data/TEL2024880_Salinity_Summary.csv').set_index(['STATION', 'EVENT', 'SAMPLE_ID'])
oxygen = pd.read_csv('data/TEL2024880_Oxygen_Rpt.csv').set_index(['STATION', 'EVENT', 'SAMPLE_ID'])
chla = pd.read_csv('data/TEL2024880_Chl_Summary.csv').set_index(['STATION', 'EVENT', 'SAMPLE_ID'])
bottle = salinity.join(oxygen, lsuffix='_PSAL', rsuffix='_DOXY')
bottle = bottle.join(chla, lsuffix='', rsuffix='_CHLA')
bottle = bottle.rename(columns={'Sal_Rep1':'PSAL', 'Chl_Rep1':'CHLA1', 'Chl_Rep2':'CHLA2', 'PRESSURE':'PRESSURE_CHLA'})
bottle['DOXY'] = bgc.unit.mL_per_L_to_umol_per_L(
    bottle['Oxy_W_Rep1'], 
    bottle['Temp_CTD_P']
    )\
        /gsw.rho_t_exact(
            gsw.SA_from_SP(
                bottle['Sal_CTD_P'].values, 
                bottle['Temp_CTD_P'].values, 
                pd.Series(bottle.shape[0]*[-61.4349]).values,
                pd.Series(bottle.shape[0]*[43.4698]).values
                ), 
            bottle['Temp_CTD_P'].values, 
            bottle['PRESSURE_PSAL'].values)*1000

# load azmp CTD data
cnv_hl = pycnv.pycnv('data/D880a025.cnv')
cnv_ll = pycnv.pycnv('data/D880a151.cnv')

hl = pd.DataFrame(dict(
    PRES=cnv_hl.data['p'],
    TEMP=cnv_hl.data['T0'],
    PSAL=cnv_hl.data['sal00'],
    DOXY=bgc.unit.mL_per_L_to_umol_per_L(cnv_hl.data['sbeox0ML/L'], cnv_hl.data['T0'])/cnv_hl.cdata['pot_rho00']*1000,
    CHLA=cnv_hl.data['flSP']
))
ll = pd.DataFrame(dict(
    PRES=cnv_ll.data['p'],
    TEMP=cnv_ll.data['T0'],
    PSAL=cnv_ll.data['sal00'],
    DOXY=bgc.unit.mL_per_L_to_umol_per_L(cnv_ll.data['sbeox0ML/L'], cnv_ll.data['T0'])/cnv_ll.cdata['pot_rho00']*1000,
    CHLA=cnv_ll.data['flSP']
))
cnv = pd.concat([hl, ll], keys=['HL_07', 'LL_09'])

# create figure, plot
fig, axes = plt.subplots(2, 5, sharex=False, sharey=True)
varnames = ['TEMP', 'PSAL', 'DOXY', 'CHLA', 'BBP700']
varlabels = [
    f'Temperature ({chr(176)}C)', 
    'Salinity', 
    'Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)', 
    'Chlorophyll (mg m$^{-3}$)', 
    'Backscatter (x1000 m$^{-1}$)'
]
figlabels = [['a', 'b', 'c', 'd', 'e'], ['f', 'g', 'h', 'i', 'j']]
colors = [
    cmo.thermal(0.8),
    cmo.haline(0.3), 
    sns.color_palette('Blues')[5], 
    sns.color_palette('Greens')[4], 
    'brown'
]
lims = [
    [1.9,11.1],
    [32.1, 35.8],
    [140, 350],
    [-0.2, 7],
    [-0.1, 6],
]
stations = ['HL_07', 'LL_09']
bottom = False

gain = {'4902622':1.027789277187524, '4902676': 1.0248297142498382}

for wmo, stn, figs, row in zip(wmos, stations, figlabels, axes):
    fn = Path(f'data/{wmo}/profiles/SR{wmo}_001.nc')
    df = nc_to_dataframe(fn, wmo=wmo)
    df['BBP700'] = df['BBP700']*1000
    df['DOXY'] = df['DOXY']*gain[f'{wmo}']
    row[0].set_ylabel('Pressure (dbar)', fontsize=9)
    row[0].set_title(f'{stn} - {wmo}', loc='left', fontweight='bold', fontsize=10)
    for v, lim, c, l, let, ax in zip(varnames, lims, colors, varlabels, figs, row):
        sl = df[v].notna()
        ax.plot(df[v].loc[sl], df.PRES.loc[sl], color=c)
        if v in ['PSAL', 'DOXY']:
            ax.plot(bottle.loc[stn,:,:][v], bottle.loc[stn,:,:][f'PRESSURE_{v}'], '*', color=c, markeredgewidth=0.5, markeredgecolor='k')
        elif v == 'CHLA':
            inset = inset_axes(ax, width='50%', height='60%', loc=7)
            inset.plot(df[v].loc[(sl) & (df.PRES < 200)], df.PRES.loc[(sl) & (df.PRES < 200)], color=c)
            inset.set_ylim((155, -5))
            inset.plot(cnv.loc[cnv.PRES < 200].loc[stn][v], cnv.loc[cnv.PRES < 200].loc[stn].PRES, linewidth=0.5, color='k')
            inset.tick_params(axis='both', which='major', labelsize=8)
            inset.set_yticks((0, 75, 150))
            for i in range(2):
                ax.plot(bottle.loc[stn,:,:][f'{v}{i+1}'], bottle.loc[stn,:,:][f'PRESSURE_{v}'], '*', color=c, markeredgewidth=0.5, markeredgecolor='k')
                inset.plot(bottle.loc[stn,:,:][f'{v}{i+1}'], bottle.loc[stn,:,:][f'PRESSURE_{v}'], '*', color=c, markeredgewidth=0.5, markeredgecolor='k')
        if v != 'BBP700':
            ax.plot(cnv.loc[stn][v], cnv.loc[stn].PRES, linewidth=0.5, color='k')
        
        if bottom:
            ax.set_xlabel(l, fontsize=9)
        else:
            ax.set_xticklabels([])
        ax.set_xlim(lim)
        at = AnchoredText(f'({let})', loc='lower left', prop=dict(size=8, fontweight='bold'), frameon=False)
        ax.add_artist(at)
    bottom = True
    
ax.set_ylim((1525, -25))
fig.tight_layout()
fig.set_size_inches(fig.get_figwidth()*15/12, fig.get_figheight()*15/12)
fig.savefig('figures/spring_initial_profiles.png', dpi=350, bbox_inches='tight')
plt.close(fig)

# calculate oxygen gain
for wmo in wmos:
    flt = bgc.sprof(wmo)
    flt.calc_gains(ref='WOA')
    print(f'{wmo} gain: {flt.gain}')

btl = bottle.loc[stations][['PRESSURE_PSAL', 'CHLA1', 'CHLA2', 'DOXY', 'PSAL']].melt(id_vars=['PRESSURE_PSAL'], value_name='Bottle', ignore_index=False)
btl = btl.loc[btl.Bottle.notna()]

argo_match = []
pres_diff = []

df = pd.concat([
    nc_to_dataframe(Path(f'data/{wmos[0]}/profiles/SR{wmos[0]}_001.nc'), wmo=wmos[0]), 
    nc_to_dataframe(Path(f'data/{wmos[1]}/profiles/SR{wmos[1]}_001.nc'), wmo=wmos[1])
]).set_index('wmo')

stn_dict = {stn:str(w) for w, stn in zip(wmos, stations)}

for i, p, v in zip(btl.index, btl.PRESSURE_PSAL, btl.variable):
    u = v if v in ['PSAL', 'DOXY'] else 'CHLA'
    stn = i[0]
    sub = df.loc[df[u].notna()].loc[stn_dict[f'{stn}']]
    ix = (sub.PRES - p).abs() == (sub.PRES - p).abs().min()
    argo_match.append(sub.loc[ix][u].iloc[0])
    pres_diff.append((sub.PRES - p).abs().min())

btl['Argo'] = argo_match
btl['PRES_DIFF'] = pres_diff
btl['VARNAME'] = [v.replace('1', '').replace('2','') for v in btl.variable]
btl = btl.loc[btl.PRES_DIFF < 5]
btl['STATION'] = [i[0] for i in btl.index]

g = sns.FacetGrid(btl, col='VARNAME', hue='STATION', sharex=False, sharey=False, despine=False, col_order=['PSAL', 'DOXY', 'CHLA'])
g = g.map(plt.plot, 'Argo', 'Bottle', marker='*', linewidth=0, markeredgewidth=0.25, markeredgecolor='k', markersize=10)
g.axes[0,0].title.set_text(varlabels[1])
g.axes[0,1].title.set_text(varlabels[2])
g.axes[0,2].title.set_text(varlabels[3])
for i, let in enumerate(['a', 'b', 'c']):
    ax = g.axes[0,i]
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    lims = [min(xlims+ylims), max(xlims+ylims)]
    ax.plot(lims, lims, 'k', linewidth=2, zorder=-1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel(ax.get_xlabel(), fontsize=10)
    ax.set_ylabel(ax.get_ylabel(), fontsize=10)
    at = AnchoredText(f'({let})', loc='lower right', prop=dict(size=8, fontweight='bold'), frameon=False)
    ax.add_artist(at)
g.axes[0,0].legend(loc=2, fontsize=10)

g.figure.savefig('figures/spring_initial_profiles_bottle_comparison.png', dpi=350, bbox_inches='tight')
plt.close(g.figure)