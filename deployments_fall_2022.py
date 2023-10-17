#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import Bbox
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import seaborn as sns
sns.set(style='ticks', context='paper', palette='colorblind')
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo

import argopandas as argo
import gsw

def label_at_edge(levels, cs, ax, fmt, side='both', pad=0.005, **kwargs):
    '''Label contour lines at the edge of plot
    Args:
        levels (1d array): contour levels.
        cs (QuadContourSet obj): the return value of contour() function.
        ax (Axes obj): matplotlib axis.
        fmt (str): formating string to format the label texts. E.g. '%.2f' for
            floating point values with 2 demical places.
    Keyword Args:
        side (str): on which side of the plot intersections of contour lines
            and plot boundary are checked. Could be: 'left', 'right', 'top',
            'bottom' or 'all'. E.g. 'left' means only intersections of contour
            lines and left plot boundary will be labeled. 'all' means all 4
            edges.
        pad (float): padding to add between plot edge and label text.
        **kwargs: additional keyword arguments to control texts. E.g. fontsize,
            color.
    '''
    collections = cs.collections
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    bbox = Bbox.from_bounds(xlim[0], ylim[0], xlim[1]-xlim[0], ylim[1]-ylim[0])
    eps = 1e-5  # error for checking boundary intersection
    # -----------Loop through contour levels-----------
    for ii, lii in enumerate(levels):
        cii = collections[ii]  # contours for level lii
        pathsii = cii.get_paths()  # the Paths for these contours
        if len(pathsii) == 0:
            continue
        for pjj in pathsii:
            # check first whether the contour intersects the axis boundary
            if not pjj.intersects_bbox(bbox, False):  # False significant here
                continue
            xjj = pjj.vertices[:, 0]
            yjj = pjj.vertices[:, 1]
            # intersection with the left edge
            if side in ['left', 'all', 'lowerleft', 'upperleft']:
                inter_idx = np.where(abs(xjj-xlim[0]) <= eps)[0]
                for kk in inter_idx:
                    inter_x = xjj[kk]
                    inter_y = yjj[kk]
                    ax.text(inter_x-pad, inter_y, fmt % lii,
                            ha='right',
                            va='center',
                            **kwargs)
            # intersection with the right edge
            if side in ['right', 'all', 'lowerright', 'upperright']:
                inter_idx = np.where(abs(xjj-xlim[1]) <= eps)[0]
                for kk in inter_idx:
                    inter_x = xjj[kk]
                    inter_y = yjj[kk]
                    ax.text(inter_x+pad, inter_y, fmt % lii,
                            ha='left',
                            va='center',
                            **kwargs)
            # intersection with the bottom edge
            if side in ['bottom', 'all', 'lowerleft', 'lowerright']:
                inter_idx = np.where(abs(yjj-ylim[0]) <= eps)[0]
                for kk in inter_idx:
                    inter_x = xjj[kk]
                    inter_y = yjj[kk]
                    ax.text(inter_x, inter_y-5*pad, fmt % lii,
                            ha='center',
                            va='top',
                            **kwargs)
            # intersection with the top edge
            if side in ['top', 'all', 'upperleft', 'upperright']:
                inter_idx = np.where(abs(yjj-ylim[-1]) <= eps)[0]
                for kk in inter_idx:
                    inter_x = xjj[kk]
                    inter_y = yjj[kk]
                    ax.text(inter_x, inter_y+pad, fmt % lii,
                            ha='center',
                            va='bottom',
                            **kwargs)
    return

# list of floats deployed
wmos = [4902598, 4902599]
# index of profiles
px = pd.concat([argo.float(wmo).synthetic_prof.subset_direction('asc').subset_date(date_end='20-10-2022') for wmo in wmos])
df = pd.concat([argo.float(wmo).synthetic_prof.subset_direction('asc').subset_date(date_end='20-10-2022').levels[['PRES', 'PSAL', 'TEMP', 'DOXY', 'CHLA', 'BBP700']] for wmo in wmos])

# get density contours
vt, vs = np.meshgrid(
    np.linspace(df.TEMP.min()-2, df.TEMP.max()+2, 100),
    np.linspace(df.PSAL.min()-0.5, df.PSAL.max()+0.5, 100)
)
pden = gsw.pot_rho_t_exact(gsw.SA_from_SP(vs, 0, px.longitude.mean(), px.latitude.mean()), vt, 0, 0) - 1000
# contour levels
levels = list(range(18, 30, 2))

# get bathymetry for map
bath_file = Path('/Users/GordonC/Documents/data/GEBCO/GEBCO_2020.nc')
bath = Dataset(bath_file)
blat = bath['lat'][:]
blon = bath['lon'][:]
elev = bath['elevation'][:]

extent = [
    px['longitude'].max() - 6,
    px['longitude'].max() + 5,
    px['latitude'].min() - 3,
    px['latitude'].max() + 6
]

ix = np.logical_and(blon > extent[0], blon < extent[1])
iy = np.logical_and(blat > extent[2], blat < extent[3])

blon = blon[ix]
blat = blat[iy]
elev = elev[iy,:]
elev = elev[:,ix]
elev = -np.ma.masked_array(elev.data, elev > 0)

# plot data
gs = GridSpec(2, 3)
fig = plt.figure()
primary_axes = [
    fig.add_subplot(gs[:,0]),
    fig.add_subplot(gs[:,1]), 
    fig.add_subplot(gs[0, 2]), 
    fig.add_subplot(gs[1, 2], projection=ccrs.PlateCarree())
]
axes = [
    primary_axes[0],
    primary_axes[0].twiny(),
    primary_axes[1],
    primary_axes[1].twiny(),
    primary_axes[1].twiny(),
    primary_axes[2],
    primary_axes[3]
]

# order of deployment
hue_order = [
    'meds/4902598/profiles/SR4902518_001.nc',
    'meds/4902599/profiles/SR4902515_001.nc',
]
# temperature profile
sns.lineplot(x='TEMP', y='PRES', hue='file', data=df, sort=False, legend=False, ax=axes[0], estimator=None)
# salinity profile
sns.lineplot(x='PSAL', y='PRES', hue='file', data=df, sort=False, legend=False, ax=axes[1], estimator=None)
# chla profile
sns.lineplot(x='CHLA', y='PRES', hue='file', data=df, sort=False, legend=False, ax=axes[2], estimator=None)
# oxygen profile
sns.lineplot(x='DOXY', y='PRES', hue='file', data=df, sort=False, legend=False, ax=axes[3], estimator=None)
# bbp profile
sns.lineplot(x='BBP700', y='PRES', hue='file', data=df, sort=False, legend=False, ax=axes[4], estimator=None)
# TS diagram w/ pot density contours 
cs = axes[-2].contour(vs, vt, pden, colors='black', levels=levels, zorder=1)
sns.scatterplot(x='PSAL', y='TEMP', hue='file', data=df, legend=False, ax=axes[-2], zorder=2, edgecolor=None)
# map of float deployment locations
im = axes[-1].contourf(
    blon, blat, elev,
    transform=ccrs.PlateCarree(),
    cmap=cmo.deep,
    vmin=0, extend='max'
)
cax = fig.add_axes([0.67, 0.03, 0.24, 0.015])
cb = plt.colorbar(im, ax=axes[-1], orientation='horizontal', cax=cax)
sns.scatterplot(x='longitude', y='latitude', hue='file', data=px, legend=False, ax=axes[-1], transform=ccrs.PlateCarree())

# format figure
fs = 8
label_at_edge(levels, cs, axes[5], '%d', side='lowerleft', pad=0.1, fontsize=fs)
cb.set_label('Depth (m)', fontsize=fs)
axes[0].set_ylim((2000, 0))
axes[0].set_xlabel(f'Temperature ({chr(176)}C)', loc='left', fontsize=fs)
axes[0].set_ylabel('Pressure (dbar)', fontsize=fs)
axes[0].legend(['HL_10\n4902598\n2022-10-09', 'LL_09\n4902599\n2022-10-13'], loc=3, fontsize=6)
axes[0].set_xlim((3,35))
axes[0].set_xticks([5, 10, 15])
axes[1].set_xlabel('Practical Salinity', loc='right', fontsize=fs)
axes[1].set_xlim((31, 37))
axes[1].set_xticks([35, 36])
axes[2].set_ylim((2000, 0))
axes[2].set_ylabel('')
axes[2].set_yticklabels([])
axes[2].set_xlabel('Chl $a$ (mg m$^{-3}$)', loc='left', fontsize=fs)
# axes[2].set_xlim((-0.2, 3))
axes[2].set_xlim((-0.2, 4))
# axes[2].set_xticks([0, 0.5, 1])
axes[2].set_xticks([0, 1])
axes[3].set_xlabel('Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)', fontsize=fs)
# axes[3].set_xlabel('Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)', loc='right', fontsize=fs)
axes[3].set_xlim((0, 400))
# axes[3].set_xlim((0, 270))
axes[3].set_xticks([140, 220])
# axes[3].set_xticks([140, 200, 260])
axes[4].xaxis.tick_bottom()
axes[4].xaxis.set_label_position('bottom')
axes[4].set_xlabel('$b_{bp}$ (m$^{-1}$)', loc='right', fontsize=fs)
axes[4].set_xlim((-0.013, 0.005))
axes[4].set_xticks([0, 0.004])
axes[-2].set_xlabel(axes[1].get_xlabel(), fontsize=fs)
axes[-2].set_ylabel(axes[0].get_xlabel(), fontsize=fs)
axes[-2].yaxis.tick_right()
axes[-2].yaxis.set_label_position('right')
axes[-2].xaxis.tick_top()
axes[-2].xaxis.set_label_position('top')
axes[-1].set_xticks(np.arange(np.ceil(extent[0]), np.floor(extent[1]), 4), crs=ccrs.PlateCarree())
axes[-1].set_yticks(np.arange(np.ceil(extent[2]), np.floor(extent[3]), 4), crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
axes[-1].xaxis.set_major_formatter(lon_formatter)
axes[-1].yaxis.set_major_formatter(lat_formatter)
axes[-1].yaxis.tick_right()
axes[-1].set_xlabel('')
axes[-1].set_ylabel('')
axes[-1].set_extent(extent)
axes[-1].add_feature(cfeature.GSHHSFeature('auto', edgecolor='black', facecolor='lightgray'))

# rotate ylabels
for ax in axes:
    for yl in ax.get_yticklabels():
        yl.update(dict(rotation=90, fontsize=fs))
        yl.set_verticalalignment('center')
    for xl in ax.get_xticklabels():
        xl.update(dict(fontsize=fs))
tls = []
for ax in axes[:-2]:
    ax.grid(axis='x')
axes[-2].grid()
for i,cl in enumerate(cb.ax.xaxis.get_ticklabels()):
    if i % 2 == 0:
        tls.append(cl.get_text())
    else:
        tls.append('')
cb.ax.xaxis.set_ticklabels(tls, fontsize=fs)
fig.savefig(Path('figures/azmp_fall_2022_floats_first_profiles_5var.png'), bbox_inches='tight', dpi=350)
plt.close(fig)
print('Done')
