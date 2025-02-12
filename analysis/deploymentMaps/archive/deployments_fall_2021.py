#!/usr/bin/python

from pathlib import Path
import regex as re

import numpy as np
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

# get argo global profile index
prof = argo.prof[:]

# list of floats deployed
wmos = ['4902515', '4902518', '4902519']
# get an index to only have those floats
ix = prof.shape[0]*[False]
for w in wmos:
    ix = ix | prof['file'].str.contains(w)
prof = prof[ix]
# don't want the downcasts - remove
upcasts = [re.match(r'.*D\.nc', f) is None for f in prof['file']]
prof = prof[upcasts]

# extract data
df = prof.levels[['PRES', 'PSAL', 'TEMP']]

# get density contours
vt, vs = np.meshgrid(
    np.linspace(df.TEMP.min()-2, df.TEMP.max()+2, 100),
    np.linspace(df.PSAL.min()-0.5, df.PSAL.max()+0.5, 100)
)
pden = gsw.pot_rho_t_exact(gsw.SA_from_SP(vs, 0, prof.longitude.mean(), prof.latitude.mean()), vt, 0, 0) - 1000
# contour levels
levels = list(range(18, 30, 2))

# get bathymetry for map
bath_file = Path('/Users/GordonC/Documents/data/GEBCO/GEBCO_2020.nc')
bath = Dataset(bath_file)
blat = bath['lat'][:]
blon = bath['lon'][:]
elev = bath['elevation'][:]

extent = [
    prof['longitude'].max() - 6,
    prof['longitude'].max() + 5,
    prof['latitude'].min() - 3,
    prof['latitude'].max() + 6
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
axes = [
    fig.add_subplot(gs[:,0]), 
    fig.add_subplot(gs[:,1]), 
    fig.add_subplot(gs[0, 2]), 
    fig.add_subplot(gs[1, 2], projection=ccrs.PlateCarree())
]

# order of deployment
hue_order = [
    'meds/4902518/profiles/R4902518_001.nc',
    'meds/4902515/profiles/R4902515_001.nc',
    'meds/4902519/profiles/R4902519_001.nc'
]
# temperature profile
sns.lineplot(x='TEMP', y='PRES', hue='file', hue_order=hue_order, data=df, sort=False, legend=False, ax=axes[0])
# salinity profile
sns.lineplot(x='PSAL', y='PRES', hue='file', hue_order=hue_order, data=df, sort=False, legend=False, ax=axes[1])
# TS diagram w/ pot density contours 
cs = axes[2].contour(vs, vt, pden, colors='black', levels=levels, zorder=1)
label_at_edge(levels, cs, axes[2], '%d', side='lowerleft', pad=0.1)
sns.scatterplot(x='PSAL', y='TEMP', hue='file', hue_order=hue_order, data=df, legend=False, ax=axes[2], zorder=2, edgecolor=None)
# map of float deployment locations
im = axes[3].contourf(
    blon, blat, elev,
    transform=ccrs.PlateCarree(),
    cmap=cmo.deep,
    vmin=0, extend='max'
)
cax = fig.add_axes([0.67, 0.03, 0.24, 0.015])
cb = plt.colorbar(im, ax=axes[3], orientation='horizontal', label='Depth (m)', cax=cax)
sns.scatterplot(x='longitude', y='latitude', hue='file', hue_order=hue_order, data=prof, legend=False, ax=axes[3], transform=ccrs.PlateCarree())

# format figure
axes[0].set_ylim((2000, 0))
axes[0].set_xlabel('Temperature ({}C)'.format(chr(176)))
axes[0].set_ylabel('Pressure (dbar)')
axes[0].legend(['HL_10: 4902518\n2021-09-05 05:17', 'HL_12: 4902515\n2021-09-25 18:59', 'LL_09: 4902519\n2021-09-07 20:59'], loc=4, fontsize=7)
axes[1].set_ylim((2000, 0))
axes[1].set_xlabel('Practical Salinity')
axes[1].set_ylabel('')
axes[1].set_yticklabels([])
axes[2].set_xlabel(axes[1].get_xlabel())
axes[2].set_ylabel(axes[0].get_xlabel())
axes[2].yaxis.tick_right()
axes[2].yaxis.set_label_position('right')
axes[2].xaxis.tick_top()
axes[2].xaxis.set_label_position('top')
axes[3].set_xticks(np.arange(np.ceil(extent[0]), np.floor(extent[1]), 4), crs=ccrs.PlateCarree())
axes[3].set_yticks(np.arange(np.ceil(extent[2]), np.floor(extent[3]), 4), crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
axes[3].xaxis.set_major_formatter(lon_formatter)
axes[3].yaxis.set_major_formatter(lat_formatter)
axes[3].yaxis.tick_right()
axes[3].set_xlabel('')
axes[3].set_ylabel('')
axes[3].set_extent(extent)
axes[3].add_feature(cfeature.GSHHSFeature('auto', edgecolor='black', facecolor='lightgray'))

# rotate ylabels
for ax in axes:
    for yl in ax.get_yticklabels():
        yl.update(dict(rotation=90))
        yl.set_verticalalignment('center')
for ax in axes[:3]:
    ax.grid()
tls = []
for i,cl in enumerate(cb.ax.xaxis.get_ticklabels()):
    if i % 2 == 0:
        tls.append(cl.get_text())
    else:
        tls.append('')
cb.ax.xaxis.set_ticklabels(tls)
fig.savefig(Path('../figures/azmp_floats_first_profiles.png'), bbox_inches='tight', dpi=350)
plt.close(fig)
print('Done')
