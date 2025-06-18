
from glob import glob

import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style='ticks', context='paper', palette='colorblind')
import cmocean

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

def add_map_features(ax):
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND)
    gl = ax.gridlines(draw_labels=True, zorder=1)
    gl.top_labels=False   # suppress top labels
    gl.right_labels=False # suppress right labels

def load_bathymetry(zip_file_url):
    """Read zip file from Natural Earth containing bathymetry shapefiles"""
    # Download and extract shapefiles
    import io
    import zipfile

    import requests
    r = requests.get(zip_file_url)
    z = zipfile.ZipFile(io.BytesIO(r.content))
    z.extractall("ne_10m_bathymetry_all/")

    # Read shapefiles, sorted by depth
    shp_dict = {}
    files = glob('ne_10m_bathymetry_all/*.shp')
    assert len(files) > 0
    files.sort()
    depths = []
    for f in files:
        depth = '-' + f.split('_')[-1].split('.')[0]  # depth from file name
        depths.append(depth)
        bbox = (-67, 41, -56, 48)  # (x0, y0, x1, y1)
        nei = shpreader.Reader(f, bbox=bbox)
        shp_dict[depth] = nei
    depths = np.array(depths)[::-1]  # sort from surface to bottom
    return depths, shp_dict

meta = pd.read_csv('../meta/ArgoFloat_DeploymentMetadata_EN728.csv')

depths_str, shp_dict = load_bathymetry('https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_bathymetry_all.zip')
# Construct a discrete colormap with colors corresponding to each depth
depths = depths_str.astype(int)
N = len(depths)
nudge = 0.01  # shift bin edge slightly to include data
boundaries = [min(depths)] + sorted(depths+nudge)  # low to high
norm = matplotlib.colors.BoundaryNorm(boundaries, N)
depth_cm = cmocean.cm.deep.resampled(N)
colors_depths = depth_cm(norm(depths))

fig = plt.figure(constrained_layout=True)
ax = fig.add_subplot(projection=ccrs.PlateCarree())
# add features, bathy
add_map_features(ax)
for i, depth_str in enumerate(depths_str):
  ax.add_geometries(
    shp_dict[depth_str].geometries(),
    crs=ccrs.PlateCarree(),
    color=colors_depths[i], zorder=0
  )
# plot profiles so far
sns.scatterplot(
  data=meta, x='Longitude', y='Latitude', 
  hue='WMO', ax=ax, palette='colorblind',
  transform=ccrs.PlateCarree(), zorder=2,
  s=100, linewidth=1.5
)
ax.set_extent((-67, -56, 41, 48))

fig.savefig('figures/2025_spring/deployment_map.png', bbox_inches='tight', dpi=350)
plt.close(fig)