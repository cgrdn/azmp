
from pathlib import Path
from netCDF4 import Dataset

import pandas as pd
import matplotlib.pyplot as plt

wmos = [4902622, 4902676]

def nc_to_dataframe(fn):
    nc = Dataset(fn)

    fdict = {}

    for v in nc.variables:
        if len(nc[v].dimensions) > 1 and nc[v].dimensions[1] == 'N_LEVELS':
            fdict[v] = nc[v][:][0,:].data
    
    df = pd.DataFrame(fdict)
    return df

for wmo in wmos:
    fn = Path(f'data/{wmo}/profiles/SR{wmo}_001.nc')
    df = nc_to_dataframe(fn)
