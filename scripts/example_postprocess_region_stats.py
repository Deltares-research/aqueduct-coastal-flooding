
# coding: utf-8

import numpy as np
import xarray as xr
import rasterio
import os
from os.path import join, abspath
import pandas as pd
from dask.diagnostics import ProgressBar
from multiprocessing.pool import ThreadPool
import dask
import dask.threaded
import itertools 
import sys

def get_tiles(bounds, resx, resy):
    """create regular grid with resolution (res);  returns a list with box polygons"""
    bounds = np.array(list(bounds))
    xmin, ymin, xmax, ymax = bounds
    width = int((xmax-xmin) // resx)
    height = int((ymax-ymin) // resy)
    result = []
    for i in range(width):
        for j in range(height):
            b = (xmin+i*resx, ymin+j*resy, xmin+(i+1)*resx, ymin+(j+1)*resy)
            result.append(b)
    return result


if __name__ == "__main__":
    num_workers=4
    bounds = -180, -90, 180, 90 # w, s, e, n
    tile_size = 60, 60 # x, y [degree] -> make sure it fits the bounds
    chunks={'lat':3600, 'lon': 3600} # -> save as tiles

    outdir = r'/home/dirk/Experiments/extreme_slr/output'
    expdir = r'/home/dirk/Datasets/exposure'
    geodir = r'/home/dirk/Datasets/WRI'
    out_dim = 'geogunit_107'

    ## open files
    # fn_inun = join(outdir, r'inuncoast_ext_slr_0.0_rp0000_0.nc')
    fn_inun = abspath(sys.argv[1])
    ds = xr.open_dataset(fn_inun, chunks=chunks)
    inun = ds.drop(['time']).squeeze()['inun']
    # exposure
    area = xr.open_dataset(join(expdir, r'cell_size_km2.nc'), chunks=chunks).drop(['time']).squeeze()['data']
    area.name = 'area'
    area['lon'] = ds['lon']
    area['lat'] = ds['lat']
    gdp = xr.open_dataset(join(expdir, 'base', 'exposure_base_2010_gdp.nc'), chunks=chunks)['GDP_(Regional)']
    gdp.name = 'gdp'
    gdp['lon'] = ds['lon']
    gdp['lat'] = ds['lat']
    pop = xr.open_dataset(join(expdir, 'base', 'exposure_base_2010_tpop.nc'), chunks=chunks)['Total population count']
    pop.name = 'pop'
    pop['lon'] = ds['lon']
    pop['lat'] = ds['lat']
    exposure = {'pop': pop, 'gdp': gdp, 'area':area}

    # geounits
    mask_da = xr.open_dataset(join(geodir, r'{}_all.nc'.format(out_dim)), chunks=chunks)['Geogunits'].astype(int)
    mask_da.name = out_dim
    mask_da['lon'] = ds['lon']
    mask_da['lat'] = ds['lat']
    region_df = pd.read_excel(join(geodir, r'{}_all_list.xlsx'.format(out_dim)), index_col=0, header=None)
    region_df.columns = ['ISO']
    region_df.index.name = out_dim

    tiles = get_tiles(bounds, *tile_size)
    with dask.set_options(get=dask.threaded.get, num_workers=num_workers):
        # for slr, rp in itertools.product(np.arange(0,4.1,0.5), [100]):
        #     # inun
        #     fn_inun = r'inuncoast_ext_slr_{:.1f}_rp{:04d}_0.nc'.format(slr, rp)
        #     print(fn_inun)
        #     inun = xr.open_dataset(join(outdir, fn_inun), chunks=chunks).drop(['time']).squeeze()['inun']
        for exp in exposure:
            region_df[exp] = 0
        for i, (w, s, e, n) in enumerate(tiles):
            print('tile {:d}/{:d}: w:{:.1f} s:{:.1f}, e:{:.1f}, n:{:.1f}'.format(i+1, len(tiles), w, s, e, n))
            # tile inun
            inun_tile = inun.sel(lat=slice(s, n), lon=slice(w, e))
            # tile & stack mask
            mask_tile = mask_da.sel(lat=slice(s, n), lon=slice(w, e))
            mask_tile = mask_tile.stack(latlon=['lat', 'lon'])
            for exp in exposure:
                print(exp)
                # tile and overlay exposure
                exp_tile = exposure[exp].sel(lat=slice(s, n), lon=slice(w, e))
                exp_tile = xr.where(inun_tile > 0, exp_tile, 0)
                # stack exposure
                exp_tile = exp_tile.stack(latlon=['lat', 'lon'])
                # sum per region
                exp_region = exp_tile.groupby(mask_tile).sum('latlon')
                exp_region.name = exp
                exp_region_df = exp_region.to_dataframe()
                # add to total
                region_df.loc[exp_region_df.index, exp] += exp_region_df[exp].values
        region_df.to_csv(join(outdir, fn_inun.replace('.nc', '_exp.csv')))

