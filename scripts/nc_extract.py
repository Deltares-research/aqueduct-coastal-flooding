#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author: Dirk Eilander (contact: dirk.eilancer@vu.nl)
# Created: June 2nd

import os, sys, glob
from os.path import basename, dirname, join
import rasterio
from rasterio import features
from rasterio.transform import Affine
import numpy as np
import pandas as pd
import math
import xarray as xr
import geopandas as gp
import fiona
import click
import datetime
# multiprocessing libraries
from dask.diagnostics import ProgressBar
from multiprocessing.pool import ThreadPool
import dask
import dask.threaded

@click.command()
@click.argument('in_fns', nargs=-1, type=click.Path(readable=True, resolve_path=True))
@click.argument('out_fn', type=click.Path(writable=True, resolve_path=True))
@click.option('-e', '--extract_fn', required=True, type=click.Path(readable=True, resolve_path=True),
    help='CSV/SHP/GTiff filename with location information where to extract time series.')
@click.option('-v', '--var_name', multiple=True, required=True,
    help='NetCDF variable(s) name.')
@click.option('-d', '--dim', default='index',
    help='NetCDF output point/region dimension name [default: index]')
@click.option('-x', '--xdim', default='lon',
    help='NetCDF x-dimension name. [default: lon]')
@click.option('-y', '--ydim', default='lat',
    help='NetCDF y-dimension name. [default: lat]')
@click.option('-t', '--tdim', default='time',
    help='NetCDF time-dimension name. [default: time]')
@click.option('--csv_xdim', default='lon',
    help='CSV x column name; only used when <extract_fn> is a csv file. [default: lon]')
@click.option('--csv_ydim', default='lat',
    help='CSV y column name; only used when <extract_fn> is a csv file. [default: lat]')
@click.option('--reducer', type=click.Choice(['mean', 'sum']), default='sum',
    help='Function used to reduce region values; only used in combination with a polygon shape pr tif- file. [default: sum]')
@click.option('--keep_attrs', is_flag=True,
    help='Flag to keep <extract_fn> CSV/SHP metadata.')
@click.option('--num_workers', default=16,
    help='Number of threads used when reading nc data. [default: 16]')
@click.option('--chunksize', default=100,
    help='NetCDF time dimension chunksize. [default: 100]')
def main(in_fns, out_fn, extract_fn, var_name, dim='id',
         xdim='lon', ydim='lat', tdim='time', # netcdf dimension names
         csv_xdim='snap_lon', csv_ydim='snap_lat', keep_attrs=False,
         reducer='mean', num_workers=16, chunksize=100, **kwargs):
    nc_extract(in_fns, out_fn, extract_fn, var_name, dim=dim,
         xdim=xdim, ydim=ydim, tdim=tdim, # netcdf dimension names
         csv_xdim=csv_xdim, csv_ydim=csv_ydim, keep_attrs=keep_attrs,
         reducer=reducer, num_workers=num_workers, chunksize=chunksize, **kwargs)


def nc_extract(in_fns, out_fn, extract_fn, var_name, dim='id',
         xdim='lon', ydim='lat', tdim='time', # netcdf dimension names
         csv_xdim='snap_lon', csv_ydim='snap_lat', keep_attrs=False,
         reducer='mean', num_workers=16, chunksize=100, **kwargs):
    """Returns a netcdf with 1D time series data extracted from gridded (2D)
    netcdf files. The 1D time series can be extracted at points or regions.
    To extract time series at point locations provide a point shapefile or csv.
    To extract time series from regions provide a polygon shapefile.
    """
    # TODO epsg conversion {'init': 'epsg:4326'}
    # csv point file
    # resolve filenames with * in windows
    if (len(in_fns)==1) and (not os.path.isfile(in_fns[0])):
        fn_tmp = in_fns[0]
        in_fns = glob.glob(fn_tmp)
        if len(in_fns) == 0:
            raise IOError('Invalid input filename {}'.format(fn_tmp))

    points, mask = False, None
    if extract_fn.endswith('.csv') or extract_fn.endswith('.txt'):
        # click.echo('reading csv data from "{}"'.format(extract_fn))
        mapping = pd.read_csv(extract_fn) #, encoding='utf-8')
        if not csv_xdim in mapping.columns:
            raise ValueError('csv_xdim not found in csv header')
        if not csv_ydim in mapping.columns:
            raise ValueError('csv_ydim not found in csv header')
        points = True
    elif extract_fn.endswith('.tif'):
        # create pandas based on tif mask
        with rasterio.open(extract_fn, 'r') as src:
            mask = np.ma.masked_equal(src.read(1), src.nodata)
            mapping = pd.DataFrame(index=np.unique(mask.data[~mask.mask]))
    else:
        # click.echo('reading shp data from "{}"'.format(extract_fn))
        if '.zip' in extract_fn:
            kwargs.update(vfs = 'zip://'+extract_fn)
            extract_fn = ''
        try:
            mapping = gp.read_file(extract_fn, **kwargs)
        except IOError as e:
            raise IOError('no valid extract_fn, should be shp or csv file: {}'.format(str(e)))
        # convert point shp to pandas dataframe with x and y col
        if np.all([g.type == 'Point' for g in mapping.geometry]):
            points = True
            x, y = zip(*[p.coords[:][0] for p in mapping.geometry])
            mapping[csv_xdim], mapping[csv_ydim] = x, y
            mapping = mapping.drop(['geometry'], axis=1)

    # set index which is used as nc dimension
    if dim in mapping.columns:
        vals = mapping[dim].values
        if vals.size != np.unique(vals).size:
            raise ValueError('Output dimension {} should contain only unique values'.format(dim))
        mapping = mapping.set_index(dim)
    else:
        mapping.index = np.arange(len(mapping), dtype=int) + 1 # zero is background value
    mapping.index.name = dim

    # set output info
    today = datetime.datetime.today().strftime('%Y-%m-%d')
    history = "Created on {} with {:s}".format(today, os.path.basename(__file__))
    global_attrs = {
        'source_files': '; '.join([basename(fn) for fn in in_fns]),
        'mapping_file': basename(extract_fn),
        'history': history
        }

    if points:
        nc_points_ts(in_fns, out_fn, mapping, var_name=var_name,
           nc_xdim=xdim, nc_ydim=ydim, nc_tdim=tdim, chunksize=chunksize,
           csv_xdim=csv_xdim, csv_ydim=csv_ydim,  keep_attrs=keep_attrs,
           global_attrs=global_attrs, num_workers=num_workers)
    else:
        nc_region_ts(in_fns, out_fn, mapping, var_name=var_name, mask=mask,
            nc_xdim=xdim, nc_ydim=ydim, nc_tdim=tdim, chunksize=chunksize,
            keep_attrs=keep_attrs, global_attrs=global_attrs,
            num_workers=num_workers)

def nc_points_ts(in_nc_fns, out_nc_fn, mapping, var_name=['outflw'],
                nc_xdim='lon', nc_ydim='lat', nc_tdim='time', chunksize=100,
                csv_xdim='snap_lon', csv_ydim='snap_lat',
                keep_attrs=False, global_attrs={},  encoding={}, num_workers=16):
    # make sure var_name is a list
    if isinstance(var_name, str):
        var_name = [var_name]
    else:
        var_name = [name for name in var_name]
    out_dim = mapping.index.name
    with dask.set_options(get=dask.threaded.get, num_workers=num_workers): #, ProgressBar():
        # read directory with netcdf output files
        # use either a string glob in the form “path/to/my/files/*.nc” or an explicit list of files to open
        # click.echo('read nc data from: \n "{}"'.format('",\n "'.join(in_nc_fns)))
        chunks = {nc_tdim: chunksize}
        ds = xr.open_mfdataset(in_nc_fns, chunks=chunks)
        # check if all var_names in dataset
        check_var_names = np.array([name in ds.variables.keys() for name in var_name])
        if not np.all(check_var_names):
            missing = np.array(var_name)[check_var_names==False].tolist()
            raise ValueError('Variable(s) {} not in netcdf dataset'.format(', '.join(missing)))
        # read ts at lat lon coordinates (snapped to the cell centers of the grid)
        x, y = mapping[csv_xdim].values, mapping[csv_ydim].values
        slice_args = {nc_xdim: x, nc_ydim: y, 'dim': out_dim}
        ds_pnt = ds[var_name].sel_points(method='nearest', **slice_args)
        # add new dimension variable
        ds_pnt[out_dim] = xr.Variable([out_dim], mapping.index.values)
        # merge csv metadata if keep_attrs
        if keep_attrs:
            rm_dict = {nc_xdim: nc_xdim+'_nc', nc_ydim: nc_ydim+'_nc'}
            ds_pnt = ds_pnt.rename(rm_dict)
            ds_out = xr.merge([mapping.to_xarray(), ds_pnt])
        else:
            ds_out = ds_pnt
        # save to netcdf (here's were all the work takes place)
        # click.echo('save time series data to "{}"'.format(out_nc_fn))
        ds_out.attrs = global_attrs
        encoding = ds.encoding
        encoding.update({name: {'zlib': True} for name in var_name})
        ds_out.sortby(nc_tdim).to_netcdf(out_nc_fn, encoding=encoding)
        # close files and finish
        ds.close()
        click.echo('finished extracting time series point data')
    return None

def nc_region_ts(in_nc_fns, out_nc_fn, region_gdf, var_name=['fldfrc'],
                 mask=None, reducer='mean', nc_xdim='lon', nc_ydim='lat', nc_tdim='time',
                 chunksize=50, keep_attrs=False, global_attrs={}, num_workers=16):
    reducers = ['mean', 'sum']
    if reducer not in reducers:
        msg = 'reducer not implemeted, select from {}'.format(', '.join(reducers))
        raise NotImplementedError(msg)
    # make sure var_name is a list
    if isinstance(var_name, str):
        var_name = [var_name]
    else:
        var_name = [name for name in var_name]

    with dask.set_options(get=dask.threaded.get, num_workers=num_workers): #, ProgressBar():
        # click.echo('read nc data from: \n "{}"'.format('",\n "'.join(in_nc_fns)))
        chunks = {nc_tdim: chunksize}
        ds = xr.open_mfdataset(in_nc_fns, chunks=chunks)
        # check if all var_names in dataset
        check_var_names = np.array([name in ds.variables.keys() for name in var_name])
        if not np.all(check_var_names):
            missing = np.array(var_name)[check_var_names==False].tolist()
            raise ValueError('Variable(s) {} not in netcdf dataset'.format(', '.join(missing)))
        # create mask
        out_dim = region_gdf.index.name
        if mask is None:
            mask = nc_rasterize(ds, region_gdf, col_name=out_dim, y_name=nc_ydim, x_name=nc_xdim)
        else:
            if not (ds[var_name[0]].isel(time=1).shape == mask.shape):
                raise IndexError('Mask and dataset x, and y dimensions do not agree')
        mask_da = xr.DataArray(mask, coords=[ds[nc_ydim], ds[nc_xdim]], dims=(nc_ydim, nc_xdim))
        mask_da.name = out_dim
        mask_da.attrs.update(long_name='region mask for {}'.format(out_dim))
        # regional mean
        xvar = ds[var_name].stack(latlon=[nc_ydim, nc_xdim])
        ds_region_grps = xvar.groupby(mask_da.stack(latlon=[nc_ydim, nc_xdim]))
        if reducer == 'mean':
            ds_region = ds_region_grps.mean('latlon')
        elif reducer == 'sum':
            ds_region = ds_region_grps.sum('latlon')
        # add index
        ds_region[out_dim] = xr.Variable([out_dim], region_gdf.index.values)
        # merge with shp metadata if keep_attrs
        if keep_attrs:
            # convert shp metadata to xarray
            ds_gdf = region_gdf.drop(['geometry'], axis=1).to_xarray()
            ds_out = xr.merge([ds_gdf, ds_region])
        else:
            ds_out = ds_region
        # add mask to output
        ds_out = ds_out.assign(mask=mask_da)
        # save to file
        # click.echo('save time series data to "{}"'.format(out_nc_fn))
        encoding=ds.encoding
        encoding.update({name: {'zlib': True} for name in var_name})
        ds_out.attrs = global_attrs
        ds_out.sortby(nc_tdim).to_netcdf(out_nc_fn, encoding=encoding)
        # close files and finish
        ds.close()
        # click.echo('finished extracting time series region data')
    return None

# utils
def latlon2transform(lat, lon, cell_center=True):
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    resX = (lon[-1] - lon[0]) / (len(lon) - 1)
    resY = (lat[-1] - lat[0]) / (len(lat) - 1)
    trans = Affine.translation(lon[0] - resX/2.*cell_center, lat[0] - resY/2.*cell_center)
    scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
    return trans * scale

def rasterize(gdf, lat=None, lon=None, fill=np.nan, col_name='index',
                cell_center=True, **kwargs):
    """Rasterize the index value geopandas dataframe onto a raster defined
    by either its lat/lon coordinates (1d arrays) or
    transform (rasterio tranfrom) and out_size (tuple).
    """
    out_size = (len(lat), len(lon))
    transform = latlon2transform(lat, lon, cell_center=cell_center)
    shapes = [(geom, value) for geom, value in zip(gdf.geometry, gdf.reset_index()[col_name])]
    raster = features.rasterize(shapes, out_shape=out_size,
                                fill=fill, transform=transform, **kwargs)
    return raster

def nc_rasterize(ds, gdf, col_name='index', x_name='lon', y_name='lat', nodata=0):
    # rasterize regions
    lat, lon = ds[y_name].data, ds[x_name].data
    mask = rasterize(gdf.reindex(), lat=lat, lon=lon, fill=nodata, dtype=int, col_name=col_name)
    return np.ma.masked_equal(mask, nodata)

if  __name__ == "__main__":
    main()
