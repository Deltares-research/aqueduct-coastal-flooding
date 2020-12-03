
# coding: utf-8

# In[1]:

import numpy as np
import rasterio
import rasterio.mask
from rasterio import features
from rasterio.transform import Affine, from_bounds
from shapely.geometry import box, mapping
import os
from os.path import join, basename, dirname
import geopandas as gpd


# In[2]:

# tiling
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

def rasterize(gdf, col_name='index', # rasterize
              transform=None, out_size=None, # option 1
              lat=None, lon=None, cell_center=True, # option 2
              fill=np.nan, dtype=None, **kwargs):
    """Rasterize the index value geopandas dataframe onto a raster defined
    by either its lat/lon coordinates (1d arrays) or
    transform (rasterio tranfrom) and out_size (tuple).
    """

    transform_option = (transform is not None) and (out_size is not None)
    latlon_option = (lat is not None) & (lon is not None)
    if (not transform_option) and (not latlon_option):
        raise ValueError("incorrect input provide either transform and out_size or lat & lon 1d arrays")
    elif latlon_option:
        out_size = (len(lat), len(lon))
        transform = latlon2transform(lat, lon, cell_center=cell_center)
    shapes = [(geom, value) for geom, value in zip(gdf.geometry, gdf.reset_index()[col_name])]
    raster = features.rasterize(shapes, out_shape=out_size,
                                fill=fill, transform=transform,
                                **kwargs)
    if dtype is not None:
        raster = np.array(raster).astype(dtype)
    return raster


# In[3]:

outdir = r'/p/11201999-eucp/42b_coastal/Coastal_Cities_EU/input'
fn_wmask = r'/p/11201999-eucp/42b_coastal/Coastal_Cities_EU/input/WATER_europe_v3.tif'
fn_ocean_shp = r'/p/1209884-aqueduct/Datasets/coastline/ne_10m_ocean.shp'
fn_dem = r'/p/11201999-eucp/42b_coastal/Coastal_Cities_EU/input/MERIT_europe_v2.tif'
#open files
fn_out = r'/p/11201999-eucp/42b_coastal/Coastal_Cities_EU/input/MERIT_europe_v2_masked.tif'


# In[5]:



gdf = gpd.read_file(fn_ocean_shp)
spatial_index = gdf.sindex


# In[11]:

with rasterio.open(fn_dem) as src, rasterio.open(fn_wmask) as src_wmask:
    prof = src.profile
    bounds = src.bounds
    res = src.res
#     tile_size = prof['blockxsize']*res[0], prof['blockysize']*res[1]
    tile_size = 4, 4
    tiles = get_tiles(bounds, *tile_size)
    chunks={'lon':int(tile_size[0]/res[0]), 'lat': int(tile_size[1]/res[1])}
    prof['blockxsize'] = chunks['lon']
    prof['blockysize'] = chunks['lat']
    prof['bigtiff'] = 'YES'
    prof['compress'] = 'lzw'
    
    with rasterio.open(fn_out, 'w', **prof) as dst:
        for i, bounds in enumerate(tiles):
            w, s, e, n = bounds
            print('tile {:d}/{:d}: w:{:.1f} s:{:.1f}, e:{:.1f}, n:{:.1f}'.format(i+1, len(tiles), w, s, e, n))
            wdw = rasterio.windows.from_bounds(w, s, e, n, dst.transform)
            bbox = box(*bounds)
            ft = mapping(bbox)

            # read dem tile from gtiff
            dem_tile, transform = rasterio.mask.mask(src, [ft], crop=True)
            dem_tile = dem_tile[0, :, :]
            if np.all(dem_tile == src.nodata):
                # ocean only tile
                dst.write(dem_tile,  window=wdw, indexes=1)
                continue
                            
            # read ocean shape for tile
            sidx = list(spatial_index.intersection(bounds))
            clipped = None
            if len(sidx) > 0:
                 # Clip the data - with these data
                clipped = gdf.iloc[sidx].copy()
                clipped["geometry"] = clipped.intersection(bbox)
                clipped = clipped[clipped.geometry.notnull()]
            if clipped is not None and len(clipped) > 0:
                # rasterize ocean
                clipped['value'] = 1
                ocean_tile = rasterize(clipped, 'value', transform=transform, out_size=dem_tile.shape)
            else:
                # land only tile
                dst.write(dem_tile,  window=wdw, indexes=1)
                continue
                
            # read wmask tile from gtiff
            wmask_tile, _ = rasterio.mask.mask(src_wmask, [ft], crop=True) 
            wmask_tile = wmask_tile[0, :, :]
            wmask_tile[wmask_tile==src_wmask.nodata] = 1
                 
            # fix mask
            dem_tile[np.logical_and(ocean_tile==1, wmask_tile==1)] = src.nodata
            
            # write output
            dst.write(dem_tile,  window=wdw, indexes=1)

