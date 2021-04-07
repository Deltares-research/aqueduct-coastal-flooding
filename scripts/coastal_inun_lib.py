# -*- coding: utf-8 -*-
"""
Created on Wed Jul 08 16:12:29 2015

@author: winsemi

$Id: coastal_inun_lib.py 445 2018-03-28 14:39:49Z eilan_dk $
$Date: 2018-03-28 16:39:49 +0200 (Wed, 28 Mar 2018) $
$Author: eilan_dk $
$Revision: 445 $
$HeadURL: https://repos.deltares.nl/repos/aqueduct/trunk/python/coastal_inundation/inun_clean/scripts/coastal_inun_lib.py $
$Keywords: $

"""

import sys
import os
import shutil
from scipy import ndimage as nd

# import admin packages
import configparser
import logging
import logging.handlers
# from optparse import OptionParser

# import general packages
import numpy as np
import math
# import pyproj

from osgeo import osr, gdal, gdalconst
from gdal_warp import gdal_warp
# from osgeo import ogr
# from hydrotools import gis
# import specific packages
import netcdf_funcs as nc_funcs
import netCDF4 as nc

def setlogger(logfilename, logReference, verbose=True):
    """
    Set-up the logging system. Exit if this fails
    """
    try:
        #create logger
        logger = logging.getLogger(logReference)
        logger.setLevel(logging.DEBUG)
        ch = logging.handlers.RotatingFileHandler(logfilename, maxBytes=10*1024*1024, backupCount=5)
        console = logging.StreamHandler()
        console.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        #create formatter
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        #add formatter to ch
        ch.setFormatter(formatter)
        console.setFormatter(formatter)
        #add ch to logger
        logger.addHandler(ch)
        logger.addHandler(console)
        logger.debug("File logging to " + logfilename)
        return logger, ch
    except IOError:
        print("ERROR: Failed to initialize logger with logfile: " + logfilename)
        sys.exit(1)


def closeLogger(logger, ch):
    logger.removeHandler(ch)
    ch.flush()
    ch.close()
    return logger, ch


def close_with_error(logger, ch, msg):
    logger.error(msg)
    logger, ch = closeLogger(logger, ch)
    del logger, ch
    sys.exit(1)


def open_conf(fn):
    config = configparser.SafeConfigParser()
    config.optionxform = str

    if os.path.exists(fn):
        config.read(fn)
    else:
        print("Cannot open config file: " + fn)
        sys.exit(1)

    return config


def configget(config, section, var, default, datatype='str'):
    """
    Gets a string from a config file (.ini) and returns a default value if
    the key is not found. If the key is not found it also sets the value
    with the default in the config-file

    Input:
        - config - python configparser object
        - section - section in the file
        - var - variable (key) to get
        - default - default value
        - datatype='str' - can be set to 'boolean', 'int', 'float' or 'str'

    Returns:
        - value (str, boolean, float or int) - either the value from the config file or the default value
    """

    try:
        if datatype == 'int':
            ret = config.getint(section, var)
        elif datatype == 'float':
            ret = config.getfloat(section, var)
        elif datatype == 'boolean':
            ret = config.getboolean(section, var)
        else:
            ret = config.get(section, var)
    except:
        ret = default

    return ret


def get_gdal_extent(filename):
    ''' Return list of corner coordinates from a dataset'''
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    if ds is None:
        print("invallid file name or file ")
    else:
        gt = ds.GetGeoTransform()
        # 'top left x', 'w-e pixel resolution', '0', 'top left y', '0', 'n-s pixel resolution (negative value)'
        nx, ny = ds.RasterXSize, ds.RasterYSize
        xmin = np.float64(gt[0])
        ymin = np.float64(gt[3]) + np.float64(ny) * np.float64(gt[5])
        xmax = np.float64(gt[0]) + np.float64(nx) * np.float64(gt[1])
        ymax = np.float64(gt[3])
        ds = None  # close dataset
    return xmin, ymin, xmax, ymax


def get_gdal_geotransform(filename):
    ''' Return geotransform of dataset'''
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    if ds is None:
        logging.warning('Could not open {:s} Shutting down'.format(filename))
        sys.exit(1)
    # Retrieve geoTransform info
    gt = ds.GetGeoTransform()
    ds = None  # close dataset
    return gt


def get_gdal_axes(filename, logging=logging):
    # TODO: check logger not used ?
    geotrans = get_gdal_geotransform(filename)
    # Retrieve geoTransform info
    originX = geotrans[0]
    originY = geotrans[3]
    resX = geotrans[1]
    resY = geotrans[5]

    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    if ds is None:
        logging.warning('Could not open {:s} Shutting down'.format(filename))
        sys.exit(1)
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX+resX/2, originX+resX/2+resX*(cols-1), cols)
    y = np.linspace(originY+resY/2, originY+resY/2+resY*(rows-1), rows)
    ds = None  # close dataset
    return x, y


def get_gdal_fill(filename, logging=logging):
    ds = gdal.Open(filename, gdal.GA_ReadOnly)
    if ds is None:
        logging.warning('Could not open {:s} Shutting down'.format(filename))
        sys.exit(1)
    # Retrieve geoTransform info
    RasterBand = ds.GetRasterBand(1)
    fill_val = RasterBand.GetNoDataValue()
    ds = None  # close dataset
    return fill_val


def get_gdal_rasterband(filename, band=1, logging=logging):
    """

    :param filename: GDAL compatible raster file to read from
    :param band: band number (default=1)
    :param logging: logging object
    :return: gdal dataset object, gdal rasterband object
    """
    ds = gdal.Open(filename)
    if ds is None:
        logging.warning('Could not open {:s} Shutting down'.format(filename))
        sys.exit(1)
    # Retrieve geoTransform info
    return ds, ds.GetRasterBand(band)   # there's only 1 band, starting from 1


def get_gt(x, y):
    # make the geotransform
    xul = x[0]-(x[1]-x[0])/2
    xres = x[1]-x[0]
    yul = y[0]+(y[0]-y[1])/2
    yres = y[1]-y[0]
    geotrans = [xul, xres, 0, yul, 0, yres]
    return geotrans


def MapToPixel(mx, my, gt):
    """Convert map to pixel coordinates
    :param  mx:    Input map x coordinate (double)
    :param  my:    Input map y coordinate (double)
    :param  gt:    Input geotransform (six doubles)
    :return: px,py Output coordinates (two ints)
    """
    if gt[2] + gt[4] == 0:  # Simple calc, no inversion required
        px = (mx - gt[0]) / gt[1]
        py = (my - gt[3]) / gt[5]
    else:
        px, py = ApplyGeoTransform(mx, my, gdal.InvGeoTransform(gt))
    return int(px), int(py)


def ApplyGeoTransform(inx, iny, gt):
    """Apply a geotransform
    :param  inx:       Input x coordinate (double)
    :param  iny:       Input y coordinate (double)
    :param  gt:        Input geotransform (six doubles)

    :return: outx,outy Output coordinates (two doubles)
    """
    outx = gt[0] + inx * gt[1] + iny * gt[2]
    outy = gt[3] + inx * gt[4] + iny * gt[5]
    return (outx, outy)


def gdal_sample_points(lat, lon, raster_file, win_size=1, func=np.mean):
    """
    sample data from a raster based on lat lon location

    :param lat: list or array with latitude coordinates
    :param lon: list or array with longitude coordinates
    :param raster_file: filename of raster file
    :param win_size: if spatial aggregate from window: give window size (default 1)
    :param func: aggregation function
    :return:
    """
    src_ds = gdal.Open(raster_file)
    if src_ds is None:
        logging.warning('Could not open {:s} Shutting down'.format(raster_file))
        sys.exit(1)
    else:
        gt = src_ds.GetGeoTransform()
        rb = src_ds.GetRasterBand(1)
        cols = src_ds.RasterXSize
        rows = src_ds.RasterYSize
        values = np.array([])
        if isinstance(lat, np.ndarray):
            lat = lat.tolist()
            lon = lon.tolist()
        if (np.size(lat) == 1) & (not isinstance(lat, list)):
            lon = [lon]
            lat = [lat]
        for my, mx in zip(lat, lon):

            # Convert from map to pixel coordinates.
            px, py = MapToPixel(mx, my, gt)
            if (px >= 0) & (px <= cols) & (py >= 0) & (py <= rows):  # check if within map extent
                if win_size > 1:
                    # limit window to raster extent
                    xoff = max(int(px - win_size / 2.), 0)
                    yoff = max(int(py - win_size / 2.), 0)
                    win_xsize = min(cols - xoff, win_size) + min(int(px - win_size / 2.), 0)
                    win_ysize = min(rows - yoff, win_size) + min(int(py - win_size / 2.), 0)
                    intval = rb.ReadAsArray(xoff, yoff, win_xsize, win_ysize)
                else:
                    intval = rb.ReadAsArray(int(px), int(py), win_size, win_size)
                if intval is None:
                    intval = [np.nan]
                if func is not None:
                    values = np.append(values, func(intval))
                else:
                    values = np.append(values, intval)
            else:
                values = np.append(values, np.nan)
        src_ds = None  # close file
        return values


def discretize_raster_bounds(fn_raster, x_tile, y_tile, x_overlap, y_overlap,
                             fn_bounds=None, x_var=None, y_var=None):
    """discretize raster into tiles with size (x_tile + 2*x_overlap, y_tile + 2*y_overlap)
    find boundary conditions within tile and check if raster has data values in tile"""
    tile_list = []
    n = 0

    # read raster
    x, y = get_gdal_axes(fn_raster)
    gt_dem = get_gdal_geotransform(fn_raster)

    # read boundary conditions
    if fn_bounds is not None:
        x_coords, y_coords = nc_funcs.get_netcdf_axes(fn_bounds, x_var=x_var, y_var=y_var)

    for x_loop in range(0, len(x), x_tile):
        start_col = np.maximum(x_loop, 0)  # cols
        end_col = np.minimum(x_loop + x_tile, len(x))
        for y_loop in range(0, len(y), y_tile):
            n += 1
            start_row = np.maximum(y_loop, 0)  # rows
            end_row = np.minimum(y_loop + y_tile, len(y))

            # get tile extent and overlap
            tile = tile_admin(x, y, gt_dem, start_col, end_col, start_row, end_row, x_overlap, y_overlap)
            tile['i'] = n
            tile_list.append(tile)

            # get indices boundary conditions per tile
            if fn_bounds is not None:
                xmin, ymin, xmax, ymax = tile['extent']
                bounds_select_bool = np.all(np.array([x_coords >= xmin, x_coords <= xmax,
                                                      y_coords >= ymin, y_coords <= ymax]), axis=0).tolist()
                tile['bounds_select'] = np.where(bounds_select_bool)[0].astype(int).tolist()
                tile['has_bounds'] = int(np.sum(bounds_select_bool) > 0)

            # check if DEM data in tile (i.e.: not ocean)
            ds, raster = get_gdal_rasterband(fn_raster)
            cols, rows = tile['size']
            origin_col, origin_row = tile['origin']
            cells_with_data = np.sum(raster.ReadAsArray(int(origin_col), int(origin_row), int(cols), int(rows)) != raster.GetNoDataValue())
            tile['has_DEM_data'] = int(np.nansum(cells_with_data) > 0)
            ds = None  # close dataset

    return tile_list


def discretize_raster(fn_raster, x_tile, y_tile, x_overlap, y_overlap):
    """discretize raster into tiles with size (x_tile + 2*x_overlap, y_tile + 2*y_overlap)
    find boundary conditions within tile and check if raster has data values in tile"""
    tile_list = []
    n = 0

    # read raster
    x, y = get_gdal_axes(fn_raster)

    # read boundary conditions
    for x_loop in range(0, len(x), x_tile):
        start_col = np.maximum(x_loop, 0)  # cols
        end_col = np.minimum(x_loop + x_tile, len(x))
        for y_loop in range(0, len(y), y_tile):
            n += 1
            start_row = np.maximum(y_loop, 0)  # rows
            end_row = np.minimum(y_loop + y_tile, len(y))

            # get tile extent and overlap
            tile = tile_admin(x, y, gt_dem, start_col, end_col, start_row, end_row, x_overlap, y_overlap)
            tile['i'] = n
            tile_list.append(tile)

    return tile_list

def tile_admin(x, y, gt, start_col, end_col, start_row, end_row, overlap_col, overlap_row):
    """
    function to read tile from a raster (GDAL object) based on rows and columns

    """
    # output dictionary
    tile = {'start_col': int(start_col), 'end_col': int(end_col), 'start_row': int(start_row), 'end_row': int(end_row)}

    # calc extent
    col_overlap_min = start_col - np.max([start_col - overlap_col, 0])
    col_overlap_max = np.min([end_col + overlap_col, len(x)]) - end_col
    row_overlap_min = start_row - np.max([start_row - overlap_row, 0])
    row_overlap_max = np.min([end_row + overlap_row, len(y)]) - end_row
    origin_col = start_col - col_overlap_min
    origin_row = start_row - row_overlap_min
    cols = (end_col + col_overlap_max) - (start_col - col_overlap_min)
    rows = (end_row + row_overlap_max) - (start_row - row_overlap_min)
    tile['overlap'] = [int(col_overlap_min), int(col_overlap_max), int(row_overlap_min), int(row_overlap_max)]

    # x & y coordinates terrain
    x_tile = x[np.arange(origin_col, origin_col+cols)]
    y_tile = y[np.arange(origin_row, origin_row+rows)]
    tile['extent'] = np.min(x_tile), np.min(y_tile), np.max(x_tile), np.max(y_tile)
    tile['size'] = (int(cols), int(rows))
    tile['origin'] = (int(origin_col), int(origin_row))
    tile['gt'] = (x[origin_col], gt[1], gt[2],
                  y[origin_row], gt[4], gt[5])
    tile['x'] = x_tile.tolist()
    tile['y'] = y_tile.tolist()

    return tile


def read_raster_tile(dem_file, tile_dict, NoDataValue = -9999):
    """
    function to read tile from a raster (GDAL object) based on rows and columns

    """
    # read tile properties
    cols, rows = tile_dict['size']
    origin_col, origin_row = tile_dict['origin']

    # read terrain
    # Open terrain data for reading
    ds, raster = get_gdal_rasterband(dem_file)
    tile = raster.ReadAsArray(origin_col, origin_row, cols, rows)
    cells_with_data = tile != raster.GetNoDataValue()

    # check if elevation data in tile
    if np.nansum(cells_with_data) > 0:
        # change NoDataValue to -9999 (see inun.pcr_preprocess)
        tile[tile == raster.GetNoDataValue()] = NoDataValue
        tile = np.ma.masked_equal(tile, NoDataValue)

        # create new gdal object in memory
        gdal_type = gdalconst.GDT_Float32
        ds_tile = gdal.GetDriverByName('MEM').Create('', cols, rows, 1, gdal_type)
        ds_tile.SetGeoTransform(tile_dict['gt'])
        ds_tile.SetProjection(ds.GetProjection())
        ds_tile.GetRasterBand(1).WriteArray(tile.data, 0, 0)
        ds_tile.GetRasterBand(1).SetNoDataValue(NoDataValue)
        ds = None  # close dataset
        return tile, ds_tile

    else:
        ds = None  # close dataset
        return None, None

def hydroldd2pcrldd(ldd):
    lddout = np.ma.masked_all_like(ldd)
    lddout[ldd == 8] = 1
    lddout[ldd == 4] = 2
    lddout[ldd == 2] = 3
    lddout[ldd == 16] = 4
    lddout[ldd == 0] = 5  # pits
    lddout[ldd == 1] = 6
    lddout[ldd == 32] = 7
    lddout[ldd == 64] = 8
    lddout[ldd == 128] = 9
    return lddout

def distance_on_unit_sphere(lat, res, R=6373000):
    """
    the local distance varies with the latitude of the earth. this function
    calculates the local resolution of a grid given its latitude

    :param lat: latitude coordinate (scalar
    :param res: local resolution in degrees
    :param R:   radius earth in local unit (R = 6373 for km)
    :return:    resolution in local unit
    """
    d = []

    # define two points along cell diagonal with res distance
    r = res/2. * math.sqrt(0.5)
    lon1, lon2 = -r, r

    for l in lat:
        lat1, lat2 = l-r, l+r

        # Convert latitude and longitude to
        # spherical coordinates in radians.
        degrees_to_radians = math.pi/180.0

        # phi = 90 - latitude
        phi1 = (90.0 - lat1)*degrees_to_radians
        phi2 = (90.0 - lat2)*degrees_to_radians

        # theta = longitude (
        theta1 = lon1*degrees_to_radians
        theta2 = lon2*degrees_to_radians

        # Compute spherical distance from spherical coordinates.
        # For two locations in spherical coordinates
        # (1, theta, phi) and (1, theta', phi')
        # cosine( arc length ) =
        #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
        # distance = rho * arc length

        cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) +
               math.cos(phi1)*math.cos(phi2))
        d.append(math.acos(cos) * R)

    # multiply arc by the radius of the earth in local units to get length.
    return d


def makeDir(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)

# def cleanDir(tempdir):
#     # Clean up the temporary folder after ourselves.
#     for filef in os.listdir(tempdir):
#         os.unlink(os.path.join(tempdir, filef))
#     os.rmdir(tempdir)

def cleanDir(tempdir, logger):
    try:
        # remove all
        shutil.rmtree(tempdir)
    except Exception as e:
        # remove per file
        try:
            for fn in os.listdir(tempdir):
                del_path = os.path.join(tempdir, fn)
                os.unlink(del_path)
            del_path = tempdir
            os.rmdir(del_path)
        except Exception as e:
            logger.error('{:s} cannot be deleted'.format(del_path))

def planar_flood(terrain, flood_level):
    return np.maximum(flood_level - terrain, 0)


def check_input_fn(fn, logger):
    if fn is None:
        return fn
    elif len(fn) == 0:
        fn = None
    else:
        if not os.path.exists(fn):
            logger.error('path to file {:s} cannot be found'.format(fn))
            sys.exit(1)
    return fn


def fill(data, invalid=None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid')
    by the value of the nearest valid data cell

    Input:
        data:    numpy array of any dimension
        invalid: a binary array of same shape as 'data'. True cells set where data
                 value should be replaced.
                 If None (default), use: invalid  = np.isnan(data)

    Output:
        Return a filled array.
    """
    if invalid is None:
        invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]


def gdal_read_slice(extent, raster_file, return_gdal=False, NNfill=False):
    """
    sample data from a raster based on lat lon location

    :param extent: list or array with min_lon, min_lat, max_lon, max_lat coordinates
    :param raster_file: filename of raster file
    :param return_gdal: returns gdal object instead of data
    :param NNfill: fills NoDataValues with NearestNeighbor
    :return: numpy array with slice of raster image or gdal object
    """
    src_ds = gdal.Open(raster_file)
    if src_ds is None:
        logging.warning('Could not open {:s} Shutting down'.format(raster_file))
        sys.exit(1)
    else:
        gt = src_ds.GetGeoTransform()
        rb = src_ds.GetRasterBand(1)
        fill_val = rb.GetNoDataValue()

        # Convert from map to pixel coordinates.
        min_lon, min_lat, max_lon, max_lat = extent
        if gt[5] > 0:
            px_left, py_bottom = MapToPixel(min_lon, max_lat, gt)
            px_right, py_top = MapToPixel(max_lon, min_lat, gt)
        else:
            px_left, py_top = MapToPixel(min_lon, max_lat, gt)
            px_right, py_bottom = MapToPixel(max_lon, min_lat, gt)

        cols = px_right-px_left+1
        rows = py_bottom-py_top+1

        # check extend in in raster
        if (py_top >= 0) & (px_left >= 0) & (py_bottom < src_ds.RasterYSize) & (px_right < src_ds.RasterXSize):
            data = rb.ReadAsArray(px_left, py_top, cols, rows)
            if gt[5] > 0:
                data = np.flipud(data)
                gt = (gt[0], gt[1], gt[2], gt[3], gt[4], -gt[5])
            data = np.ma.masked_equal(data, fill_val)

            if NNfill:
                data = fill(data, data.mask)
            if not return_gdal:
                src_ds = None
                return data

            else:
                # create gdal object
                src_mem = gdal.GetDriverByName('MEM').Create('', cols, rows, 1, gdalconst.GDT_Float32)
                new_geo = (min_lon, gt[1], gt[2],
                           max_lat, gt[4], gt[5])
                src_mem.SetGeoTransform(new_geo)
                src_mem.SetProjection(src_ds.GetProjection())
                src_mem.GetRasterBand(1).WriteArray(data, 0, 0)
                src_mem.GetRasterBand(1).SetNoDataValue(fill_val)
                src_ds = None  # close dataset
                return src_mem

        # extent outside grid
        else:
            src_ds = None  # close dataset
            return None


def read_warp_fill_grid_slice(fn, extent, ds_clone, fillNN=True, fill_constant=0, test=False):
    """

    :param fn:       filename of GeoTIFF raster to read
    :param extent:   list or array with min_lon, min_lat, max_lon, max_lat coordinates
    :param ds_clone: gdal clone object with terrain data
    :param fillNN:   if True gaps are filled based on nearest neighbor
    :param fill_constant constant value to use if fillNN == False or no subs data in tile
    :return:
    """
    # read slice from grid based on extent in lat lon
    ds_in = gdal_read_slice(extent, fn, return_gdal=True)

    # check if data in slice -> extent within raster
    if ds_in is not None:
        # read clone data
        clone_data = ds_clone.GetRasterBand(1).ReadAsArray(0, 0)
        fill_val = ds_clone.GetRasterBand(1).GetNoDataValue()
        clone_data = np.ma.masked_equal(clone_data, fill_val)

        # reproject if different slice
        if ds_in.GetRasterBand(1).ReadAsArray(0, 0).shape != clone_data.shape:
            ds_proj = gdal_warp('', ds_clone, '', format='MEM', ds_in=ds_in)
        else:
            ds_proj = ds_in

        # fill pixels with terrain data but no data in data
        data = ds_proj.GetRasterBand(1).ReadAsArray(0, 0)
        fill_val = ds_proj.GetRasterBand(1).GetNoDataValue()
        data[data == fill_val] = np.nan
        mask = (~clone_data.mask & np.isnan(data))

        if np.nansum(mask) > 0:
            if fillNN:
                data = fill(data)
                data[clone_data.mask] = np.nan
                mask = (~clone_data.mask & np.isnan(data))

            data[mask] = fill_constant
        ds_in = None
        ds_proj = None
        return data

    # slice extent outside raster
    else:
        return None
