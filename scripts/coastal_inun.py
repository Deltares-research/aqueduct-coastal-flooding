# -*- coding: utf-8 -*-
"""
Created on Wed Jul 08 16:12:29 2015

@author: winsemi

$Id: coastal_inun.py 528 2018-06-19 08:41:05Z eilan_dk $
$Date: 2018-06-19 10:41:05 +0200 (Tue, 19 Jun 2018) $
$Author: eilan_dk $
$Revision: 528 $
$HeadURL: https://repos.deltares.nl/repos/aqueduct/trunk/python/coastal_inundation/inun_clean/scripts/coastal_inun.py $
$Keywords: $

"""

import shutil
# import admin packages
from optparse import OptionParser
import datetime
import sys
import time
from multiprocessing import cpu_count, Pool
from shutil import copyfile
import os as os
import traceback
import uuid

# import general packages
import numpy as np
import netCDF4 as nc
import json
from osgeo import osr

# import specific packages
import coastal_inun_lib as cl
from inun import flood_tile
import netcdf_funcs as nc_funcs
# import matplotlib.pyplot as plt
from gdal_writemap import gdal_writemap
from gdal_readmap import gdal_readmap
import pdb

def process_tile(tile):
    """
    per tile with raster data and boundary conditions calculate flood extent
    first read data for scenario
    then process flood routine
    this part can be processed parallel on multiple cores

    :param tile:            list with tile admin [list]
    :param options:         options from ini file and command line
    :param boundaries:      boundary water levels
    :param x_boundaries:    x coordinates boundary conditions [1d array]
    :param y_boundaries:    y coordinates boundary conditions [1d array]
    :param output:          queue for output = flood 2d array

    """
    time0 = time.time()
    options = tile['options']
    print('processing tile {:03d}:'.format(tile['i']))
    print('xmin: {:d} xmax: {:d} ymin: {:d} ymax: {:d}; top-left corner (lon, lat): {:.2f}, {:.2f}'.format(
        tile['start_col'], tile['end_col'], tile['start_row'], tile['end_row'],
        tile['extent'][0], tile['extent'][3]))

    if tile['has_DEM_data']:
        # read tile DEM data
        terrain, ds_terrain = cl.read_raster_tile(options.dem_file, tile)

        if tile['has_bounds']:
            # read boundary conditions
            bounds_select = tile['bounds_select']
            with nc.Dataset(options.boundary, 'r') as a:
                # read the boundary conditions. This may also be a time series if a full event simulation is passed.
                boundaries = np.atleast_1d(np.squeeze(a.variables[options.boundary_variable][:, bounds_select]))
                x_boundaries = np.atleast_1d(np.squeeze(a.variables[options.x_var][:][bounds_select]))
                y_boundaries = np.atleast_1d(np.squeeze(a.variables[options.y_var][:][bounds_select]))
                assert boundaries.size == x_boundaries.size == y_boundaries.size != 0, "wrong input format"
        else:
            boundaries = []

        # if one or more boundaries with finite values
        if len(boundaries) > 0:
            if options.ldd_file is not None:
                # ldd needs exactly same raster definition as DEM for this to work!
                # here the map is read, if different in shape warped to clone and gaps where terrain are filled
                ldd_map = cl.read_warp_fill_grid_slice(options.ldd_file, tile['extent'], ds_terrain, fill_constant=np.nan)
            else:
                ldd_map = None  # ldd will be calculated on the fly (expensive)

            if options.subsidence_map is not None:
                # resamples and smooths (using Cubic interpolation) the grid to the dtm resolution
                subs2030_map = cl.read_warp_fill_grid_slice(options.subsidence_map, tile['extent'], ds_terrain,
                                                            fillNN=False, fill_constant=0)
                if subs2030_map is not None:
                    terrain = terrain - subs2030_map
                    del subs2030_map  # cleanup memory

            if options.water_perc_file is not None:
                # ldd needs exactly same raster definition as DEM for this to work!
                # here the map is read, if different in shape warped to clone and gaps where terrain are filled
                water_perc = cl.gdal_read_slice(tile['extent'], options.water_perc_file)
            else:
                water_perc = None  # ldd will be calculated on the fly (expensive)

            # update boundaries with EGM correction and SLR scenario
            if options.egm_file is not None:
                print('Updating boundary conditions with EGM correction in tile {:03d}'.format(tile['i']))
                egm_correction = cl.gdal_sample_points(y_boundaries, x_boundaries, options.egm_file)
                # add egm correction to boundary conditions
                boundaries += egm_correction

            # update boundaries with SLR corrections
            if options.sea_level_rise_map is not None:
                print('Updating boundary conditions sea level rise in tile {:03d}'.format(tile['i']))
                if isinstance(options.sea_level_rise_map,float):
                    slr_correction = options.sea_level_rise_map
                else:
                    slr_correction = cl.gdal_sample_points(y_boundaries, x_boundaries, options.sea_level_rise_map)
                # add sea level rise to boundary conditions
                boundaries += slr_correction

            # calculate x en y arrays
            cols, rows = tile['size']
            gt = tile['gt']
            x_tile = np.linspace(gt[0]+gt[1]/2, gt[0]+gt[1]/2+gt[1]*(cols-1), cols)
            y_tile = np.linspace(gt[3]+gt[5]/2, gt[3]+gt[5]/2+gt[5]*(rows-1), rows)

            # flood routine if boundary points in tile
            print('Computing flooding in tile {:03d}:'.format(tile['i']))
            flood = flood_tile(terrain, x_tile, y_tile, tile['i'],
                               boundaries, x_boundaries, y_boundaries,
                               ldd=ldd_map, resistance=options.resistance,
                               water_perc=water_perc, zero_resistance_waterp=options.zrw,
                               tempdir=options.tempdir, test=(options.test is not None),
                               nodatavalue=options.nodatavalue, dist_method=options.dist_method)

        # no boundary conditions
        else:

            flood = np.ones(np.shape(terrain))*options.nodatavalue  # new nodatavalue grid
            flood[terrain > -100] = 0  # zero flooding where terrain data

        ds_terrain = None  # close dataset
        del terrain, ds_terrain  # cleanup memory

    # no dem
    else:
        cols, rows = tile['size']
        flood = np.ones((rows, cols))*options.nodatavalue

    # mask out nodatavalue
    flood = np.ma.masked_equal(flood, options.nodatavalue)
    # cut relevant part and swap y-axis
    col_overlap_min, col_overlap_max, row_overlap_min, row_overlap_max = tile['overlap']
    if col_overlap_max == 0:
        col_overlap_max = -flood.shape[1]
    if row_overlap_max == 0:
        row_overlap_max = -flood.shape[0]
    flood_cut = flood[row_overlap_min:-row_overlap_max, col_overlap_min:-col_overlap_max]
    x, y = tile['x'], tile['y']
    x_cut = x[col_overlap_min:-col_overlap_max]
    y_cut = y[row_overlap_min:-row_overlap_max]

    # save to tiff
    gdal_writemap(tile['fn'], 'GTiff', x_cut, y_cut, flood_cut, options.nodatavalue, srs=options.srs)

    # done with tile
    seconds = float(time.time()-time0)
    print('finished flood routine for tile {:05d} with {:d} boundary conditions  in {:.2f} sec'.format(
    tile['i'], len(tile['bounds_select']), seconds))


def main():
    first_time = time.time()
    ### Read input arguments #####
    # parser = OptionParser()
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option('-q', '--quiet',
                      dest='verbose', default=True, action='store_false',
                      help='do not print status messages to stdout')
    parser.add_option('-i', '--ini', dest='inifile',
                      default='coastal_inun.ini', nargs=1,
                      help='ini configuration file')
    parser.add_option('-b', '--boundary',
                      nargs=1, dest='boundary',
                      help='boundary conditions file (NetCDF point time series file')
    parser.add_option('-v', '--boundary_variable',
                      nargs=1, dest='boundary_variable',
                      default='waterlevel',
                      help='variable name of boundary conditions')
    parser.add_option('-s', '--sea_level_rise',
                      dest='sea_level_rise_map', default='',
                      help='Sea level rise map (GeoTIFF)')
    parser.add_option('-g', '--subsidence',
                      dest='subsidence_map', default='',
                      help='Subsidence map (GeoTIFF)')
    parser.add_option('-d', '--destination',
                      dest='destination', default='',
                      help='Destination file')
    parser.add_option('-t', '--time',
                      dest='time', default='2010-01-01 00:00:00',
                      help='Time stamp of flood condition')
    parser.add_option('-x', '--test',
                      dest='test', default=None,
                      help='test specific tile number; report intermediate outputs')
    # for testing: tile_settings, tempdir, resistance and dist_method options as command line options.
    # if not set, these options are read from ini (default)
    parser.add_option('-y', '--tile_settings',
                      dest='tiles', default=None,
                      help='filename of JSON tile settings')
    parser.add_option('-m', '--dist_method',
                      dest='dist_method', default=None,
                      help="calculate distance along 'ldd' or use the 'eucledian' distance ")
    parser.add_option('-r', '--resistance',
                      dest='resistance', default=None,
                      help="decrease in water level as function of distance from coast [m/m]")
    parser.add_option('--zrw',
                      dest='zrw', default=None,
                      help="zero resistance water percentage thresho")
    parser.add_option('-z', '--tempdir',
                      dest='tempdir', default=None,
                      help="output directory for temporary data")
    parser.add_option('-w', '--nworkers',
                      dest='nworkers', default=cpu_count()-1,
                      help="number of parallel workers; if 1 it runs sequential") 
    (options, args) = parser.parse_args()

    if not os.path.exists(options.inifile):
        print('path to ini file cannot be found: {:s}'.format(options.inifile))
        sys.exit(1)

    # file names and directory bookkeeping
    options.destination = os.path.abspath(options.destination)
    options.dest_path = os.path.split(options.destination)[0]
    logfilename = options.destination[:-3] + '.log'

    # create dir if not exist
    if not os.path.isdir(options.dest_path):
        os.makedirs(options.dest_path)
    # delete old destination and log files
    else:
        if os.path.isfile(options.destination):
            os.unlink(options.destination)
        if os.path.isfile(logfilename):
            os.unlink(logfilename)

    # set up the logger
    logger, ch = cl.setlogger(logfilename, 'COASTAL_INUN', options.verbose)
    logger.info('$Id: coastal_inun.py 528 2018-06-19 08:41:05Z eilan_dk $')

    ### READ CONFIG FILE
    # open config-file
    config = cl.open_conf(options.inifile)

    # read settings

    options.dem_file = os.path.abspath(cl.configget(config, 'maps', 'dem_file', True))
    options.ldd_file = cl.configget(config, 'maps', 'ldd_file', None)
    if options.ldd_file is not None:
        options.ldd_file = os.path.abspath(options.ldd_file)
    options.water_perc_file = cl.configget(config, 'maps', 'water_perc_file', None)
    if options.water_perc_file is not None:
        options.water_perc_file = os.path.abspath(options.water_perc_file)
    options.egm_file = os.path.abspath(cl.configget(config, 'maps', 'egm_file', ''))
    options.x_var = cl.configget(config, 'boundary', 'x_var', 'station_x_coordinate')
    options.y_var = cl.configget(config, 'boundary', 'y_var', 'station_y_coordinate')
    options.x_tile = cl.configget(config, 'tiling', 'x_tile', 600, datatype='int')
    options.y_tile = cl.configget(config, 'tiling', 'y_tile', 600, datatype='int')
    options.x_overlap = cl.configget(config, 'tiling', 'x_overlap', 60, datatype='int')
    options.y_overlap = cl.configget(config, 'tiling', 'y_overlap', 60, datatype='int')
    if options.tiles is None:
        options.tiles = cl.configget(config, 'tiling', 'tiles', None)
        if options.tiles is not None:
            options.tiles = os.path.abspath(options.tiles)
    if options.resistance is None:
        options.resistance = cl.configget(config, 'flood_routine', 'resistance', 0.00050, datatype='float')
    else:
        options.resistance = float(options.resistance)
    if options.zrw is None:
        options.zrw = cl.configget(config, 'flood_routine', 'waterp_thresh', 1.0, datatype='float')
    else:
        options.zrw = float(options.zrw)
    if options.dist_method is None:
        options.dist_method = cl.configget(config, 'flood_routine', 'dist_method', 'eucledian')
    if options.tempdir is None:
        options.tempdir = os.path.abspath(cl.configget(config, 'flood_routine', 'tempdir', os.path.join(options.dest_path, 'temp{:s}'.format(str(uuid.uuid4())))))
    options.nodatavalue = cl.configget(config, 'flood_routine', 'nodatavalue', -9999, datatype='float')
    options.srs = osr.GetWellKnownGeogCSAsWKT(cl.configget(config, 'flood_routine', 'srs', 'EPSG:4326'))

    # required input
    if not options.destination:   # if destination is not given
        parser.error('destination not given')
    if not options.boundary:   # if boundary conditions argument is not given
        #options.boundary='global_etc_rp_database.nc'
        parser.error('boundary conditions not given')
    if not os.path.exists(options.dem_file):
        logger.error('path to dem file {:s} cannot be found'.format(options.dem_file))
        sys.exit(1)
    if options.dist_method not in ['ldd', 'eucledian']:
        logger.error("unknown value for distance method use 'ldd' or 'eucledian'")
        sys.exit(1)

    # check paths and set default to None if not given
    options.ldd_file = cl.check_input_fn(options.ldd_file, logger)
    options.egm_file = cl.check_input_fn(options.egm_file, logger)
    try:
        options.sea_level_rise_map = float(options.sea_level_rise_map) #constant SLR
    except:
        options.sea_level_rise_map = cl.check_input_fn(options.sea_level_rise_map, logger)
    options.subsidence_map = cl.check_input_fn(options.subsidence_map, logger)
    options.water_perc_file = cl.check_input_fn(options.water_perc_file, logger)

    # make sure tempdir is new empty folder
    if not os.path.isdir(options.tempdir):
        os.makedirs(options.tempdir)
    elif os.listdir(options.tempdir) == "":
        options.tempdir = options.tempdir  # do nothing
    else:
        n = 1
        options.tempdir = options.tempdir + '_{:03d}'.format(n)
        while os.path.isdir(options.tempdir):
            n += 1
            options.tempdir = options.tempdir[:-4] + '_{:03d}'.format(n)
        os.makedirs(options.tempdir)

    # write info to logger
    logger.info('Destination file: {:s}'.format(options.destination))
    logger.info('Temporary directory: {:s}'.format(options.tempdir))
    logger.info('Time of flood conditions: {:s}'.format(options.time))
    logger.info('DEM file: {:s}'.format(options.dem_file))
    logger.info('LDD file: {:s}'.format(options.ldd_file))
    logger.info('EGM file: {:s}'.format(options.egm_file))
    logger.info('Water mask file: {:s}'.format(options.water_perc_file))
    logger.info('Sea level rise map: {}'.format(options.sea_level_rise_map))
    logger.info('Subsidence map: {:s}'.format(options.subsidence_map))
    logger.info('Using tiling from json file: {:s}'.format(str(options.tiles is not None)))
    if options.tiles is None:
        logger.info('Columns per tile: {:d}'.format(options.x_tile))
        logger.info('Rows per tile: {:d}'.format(options.y_tile))
        logger.info('Columns overlap: {:d}'.format(options.x_overlap))
        logger.info('Rows overlap: {:d}'.format(options.y_overlap))
    logger.info('Flood resistance: {:.6f}'.format(options.resistance))
    logger.info('Distance method: {:s}'.format(options.dist_method))
    if options.test is None:
        logger.info('Running test: False')
    else:
        logger.info('Running test: True')

    #########################################################################
    # PREPARE TILES AND OUTPUT
    #########################################################################

    # add metadata from the section [metadata]
    metadata_global = {}
    meta_keys = config.options('metadata')
    for key in meta_keys:
        metadata_global[key] = config.get('metadata', key)
    # add a number of metadata variables that are mandatory
    metadata_global['config_file'] = os.path.abspath(options.inifile)
    metadata_var = {}
    metadata_var['units'] = 'm'
    metadata_var['standard_name'] = 'water_surface_height_above_reference_datum'
    metadata_var['long_name'] = 'Coastal flooding'
    metadata_var['comment'] = 'water_surface_reference_datum_altitude is given in file {:s}'.format(options.dem_file)

    # copy inifile to tempdir for reproducibility
    inifilename = options.destination[:-3] + '.ini'
    copyfile(options.inifile, inifilename)

    # Read extent from a GDAL compatible file
    try:
        x, y = cl.get_gdal_axes(options.dem_file, logging=logger)
    except:
        msg = 'Input file {:s} not a gdal compatible file'.format(options.dem_file)
        cl.close_with_error(logger, ch, msg)
        sys.exit(1)

    # open the variable with boundary conditions as preparation to read array parts
    try:
        with nc.Dataset(options.boundary, 'r') as a:
            # read history from boundary conditions
            try:
                history = a.history
            except:
                history = 'not provided'
            metadata_global['history'] = """Created by: $Id: coastal_inun.py 528 2018-06-19 08:41:05Z eilan_dk $,
                                            boundary conditions from {:s},\nhistory: {:s}""".format(
                os.path.abspath(options.boundary), history)
    except:
        msg = 'Input file {:s} not a netcdf compatible file'.format(options.boundary)
        cl.close_with_error(logger, ch, msg)
        sys.exit(1)

    # first -setup a NetCDF file
    nc_funcs.prepare_nc(options.destination, x, np.flipud(y),
                        [datetime.datetime.strptime(options.time, '%Y-%m-%d %H:%M:%S')],
                        metadata_global, units='Days since 1960-01-01 00:00:00')
    nc_funcs.append_nc(options.destination, 'inun', chunksizes=(1, min(options.y_tile, len(y)), min(options.x_tile, len(x))),
                       fill_value=options.nodatavalue, metadata=metadata_var)

    # read tile settings is json file given
    if options.tiles is not None:
        with open(options.tiles, 'r') as data_file:
            tile_list = json.load(data_file)
        logger.info('raster tile setting read from {:s}'.format(options.tiles))
    #  otherwise, start discretizing
    else:
        logger.info('discretizing raster...')
        tile_list = cl.discretize_raster_bounds(options.dem_file,
                                                options.x_tile, options.y_tile, options.x_overlap, options.y_overlap,
                                                options.boundary, options.x_var, options.y_var)
        # write output to tempdir
        tiles_fn = options.destination[:-3] + '_tiles.json'
        with open(tiles_fn, 'w') as outfile:
            json.dump(tile_list, outfile)

    # add options to tiles list
    tile_list_out = []
    for t in tile_list:
        t['options'] = options
        t['fn'] = os.path.join(options.tempdir, 'flood_{:05d}.tif'.format(t['i']))
        tile_list_out.append(t)
    tile_list = tile_list_out

    # discretizing finished
    n_tiles = len(tile_list)
    logger.info('raster discretized in {:d} tiles'.format(n_tiles))

    # find tiles with SRTM data and boundary conditions
    flood_tiles = [i for i, tile in enumerate(tile_list) if tile['has_DEM_data'] and tile['has_bounds']]
    logger.info('run {:d} tiles with SRTM data and boundary conditions'.format(len(flood_tiles)))

    # test with random subset of tiles
    if options.test is not None:
        tile_list = [tile for tile in tile_list if tile['i'] in [int(options.test)]]
        n_tiles = len(tile_list)
        logger.info('test with subset of {:d} tiles'.format(n_tiles))

    #########################################################################
    # PROCESS TILES
    #########################################################################

    # initialize multiprocessing parameters
    nworkers = np.min([int(options.nworkers), n_tiles])   # number of cores to work with
    logger.info('start processing tiles with {:d} parallel processes'.format(nworkers))
    time0 = time.time()

    # process tiles multicore
    if nworkers > 1:
       p = Pool(nworkers)
       try:
           p.map(process_tile, tile_list)
       except:
           traceback.print_exc()
       p.close()
    else:
       [process_tile(t) for t in tile_list] #[process_tile(t) for t in tile_list] 

    # ##############################################
    seconds = float(time.time()-time0)
    hours, seconds = seconds // 3600, seconds % 3600
    minutes, seconds = seconds // 60, seconds % 60
    logger.info('finished processing {:d} tiles in {:02d}:{:02d}:{:02d}'.format(
        n_tiles, int(hours), int(minutes), int(seconds)))

    #########################################################################
    # WRITE DATA TO NETCDF
    #########################################################################
    # open the prepared netCDF file for appending
    nc_obj = nc.Dataset(options.destination, 'a')
    nc_var = nc_obj.variables['inun']
    logger.info('start appending data to netcdf file')
    time0 = time.time()

    # read data from tempdir and append
    for tile in tile_list:
        # read data
        _, _, flood_cut, fill_val = gdal_readmap(tile['fn'], 'GTiff')
        flood_cut = np.flipud(flood_cut)
        if tile['start_row'] == 0:
            nc_var[0, -tile['end_row']:, tile['start_col']:tile['end_col']] = flood_cut
        else:
            nc_var[0, -tile['end_row']:-tile['start_row'], tile['start_col']:tile['end_col']] = flood_cut

    # now close nc file
    nc_obj.sync()
    nc_obj.close()

    seconds = float(time.time()-time0)
    hours, seconds = seconds // 3600, seconds % 3600
    minutes, seconds = seconds // 60, seconds % 60
    logger.info('finished writing {:d} tiles in {:02d}:{:02d}:{:02d}'.format(
        n_tiles, int(hours), int(minutes), int(seconds)))

    # cleanup temp dir
    if options.test is None:
        cl.cleanDir(options.tempdir, logger)

    # log total processing time
    seconds = float(time.time()-first_time)
    hours, seconds = seconds // 3600, seconds % 3600
    minutes, seconds = seconds // 60, seconds % 60
    logger.info('total processing time {:02d}:{:02d}:{:02d}'.format(
        int(hours), int(minutes), int(seconds)))

    # close logger
    logger, ch = cl.closeLogger(logger, ch)
    del logger, ch
    sys.exit(0)

if __name__ == "__main__":
    main()
