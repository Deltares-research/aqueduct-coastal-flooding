#!/usr/bin/env python
"""
Created on Tue Jul 16 11:41:58 2013

@author: Hessel Winsemius

$Id: gdal_writemap.py 445 2018-03-28 14:39:49Z eilan_dk $
$Date: 2018-03-28 16:39:49 +0200 (Wed, 28 Mar 2018) $
$Author: eilan_dk $
$Revision: 445 $
$HeadURL: https://repos.deltares.nl/repos/aqueduct/trunk/python/coastal_inundation/inun_clean/scripts/gdal_writemap.py $
$Keywords: $

"""

import numpy as np
from osgeo import gdal
import os
import sys
import logging

def gdal_writemap(file_name, file_format, x, y, data, fill_val, zlib=False,
                  gdal_type=gdal.GDT_Float32, resolution=None, srs=None):
    """ Write geographical file from numpy array
    Dependencies are osgeo.gdal and numpy
    Input:
        file_name: -- string: reference path to GDAL-compatible file
        file_format: -- string: file format according to GDAL acronym
        (see http://www.gdal.org/formats_list.html)
        x: -- 1D np-array: x-axis, or (if only one value), top-left x-coordinate
        y: -- 1D np-array: y-axis, or (if only one value), top-left y-coordinate
        data: -- 2D np-array: raster data
        fill_val: -- float: fill value
        --------------------------------
    optional inputs:
        zlib=False: -- boolean: determines if output file should be internally 
                        zipped or not
        gdal_type=gdal.GDT_Float32: -- gdal data type to write
        resolution=None: -- resolution of dataset, only needed if x and y are given as upperleft coordinates
        srs=None: -- projection object (imported by osgeo.osr)
    """
    # make the geotransform
    # Give georeferences
    if hasattr(x, '__len__'):
        # x is the full axes
        xul = x[0]-(x[1]-x[0])/2
        xres = x[1]-x[0]
    else:
        # x is the top-left corner
        xul = x
        xres = resolution
    if hasattr(y, '__len__'):
        # y is the full axes
        yul = y[0]+(y[0]-y[1])/2
        yres = y[1]-y[0]
    else:
        # y is the top-left corner
        yul = y
        yres = -resolution
    geotrans = [xul, xres, 0, yul, 0, yres]
    
    gdal.AllRegister()
    driver1 = gdal.GetDriverByName('GTiff')
    driver2 = gdal.GetDriverByName(file_format)
    # Processing
    temp_file_name = str('{:s}.tif').format(file_name)
    logging.info(str('Writing to temporary file {:s}').format(temp_file_name))
    if zlib:
        TempDataset = driver1.Create(temp_file_name, data.shape[1],
                                     data.shape[0], 1, gdal_type,
                                     ['COMPRESS=DEFLATE'])
    else:
        TempDataset = driver1.Create(temp_file_name, data.shape[1],
                                     data.shape[0], 1, gdal_type)
    TempDataset.SetGeoTransform(geotrans)
    if srs is not None:
        if type(srs) == str:
            TempDataset.SetProjection(srs)
        else:
            TempDataset.SetProjection(srs.ExportToWkt())

    # get rasterband entry
    TempBand = TempDataset.GetRasterBand(1)
    # fill rasterband with array
    TempBand.WriteArray(data, 0, 0)
    TempBand.FlushCache()
    TempBand.SetNoDataValue(fill_val)
    # Create data to write to correct format (supported by 'CreateCopy')
    logging.info(str('Writing to {:s}').format(file_name))
    if zlib:
        driver2.CreateCopy(file_name, TempDataset, 0, ['COMPRESS=DEFLATE'])
    else:
        driver2.CreateCopy(file_name, TempDataset, 0)
    TempDataset = None
    os.remove(temp_file_name)
