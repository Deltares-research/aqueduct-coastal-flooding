#!/usr/bin/env python
"""
Created on Tue Jul 16 11:41:58 2013

@author: Hessel Winsemius

$Id: gdal_warp.py 445 2018-03-28 14:39:49Z eilan_dk $
$Date: 2018-03-28 16:39:49 +0200 (Wed, 28 Mar 2018) $
$Author: eilan_dk $
$Revision: 445 $
$HeadURL: https://repos.deltares.nl/repos/aqueduct/trunk/python/coastal_inundation/inun_clean/scripts/gdal_warp.py $
$Keywords: $

"""

import numpy as np
from osgeo import gdal, gdalconst, osr
from scipy.interpolate import griddata
import pdb

def gdal_warp(src_filename, clone_filename, dst_filename, gdal_type=gdalconst.GDT_Float32,
              gdal_interp=gdalconst.GRA_Cubic, format='GTiff', ds_in=None, override_src_proj=None):
    """
    Equivalent of the gdalwarp executable, commonly used on command line.
    The function prepares from a source file, a new file, that has the same 
    extent and projection as a clone file.
    The clone file should contain the correct projection. 
    The same projection will then be produced for the target file.
    If the clone does not have a projection, EPSG:4326 (i.e. WGS 1984 lat-lon)
    will be assumed.

    :param src_filename: string - file with data that will be warped
    :param clone_filename: string - containing clone file (with projection information)
    :param dst_filename: string - destination file (will have the same extent/projection as clone)
    :param gdal_type: - data type to use for output file (default=gdalconst.GDT_Float32)
    :param gdal_interp: - interpolation type used (default=gdalconst.GRA_Bilinear)
    :param format: - GDAL data format to return (default='GTiff')
    :return: No parameters returned, instead a file is prepared
    """
    if ds_in is None:
        src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
    else:
        src = ds_in
    src_proj = src.GetProjection()
    if override_src_proj is not None:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(override_src_proj)
        src_proj = srs.ExportToWkt()        
    src_nodata = src.GetRasterBand(1).GetNoDataValue()

    # We want a section of source that matches this:
    if type(clone_filename) == str:
        clone_ds = gdal.Open(clone_filename, gdalconst.GA_ReadOnly)
    else:
        clone_ds = clone_filename
    clone_proj = clone_ds.GetProjection()
    if not clone_proj:
        # assume a WGS 1984 projection
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        clone_proj = srs.ExportToWkt()
    clone_geotrans = clone_ds.GetGeoTransform()
    wide = clone_ds.RasterXSize
    high = clone_ds.RasterYSize

    # Output / destination
    dst_mem = gdal.GetDriverByName('MEM').Create('', wide, high, 1, gdal_type)
    dst_mem.SetGeoTransform(clone_geotrans)
    dst_mem.SetProjection(clone_proj)
    if not(src_nodata is None):
        dst_mem.GetRasterBand(1).SetNoDataValue(src_nodata)
    # fill dataset with NoDataValue rather than zeros (default)
    dst_mem.GetRasterBand(1).Fill(src_nodata, 0)

    # reproject
    gdal.ReprojectImage(src, dst_mem, src_proj, clone_proj, gdal_interp)

    if format == 'MEM':
        return dst_mem
    else:
        # retrieve numpy array of interpolated values
        # write to final file in the chosen file format
        gdal.GetDriverByName(format).CreateCopy(dst_filename, dst_mem, 0)
