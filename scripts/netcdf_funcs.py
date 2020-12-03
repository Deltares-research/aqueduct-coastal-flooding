# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 08:11:19 2015

@author: winsemi
"""

### EXAMPLE FOR GFP COURSE ON NETCDF
import sys, os
import logging
import netCDF4 as nc
import numpy as np

def prepare_nc(nc_file, x, y, time_list, metadata, units='Days since 2012-01-01 00:00:00', calendar='gregorian', format="NETCDF4", zlib=True, logging=logging):
    """
    This function prepares a NetCDF file with given metadata, for a certain year, daily basis data
    The function assumes a gregorian calendar and a time unit 'Days since 2012-01-01 00:00:00'
    inputs:
        trg_file:     path to new netcdf file
        time_list:    list with times (in datetime format) which will give the time axis
        x:            xaxis
        y:            yaxis
        metadata:     dictionary with global attributes
        units:        time units to use in time axis
    """
    
    print('Setting up "' + nc_file + '"')
    startDayNr = nc.date2num(time_list[0], units=units, calendar=calendar)
    endDayNr   = nc.date2num(time_list[-1], units=units, calendar=calendar)
    time       = np.arange(startDayNr, endDayNr+1)
    try:
        nc_trg     = nc.Dataset(nc_file, 'w', format=format, zlib=zlib)
    except IOError:
        logging.error('Cannot write to {:s}'.format(nc_file))
        sys.exit(1)

    print('Setting up dimensions and attributes. lat: ' + str(len(y))+ " lon: " + str(len(x)))
    nc_trg.createDimension('time', 0)
    nc_trg.createDimension('lat', len(y))
    nc_trg.createDimension('lon', len(x))
    DateHour = nc_trg.createVariable('time', 'f8', ('time',))
    DateHour.units = units
    DateHour.calendar = calendar
    DateHour.standard_name = 'time'
    DateHour.long_name = 'time'
    DateHour.axis = 'T'
    DateHour[:] = time
    y_var   = nc_trg.createVariable('lat','f8',('lat',))
    y_var.standard_name = 'latitude'
    y_var.long_name = 'latitude'
    y_var.units = 'degrees_north'
    y_var.axis = 'Y'
    x_var = nc_trg.createVariable('lon','f8',('lon',))
    x_var.standard_name = 'longitude'
    x_var.long_name = 'longitude'
    x_var.units = 'degrees_east'
    x_var.axis = 'X'
    y_var[:] = y
    x_var[:] = x
    projection= nc_trg.createVariable('projection','c')
    projection.long_name = 'wgs84'
    projection.EPSG_code = 'EPSG:4326'
    projection.proj4_params = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    projection.grid_mapping_name = 'latitude_longitude'

    # now add all attributes from user-defined metadata
    for attr in metadata:
        nc_trg.setncattr(attr, metadata[attr])
    nc_trg.sync()
    nc_trg.close()

def append_nc(nc_file, var_name, dtype='f4', chunksizes=(1, 128, 128), fill_value=-9999, metadata={}, logging=logging):
    """
    Write a new (empty) variable to target NetCDF file. 
    input:
    ::

        nc_file:         NetCDF object, referring to target file
        var_name:       String, variable name of source NetCDF
        metadata:       dictionary of attributes, belonging to the variable. 

    """
    
    # add the variable
    try:
        nc_obj = nc.Dataset(nc_file, 'a')
    except IOError:
        logging.error('Cannot write to {:s}'.format(nc_file))
        sys.exit(1)

    variab = nc_obj.createVariable(var_name, dtype,
                                    ('time', 'lat', 'lon',),
                                    chunksizes=chunksizes,
                                    fill_value=fill_value,
                                    zlib=True)
    # add some general attributes usually used in lat lon data
    variab.coordinates   = 'lat lon'
    # if a attributes dictionary exists, then append attributes from this dictionary
    if metadata:
        for attribute in metadata:
            variab.setncattr(attribute, metadata[attribute])
    nc_obj.sync()
    nc_obj.close()

def get_netcdf_axes(filename, x_var='x', y_var='y', logging=logging):
    try:
        nc_obj = nc.Dataset(filename, 'r')
    except IOError:
        logging.error('Boundary condition file {:s} not existent or not a NetCDF file'.format(filename))
        sys.exit(1)
    try:
        x_coords = nc_obj.variables[x_var][:]
    except NameError:
        logging.error('Variable containing x-coordinates {:s} does not exist'.format(x_var))
        sys.exit(1)
    try:
        y_coords = nc_obj.variables[y_var][:]
    except NameError:
        logging.error('Variable containing y-coordinates {:s} does not exist'.format(y_var))
        sys.exit(1)
    return x_coords, y_coords
    nc_obj.close()