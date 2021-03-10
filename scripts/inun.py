# -*- coding: utf-8 -*-
"""
Created on 2016-March-24
@author: Dirk Eilander (dirk.eilander@deltares.nl)

$Id: inun.py 528 2018-06-19 08:41:05Z eilan_dk $
$Date: 2018-06-19 10:41:05 +0200 (Tue, 19 Jun 2018) $
$Author: eilan_dk $
$Revision: 528 $
$HeadURL: https://repos.deltares.nl/repos/aqueduct/trunk/python/coastal_inundation/inun_clean/scripts/inun.py $
$Keywords: $

"""

import numpy as np
import pcraster as pcr
from gdal_writemap import gdal_writemap
import coastal_inun_lib as cl
import os
# from scipy.interpolate import NearestNDInterpolator
import warnings
warnings.filterwarnings("ignore", '.*val in CSF.*', )  # ignore read error PCRaster

def flood_tile(terrain, x_terrain, y_terrain, itile,
               boundaries, x_boundaries, y_boundaries, ldd=None,
               resistance=0.05, water_perc=None, zero_resistance_waterp=1.0,
               tempdir=r'd:/temp/inun', test=False,
               nodatavalue=-9999, dist_method='ldd'):
    """
    function for rapid coast flood mapping based on topography and resistance

    :param terrain:         masked 2d numpy array with elevation data
    :param x_terrain:       1d numpy array with x coordinates
    :param y_terrain:       1d numpy array with y coordinates
    :param boundaries:      1d numpy array water level boundary conditions
    :param x_boundaries:    1d numpy array with x coordinates boundary conditions
    :param y_boundaries:    1d numpy array with y coordinates boundary conditions
    :param ldd              2 numpy array with lld data
                            if NONE an ldd is calculated based on the DEM file (default=None)
    :param resistance:      resistance to flooding (decrease in water depth over ldd distance) [km-1]
    :param tempdir:         directory to save temporary pcRaster clone file
    :param test:            if True, some intermediate steps are saved and clone map is not removed

    :return:                2d numpy array with flood depth [m]
    """
    # create geo-transform and get pixel coordinates of boundaries
    if dist_method=='ldd':
        create_ldd = True
    else:
        create_ldd = False
    # import pdb; pdb.set_trace()
    boundaries = np.array(boundaries).flatten()
    gt = cl.get_gt(x_terrain, y_terrain)
    ids = np.arange(1, np.size(boundaries)+1)  # ids start at one!
    locs_px = [cl.MapToPixel(x_b, y_b, gt) for x_b, y_b in zip(x_boundaries, y_boundaries)]

    # pcr set clone
    # and calc ldd if ldd is None OR translate to pcr map if a 2d array is given
    print('Preprocessing dem')
    dem, ldd = pcr_preprocess(terrain, x_terrain, y_terrain, itile, tempdir,
                                ldd_in=ldd, test=test, create_ldd=create_ldd)
    points, _ = val2pcrmap(np.shape(terrain), locs_px, ids)

    # cell resolution in km at latitude degree
    cellres_array = cl.distance_on_unit_sphere(np.array(y_terrain), np.abs(y_terrain[1]-y_terrain[0]))
    cellres_np = np.tile(np.array(cellres_array), (len(x_terrain), 1)).T
    # cellres = pcr.numpy2pcr(pcr.Scalar, cellres_np, 0) x an y axis get mixed up using this function??
    # work around -> save map to disk using gdal_write and read back...
    # TODO: fix work around for mixing up axis by numpy2pcr
    fn = os.path.join(tempdir, '_{:03d}_cellres1.map'.format(itile))
    gdal_writemap(fn, 'PCRaster', x_terrain, y_terrain, cellres_np, 0)
    cellres = pcr.readmap(fn)
    if not test:
        os.unlink(fn)
    if os.path.isfile(fn + '.aux.xml'):
        os.unlink(fn + '.aux.xml')

    # water_perc
    if water_perc is not None:
        fn_wperc = os.path.join(tempdir, '_{:03d}_wperc.map'.format(itile))
        water_perc_np = np.copy(water_perc)
        gdal_writemap(fn_wperc, 'PCRaster', x_terrain, y_terrain, water_perc, -9999)
        water_perc = pcr.readmap(fn_wperc)
        # cleanup
        if not test:
            os.unlink(fn_wperc)
        if os.path.isfile(fn_wperc + '.aux.xml'):
            os.unlink(fn_wperc+'.aux.xml')

    # find coastal pixels
    ids_coastline, ids_sea = pcr_coast(dem, points)  # returns 2D array with ones at coastal cells

    # find points that have been mapped
    ids_coastline_np = pcr.pcr2numpy(pcr.cover(ids_coastline, 0), 0)
    ids_projected = np.unique(ids_coastline_np).astype(np.int).tolist()
    idx_mapped = np.asarray([idx for idx, i in enumerate(ids) if i in ids_projected])
    if idx_mapped.size > 0:
        h_bounds = boundaries[idx_mapped]
        ids0 = ids[idx_mapped]
        print('Inundating tile')
        # calculate flood depth with resistance along ldd from coastal pixel with boundary
        flood_depth, dist2coast, dem_adjust = pcr_inun(dem, ids0, h_bounds, ids_coastline,
                                                    resistance, water_perc, zero_resistance_waterp,
                                                    cellres, dist_method, ldd)
    else:
        print('Unable to map boundary conditions. check intermediate outputs')
        flood_depth = pcr.scalar(-9999)
        dem_adjust, dist2coast = None, None
        test = True
    # if in test mode report some intermediate steps
    if test:
        if idx_mapped.size > 0:
            np.savetxt(os.path.join(tempdir, '_{:03d}_bounds.csv'.format(itile)), 
                np.stack([idx_mapped+1, x_boundaries[idx_mapped], y_boundaries[idx_mapped], h_bounds], axis=1), delimiter=",", fmt='%.4f')
        if dist2coast is not None:
            pcr.report(dist2coast * pcr.scalar(resistance),
                        os.path.join(tempdir, '_{:03d}_resistance.map'.format(itile)))
        if dem_adjust is not None:
            pcr.report(dem_adjust, os.path.join(tempdir, '_{:03d}_dem_adjust.map'.format(itile)))
        if ldd is not None:
            pcr.report(ldd, os.path.join(tempdir, '_{:03d}_ldd.map'.format(itile)))
        pcr.report(pcr.ifthen(flood_depth > 0, flood_depth),
                   os.path.join(tempdir, '_{:03d}_flood_depth.map'.format(itile)))
        pcr.report(pcr.nominal(pcr.ifthen(ids_coastline > 0, ids_coastline)),
                   os.path.join(tempdir, '_{:03d}_ids_coastline.map'.format(itile)))
        pcr.report(pcr.nominal(pcr.ifthen(pcr.scalar(ids_sea) > 0, ids_sea)),
                   os.path.join(tempdir, '_{:03d}_ids_sea.map'.format(itile)))
        pcr.report(pcr.nominal(pcr.ifthen(pcr.scalar(points) > 0, points)),
                   os.path.join(tempdir, '_{:03d}_ids.map'.format(itile)))
    # translate flood depth pcraster map to numpy 2d array and return
    flood_depth_np = pcr.pcr2numpy(flood_depth, terrain.fill_value.item()).astype(np.float32)  # 2numpy 2d array
    # set dem nodata and permanent water to nodata values in output
    # flood_depth_np[np.logical_or(flood_depth_np == terrain.fill_value.item(), water_perc_np >= float(zero_resistance_waterp))] = nodatavalue
    flood_depth_np[flood_depth_np == terrain.fill_value.item()] = nodatavalue
    flood_depth_np = np.ma.masked_equal(flood_depth_np, nodatavalue, copy=True)

    return flood_depth_np


def pcr_inun(dem, ids, h_bounds, ids_coastline,
                resistance=0.,  water_perc=None, zero_resistance_waterp=1.0,
                cellres=1, dist_method='eucledian', ldd=None):
    """ planar inundation routine per segment

    :param dem:             pcr dem
    :param ids:             local ids of boundary conditions, starting a 1 (not zero!)
    :param h_bounds:        water level boundary at diva segment
    :param ids_coastline:   pcraster map with coastal segments ids
    :param resistance:      constant or pcrmap unit km-1; (default 0: no resistance is calculated)
    :param cellres:         cell resolution in km, varies with latitude degrees
    :param ldd:             pcraster map with local drainage direction to calculate resistance along ldd;
                            if None (default) resistance is calculated using 'plain' nearest neighbour


    :return:                pcrmap with flood depth
    """
    pcr.setglobaloption("unitcell")
    if resistance > 0:
        coastline = pcr.cover(pcr.ifthenelse(pcr.scalar(ids_coastline) > 0, pcr.boolean(1), 0), pcr.boolean(0))
        mask = pcr.ifthen(dem > -9999, pcr.scalar(1))
        if dist_method == 'ldd':
            # Distance to coast along ldd
            dist2coast0 = pcr.ldddist(ldd, coastline, cellres)
            # find edge of area with distances -> water divide
            dist2coast_mask = pcr.cover(pcr.ifthenelse(dist2coast0 > 0, pcr.boolean(0), pcr.boolean(1)),
                                        pcr.boolean(1))
            start = pcr.ifthenelse(
                ((pcr.window4total(pcr.scalar(dist2coast_mask)) > 0) & (dist2coast_mask == pcr.boolean(0))) |
                coastline,
                pcr.boolean(1), pcr.boolean(0))
            # continue distance beyond water divide with eucledian dist
            dist2coast1 = pcr.spread(start, dist2coast0, cellres*mask)
            dist2coast = pcr.ifthenelse(dist2coast_mask, dist2coast1, dist2coast0)

        elif dist_method == 'eucledian':
            # dist to coast using nearest neighbor
            if water_perc is None:
                dist2coast = pcr.spread(coastline, 0, cellres*mask)
            else:
                # zero resistance for cells with water_perc >= zero_resistance_waterp
                zrw = float(zero_resistance_waterp)
                water_perc = pcr.ifthenelse(water_perc >= zrw,
                                            pcr.scalar(1),
                                            water_perc / zrw)
                dist2coast = pcr.spread(coastline, 0, cellres*mask*(1 - water_perc))

        dem_adjust = dem + pcr.cover(dist2coast, 0) * pcr.scalar(resistance)   # raise the elevation using a damping factor
    else:
        dem_adjust = dem
        dist2coast = pcr.scalar(1)

    fld_depth = pcr.ifthen(dem > -9999, pcr.scalar(0))

    for i, h in zip(ids, h_bounds):
        coast_segment = pcr.ifthenelse(ids_coastline == int(i), pcr.boolean(1), pcr.boolean(0))
        # find area below flood_level
        fld_prone = pcr.ifthenelse(dem_adjust <= pcr.scalar(float(h)), pcr.boolean(1), pcr.boolean(0))
        # make contiguous groups of cells which are below flood level
        fld_clump = pcr.clump(fld_prone)
        # find flooded area connected to diva segment
        fld_coast = pcr.ifthenelse(pcr.areamaximum(pcr.scalar(fld_clump) * pcr.scalar(coast_segment), fld_clump) > 0,
                                   pcr.boolean(1), pcr.boolean(0))
        # get max fld depth map
        fld_depth = pcr.max(fld_depth, pcr.ifthenelse(fld_coast, pcr.scalar(pcr.scalar(float(h)) - dem_adjust), 0))

    return fld_depth, dist2coast, dem_adjust


def pcr_preprocess(dem_in, x, y, itile, tempdir, ldd_in=None,
                    test=False, create_ldd=True):
    """
    function to set pcr clone and translate DEM (terrain) and ldd numpy 2d arrays to pcr maps

    :param terrain:     masked numpy 2d array with elevation data
    :param x:           numpy 1d array with x coordinates of elevation grid
    :param y:           numpy 1d array with Y coordinates of elevation grid
    :param tempdir:     string with directory to temporary save clone pcrmap
    :param ldd:         numpy 2d array with ldd grid, make sure it uses the pcrmap definition of ldd
    :param test:        if True do not remove clone maps

    :return:            pcr maps for dem and ldd
    """
    # create clone in temp_dir
    fn_clone = os.path.join(tempdir, '_{:03d}_dem.map'.format(itile))
    cl.makeDir(tempdir)  # make dir if not exist
    # DEM
    gdal_writemap(fn_clone, 'PCRaster', x, y, dem_in, -9999)  # note: missing value needs conversion in python item
    pcr.setclone(fn_clone)
    pcr.setglobaloption("unitcell")
    dem = pcr.readmap(fn_clone)
    # cleanup
    if not test:
        os.unlink(fn_clone)  # cleanup clone file
    os.unlink(fn_clone+'.aux.xml')
    # LDD
    if create_ldd:
        if ldd_in is None:
            print('Calculating LDD')
            ldd = pcr.lddcreate(dem, 1E31, 1E31, 1E31, 1E31)
        else:
            # TODO note that np.nan is default NoDataValue for ldd file. check this when reading data
            # TODO: x and y axis got mixed up when translating with numpy2pcr. check if it works here!
            # in process_tile function in coastal_inun.py
            ldd = pcr.lddrepair(pcr.ldd(pcr.numpy2pcr(pcr.Ldd, ldd_in, np.nan)))
    else:
        ldd = None
    return dem, ldd


def val2pcrmap(dims, locs_px, val):

    # make map with points
    points = np.zeros(dims)
    for xy, i in zip(locs_px, val):
        x, y = xy
        points[y, x] = i
    return pcr.numpy2pcr(pcr.Nominal, points, -9999), points


def pcr_coast(dem, points):
    """  project points to coast with nearest neighbourhood
    finds coastal cells based on dem with NoDataValues and the locations of boundary conditions at sea
    using pcr spread the a nearest neighbor interpolation of the point ids is done for coastal cells

    :param dem: pcr dem
    :param points: pcrmap with ids in cells

    :returns: pcr map with location ids projected to coastline
    """
    # clump areas based on NoDataValues in dem
    dem_NoDataValues = pcr.cover(pcr.ifthenelse(dem > -9999, pcr.boolean(0), pcr.boolean(1)), pcr.boolean(1))
    # find number of boundary conditions in area where dem_novalue
    pcr.setglobaloption("nondiagonal")  # only top, bottom, left, right
    area_nbounds = pcr.areatotal(pcr.scalar(points), pcr.clump(dem_NoDataValues)) * pcr.scalar(dem_NoDataValues)
    pcr.setglobaloption("diagonal")  # diagonal again
    # make sea (True) and land (False) mask
    if np.any(pcr.pcr2numpy(area_nbounds,-9999) > 0):
        sea = pcr.ifthenelse(area_nbounds > 0, pcr.boolean(1), pcr.boolean(0))
    else:
        sea = dem_NoDataValues
    # find coast based on sea in neighboring cells and at land (sea = 0)
    coast = pcr.ifthenelse((pcr.window4total(pcr.scalar(sea)) > pcr.scalar(0)) & (sea == pcr.boolean(0)),
                           pcr.boolean(1), pcr.boolean(0))

    # move points to nearest sea cell(s)
    point_dist = pcr.ifthenelse(sea, pcr.spread(points, 0, 1), 1E31)  # distance from each point for sea cells
    nnpoints = pcr.ifthenelse(sea, pcr.spreadzone(points, 0, 1), 0)  # closest point for sea cells
    dist2sea = pcr.areaminimum(point_dist, nnpoints)  # shortest distance to each point to sea
    points_in_sea = pcr.nominal(pcr.ifthenelse(dist2sea == point_dist, nnpoints, 0))  # map points to nearest sea cell

    # map point at sea to coastline according to shortest distance over sea
    res = pcr.ifthenelse((pcr.scalar(sea) + pcr.scalar(coast)) >= 1, pcr.scalar(1), 1E31)  # mask out non sea or coast cells
    ids_coastline = pcr.scalar(pcr.spreadzone(points_in_sea, 0, res)) * pcr.scalar(coast)

    return ids_coastline, points_in_sea

# deprecated!
# def ids2coastNN(locs_px, vals, coast):
#     """
#     deprecated! replace with ids coastline -> pcr.scalar(pcr.spreadzone(points, 0, 1)) * coast
#     """
#     # find coastal cells
#     coast_np = pcr.pcr2numpy(coast, -9999).astype(np.float32)
#     py, px = np.where(coast_np == 1)
#     coastal_cells = [(int(x), int(y)) for x, y in zip(px, py)]
#
#     # NN interpolation
#     ip = NearestNDInterpolator(locs_px, vals)
#     coastal_id = [ip(cc) for cc in coastal_cells]
#
#     # plot on map
#     ids_coastline_np = np.zeros(np.shape(coast_np)).astype(np.float32)
#     for xy, i in zip(coastal_cells, coastal_id):
#         x, y = xy
#         ids_coastline_np[y, x] = i
#
#     return pcr.numpy2pcr(pcr.Scalar, ids_coastline_np, np.nan), ids_coastline_np
