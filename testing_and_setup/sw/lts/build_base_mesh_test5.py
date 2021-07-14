#!/usr/bin/env python

import numpy as np
import argparse as ap
from mpas_tools.ocean import build_spherical_mesh
from mpas_tools.mesh.creation.util import lonlat2xyz
from numpy import radians, cos, sin, arcsin, sqrt
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


def cellWidthVsLatLon(coarse_res, fine_res, fine_radius, plots):
    """
    Create cell width array for this mesh on a regular latitude-longitude grid.

    Returns
    -------
    cellWidth : ndarray
        m x n array of cell width in km

    lon : ndarray
        longitude in degrees (length n and between -180 and 180)

    lat : ndarray
        longitude in degrees (length m and between -90 and 90)
    """

    ddeg = 0.5  # computed grid resolution, in degrees
    latCenter = 30.0  # center point in degrees
    lonCenter = 270.0  # center point in degrees

    # NOTE: this script produces a mesh that is already earth sized, 
    # so COMMENT the scaling on init in the sw core
    
    # coarseResolution was set to
    # 200 for convergence test in paper, 
    # 21 for CPU ratio test, 
    # 40 for speed up test
    # 250 to have a very coarse mesh for testing
    coarseResolution = coarse_res  # km

    # fineResolution was set to
    # was 50 for convergence tests in paper,
    # 7 for CPU ratio test,
    # 2.5 for speed up test,
    # 80 to have a very coarse mesh for testing
    fineResolution = fine_res  # km

    # fineRadius was set to
    # 6000 for convergence test in paper,
    # 1500 for CPU ratio test,
    # 400 for speed up test,
    # 1500 to have very coarse mesh for testing
    fineRadius = fine_radius  # km

    transitionWidth = 100 # km
    earthRadius = 6.371e3 # km
    refinementFactor = 2
    
    lat = np.arange(-90, 90.01, ddeg)
    lon = np.arange(-180, 180.01, ddeg)
    latCenter = np.deg2rad(latCenter)
    lonCenter = np.deg2rad(lonCenter)

    # create meshed 2D grid in radians
    lonGrid, latGrid = np.meshgrid(np.deg2rad(lon), np.deg2rad(lat))

    # Halversine formula for distance
    distance =  ( np.sin((latGrid - latCenter) * 0.5) ** 2
                + (np.cos(latCenter)*np.cos(latGrid)
                * np.sin((lonGrid - lonCenter) * 0.5) ** 2) )
    distance = ( 2.0 * earthRadius
               * np.arctan2(np.sqrt(distance), np.sqrt(1.0 - distance)) )

    tanhDistance = 0.5 * ( 1 + np.tanh((fineRadius - distance) 
                                       / transitionWidth) )
     
    cellWidth = ( coarseResolution + (fineResolution - coarseResolution)
                * tanhDistance )

    if plots:
        varList = ['latGrid', 'lonGrid', 'distance', 'tanhDistance',
                   'cellWidth']
        fig = plt.gcf()
        plt.clf()
        fig.set_size_inches(20.0, 20.0)
        iPlt = 1
        for varName in varList:
            plt.subplot(3,2,iPlt)
            plt.imshow(vars()[varName])
            iPlt += 1
            plt.title(varName)
            plt.xlabel('lon index')
            plt.ylabel('lat index')
            plt.colorbar()
        plt.savefig('cellWidth.png')

    return cellWidth, lon, lat


def main(output_file, coarse_res, fine_res, fine_radius, plots):
    """
    Python script to build spatial mesh for test case 5 of Williamson et al.

    Run `./build_base_mesh_for_test5.py --help` for usage information.
    """

    cellWidth, lon, lat = cellWidthVsLatLon(coarse_res, fine_res, fine_radius, 
                                            plots)
    build_spherical_mesh(cellWidth, lon, lat, out_filename=output_file)


if __name__ == '__main__':
    # If this is the primary module, run main

    parser = ap.ArgumentParser(description='Python script to build spatial mesh \
                                        for test case 5 from Williamson et al.')

    parser.add_argument('-o', '--output-file', dest='output_file', type=str,
                    default='base_mesh.nc',
                    help='Name for output file. Default is `base_mesh.nc`.')

    parser.add_argument('-c', '--coarse-res', dest='coarse_res', type=float,
                        default=40,
                        help='Size in km for coarse cells. Default is \
                        40.')
    
    parser.add_argument('-f', '--fine-res', dest='fine_res', type=float,
                        default=2.5,
                        help='Size in km for coarse cells. Default is \
                        2.5.')

    parser.add_argument('-s', '--fine-radius', dest='fine_radius', type=float, 
                        default=400,
                        help='Radius around the mountain to place cells with \
                        the fine resolution. Default is 400.')

    parser.add_argument('-p', '--plots', dest='plots', action="store_true",
                    help='Produce plots of the mesh using matplotlib.')

    args = parser.parse_args()


    main(args.output_file,
         args.coarse_res,
         args.fine_res,
         args.fine_radius,
         args.plots)

