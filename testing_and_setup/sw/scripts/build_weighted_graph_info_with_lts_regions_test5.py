#!/usr/bin/env python

import os
import shutil
import numpy as np
import xarray as xr
import argparse as ap
import math
import time


parser = ap.ArgumentParser(description='Python script to label the cells \
        from a base mesh for test case 5 from \
        Williamson et al. for LTS. This script \
        MODIFIES the graph.info file specified in its \
        command-line arguments. Note that \
        before running this script on a base \
        mesh, one should convert it to a valid \
        MPAS mesh with `MpasMeshConverter.x \
        <base_mesh>` -- doing This also produces \
        a corresponding graph.info file.')

parser.add_argument('-b', '--base-mesh', dest='base_mesh',
        default='base_mesh.nc',
        help='File containing the base mesh. Default is \
                `base_mesh.nc`.')

parser.add_argument('-g', '--graph-info', dest='graph_info',
        default='graph.info',
        help='The graph.info file corresponding to the base \
                mesh. Default is `graph.info`.')

parser.add_argument('-c', '--coarse-region-dist', 
        dest='coarse_region_dist', default=0.55,
        help='Cells more than this distance way from the mountain at \
                the north pole will be part of the coarse region.')

parser.add_argument('--lts2', action='store_true',
        help='Prepare the mesh for LTS2 rather than LTS3 which is the \
                default.')

args = parser.parse_args()


def main(base_mesh, graph_info, coarse_region_dist, lts2):
    timeStart = time.time()

    ds = xr.open_dataset(base_mesh)

    # Set to 1 for LTS2, 3 for LTS3
    if lts2:
        nLTSHalos = 1
    else:
        nLTSHalos = 3

    # THIS FLAG HAS TO BE THE SAME AS THE ONE SET IN 
    # THE `mpas_init_LTS` ROUTINE in `src/core_sw/mpas_sw_test_cases.F`
    moreCellsOnInterface = 1

    nLTSHalosCopy = nLTSHalos
    if moreCellsOnInterface == 1:
        # THE NUMBER SET HERE IS THE NUMBER OF EXTRA LTS HALO LAYERS 
        # AND HAS TO BE THE SAME SET IN init_LTS 
        nLTSHalosCopy = 40 

    nCells = ds['nCells'].size
    nEdges = ds['nEdges'].size
    areaCell = ds['areaCell'].values
    nEdgesOnCell = ds['nEdgesOnCell'].values
    cellsOnEdge = ds['cellsOnEdge'].values
    edgesOnCell = ds['edgesOnCell'].values
    latCell = ds['latCell']
    lonCell = ds['lonCell']
    latPoint = math.pi / 6.0
    lonPoint = 3.0 * math.pi * 0.5

    LTSRegion = [1] * nCells
    LTSRegionLocal = [1] * nCells

    for iCell in range(0, nCells):
        arg1 = np.sqrt(np.sin(0.5 * (latPoint - latCell[iCell])) ** 2 
                + np.cos(latCell[iCell]) * np.cos(latPoint) 
                * np.sin(0.5 * (lonPoint - lonCell[iCell])) ** 2)
        sphereDistance = 2 * np.arcsin(arg1)

        #CAREFUL: THIS HAS TO MATCH WHAT'S ON `src/core_sw/mpas_sw_test_cases.F` 
        if sphereDistance > coarse_region_dist: 
            LTSRegion[iCell] = 2
            LTSRegionLocal[iCell] = 2

    minMaxLTSRegion = [1] * 2
    minMaxLTSRegion[0] = min(LTSRegionLocal)
    minMaxLTSRegion[1] = max(LTSRegionLocal)

    for iEdge in range(0, nEdges):
        cell1 = cellsOnEdge[iEdge, 0] - 1 
        cell2 = cellsOnEdge[iEdge, 1] - 1

        if ( (LTSRegionLocal[cell1] == minMaxLTSRegion[0]) 
                and (LTSRegionLocal[cell2] == minMaxLTSRegion[1]) ): 
            LTSRegionLocal[cell1] = minMaxLTSRegion[0] + 2
            LTSRegionLocal[cell2] = minMaxLTSRegion[1] + 2

        elif ( (LTSRegionLocal[cell2] == minMaxLTSRegion[0]) 
                and (LTSRegionLocal[cell1] == minMaxLTSRegion[1]) ):
            LTSRegionLocal[cell2] = minMaxLTSRegion[0] + 2
            LTSRegionLocal[cell1] = minMaxLTSRegion[1] + 2

        elif ( (LTSRegionLocal[cell1] == minMaxLTSRegion[0] + 2) 
                and (LTSRegionLocal[cell2] == minMaxLTSRegion[1]) ):
            LTSRegionLocal[cell2] = minMaxLTSRegion[1] + 2

        elif ( (LTSRegionLocal[cell2] == minMaxLTSRegion[0] + 2) 
                and (LTSRegionLocal[cell1] == minMaxLTSRegion[1]) ):
            LTSRegionLocal[cell1] = minMaxLTSRegion[1] + 2

        elif ( (LTSRegionLocal[cell1] == minMaxLTSRegion[0]) 
                and (LTSRegionLocal[cell2] == minMaxLTSRegion[1] + 2) ):
            LTSRegionLocal[cell1] = minMaxLTSRegion[0] + 2

        elif ( (LTSRegionLocal[cell2] == minMaxLTSRegion[0]) 
                and (LTSRegionLocal[cell1] == minMaxLTSRegion[1] + 2) ):
            LTSRegionLocal[cell2] = minMaxLTSRegion[0] + 2

    for iRegion in range(0, 2):
        for iHalo in range(0, nLTSHalosCopy-1):
            for iCell in range(0, nCells):
                if (LTSRegionLocal[iCell] == ( minMaxLTSRegion[iRegion] 
                    + 2 * (iHalo+1)) ):
                    for i in range(0, nEdgesOnCell[iCell]):
                        iEdge = edgesOnCell[iCell, i] - 1
                        cell1 = cellsOnEdge[iEdge, 0] - 1
                        cell2 = cellsOnEdge[iEdge, 1] - 1

                        if (LTSRegionLocal[cell1] == minMaxLTSRegion[iRegion]):
                            LTSRegionLocal[cell1] = ( LTSRegionLocal[cell1] 
                                +  2 * (iHalo + 2) )
                        elif (LTSRegionLocal[cell2] == minMaxLTSRegion[iRegion]):
                            LTSRegionLocal[cell2] = ( LTSRegionLocal[cell2] 
                                +  2 * (iHalo + 2) )

    if moreCellsOnInterface == 1:
        # this means we are using LTS2
        if (nLTSHalos == 1): 
            for iCell in range(0, nCells):
                
                # if it's odd and not 1 (fine) then must
                # be 3 (interface layer 1)
                if ( (LTSRegionLocal[iCell] % 2 == 1) 
                        and (LTSRegionLocal[iCell] != 1) ):
                    LTSRegionLocal[iCell] = 3
                
                # if it's even and not 2 (coarse) then 
                # must be 4 (interface layer 2)
                elif ( (LTSRegionLocal[iCell] % 2 == 0) 
                        and (LTSRegionLocal[iCell] != 2) ): 
                    LTSRegionLocal[iCell] = 4
        
        # this means we are using LTS3
        elif (nLTSHalos == 3): 
            for iCell in range(0, nCells) :
                if( (LTSRegionLocal[iCell] % 2 == 1) 
                        and (LTSRegionLocal[iCell] != 1) ):
                    # if we are here it could be either 
                    # interface 1 or those two layers of 
                    # fine we need for the third order
                    if ( LTSRegionLocal[iCell] 
                            == 2 * nLTSHalosCopy + 1 ):
                        LTSRegionLocal[iCell] = 7
                    elif ( LTSRegionLocal[iCell] 
                            == 2 * nLTSHalosCopy - 1 ):
                        LTSRegionLocal[iCell] = 5 
                    else:
                        LTSRegionLocal[iCell] = 3 
                
                elif ( (LTSRegionLocal[iCell] % 2 == 0) 
                        and (LTSRegionLocal[iCell] != 2) ):
                    # if we are here it could be either interface 2 
                    # or those two layers of coarse
                    if (LTSRegionLocal[iCell] == 2 * nLTSHalosCopy + 2):
                        LTSRegionLocal[iCell] = 8 
                    elif (LTSRegionLocal[iCell] == 2 * nLTSHalosCopy):
                        LTSRegionLocal[iCell] = 6 
                    else :
                        LTSRegionLocal[iCell] = 4 

    fineCells = 0
    coarseCells = 0
    
    newf = ""
    with open(graph_info,'r') as f:
        lines = f.readlines()
        newf += lines[0].strip() + " 010 3 \n"
        for iCell in range(1,len(lines)):
            if ( LTSRegionLocal[iCell-1] == 1 
                    or LTSRegionLocal[iCell-1] == 5 
                    or LTSRegionLocal[iCell-1] == 7 ):  # fine 
                newf+= "0 1 0 " + lines[iCell].strip() + "\n"
                fineCells = fineCells + 1
            elif ( LTSRegionLocal[iCell-1] == 2 
                    or LTSRegionLocal[iCell-1] == 6 
                    or LTSRegionLocal[iCell-1] == 8 ):  # coarse
                newf+= "1 0 0 " + lines[iCell].strip() + "\n"
                coarseCells = coarseCells + 1
            else:  # interface 1 and 2
                newf+= "0 0 1 " + lines[iCell].strip() + "\n"
                coarseCells = coarseCells + 1 
        
    with open(graph_info, 'w') as f:
        f.write(newf)

    print(fineCells)
    print(coarseCells)
    ratio = max(areaCell) / min(areaCell)
    print(ratio)

# END of main()


if __name__ == '__main__':
    # If called as a primary module, run main
    main(args.base_mesh, args.graph_info, args.coarse_region_dist, args.lts2)

