#!/usr/bin/env python

import os
import subprocess as sp
import shutil
import numpy as np
import xarray as xr
import argparse as ap
import math
import time
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


def main(base_mesh, graph_info, num_interface):

    # read in mesh data
    ds = xr.open_dataset(base_mesh)
    nCells = ds['nCells'].size
    nEdges = ds['nEdges'].size
    areaCell = ds['areaCell'].values
    nEdgesOnCell = ds['nEdgesOnCell'].values
    cellsOnEdge = ds['cellsOnEdge'].values
    edgesOnCell = ds['edgesOnCell'].values
    latCell = ds['latCell']
    lonCell = ds['lonCell']
    

    #####
    # label cells for LTS
    #####

    # [lat, lon] points defining the fine region
    fineRegionPts = np.array([[0.481, -1.737 + 2 * math.pi],
                              [0.311, -1.701 + 2 * math.pi],
                              [0.234, -1.508 + 2 * math.pi],
                              [0.148, -1.430 + 2 * math.pi],
                              [0.151, -1.397 + 2 * math.pi],
                              [0.163, -1.383 + 2 * math.pi],
                              [0.120, -1.320 + 2 * math.pi],
                              [0.077, -0.921 + 2 * math.pi],
                              [0.199, -0.784 + 2 * math.pi],
                              [0.496, -0.750 + 2 * math.pi],
                              [0.734, -0.793 + 2 * math.pi],
                              [0.826, -0.934 + 2 * math.pi],
                              [0.871, -1.001 + 2 * math.pi],
                              [0.897, -0.980 + 2 * math.pi],
                              [0.914, -1.012 + 2 * math.pi],
                              [0.850, -1.308 + 2 * math.pi],
                              [0.743, -1.293 + 2 * math.pi],
                              [0.638, -1.781 + 2 * math.pi],
                              [0.481, -1.737 + 2 * math.pi]])
    fineRegion = Polygon(fineRegionPts)

    # start by assuming all cells set to coarse
    LTSRegion = [2] * nCells

    # check each cell, if in the fine region, label as fine
    print('Labeling fine cells...')
    for iCell in range(0, nCells):
        cellPt = Point(latCell[iCell], lonCell[iCell])
        if fineRegion.contains(cellPt):
            LTSRegion[iCell] = 1
        # END if
    # END for

    # first layer of cells with label 5
    changedCells = [[], []]
    for iEdge in range(0, nEdges):
        cell1 = cellsOnEdge[iEdge, 0] - 1
        cell2 = cellsOnEdge[iEdge, 1] - 1
        
        if (cell1 != -1 and cell2 != -1):
            if (LTSRegion[cell1] == 1 and
                LTSRegion[cell2] == 2):
                
                LTSRegion[cell2] = 5
                changedCells[0].append(cell2)

            elif (LTSRegion[cell1] == 2 and
                  LTSRegion[cell2] == 1):
                
                LTSRegion[cell1] = 5
                changedCells[0].append(cell1)
            # END if
        # END if
    # END for

    # second and third layer of cells with label 5
    # only looping over cells changed during loop for previous layer
    # at the end of this loop, changedCells[0] will have the list of cells
    # sharing edegs with the coarse cells
    print('Labeling interface-adjacent fine cells...')
    for i in range(0, 2):  # this loop creates 2 layers
        changedCells[(i+1)%2] = []

        for iCell in changedCells[i%2]:
            edges = edgesOnCell[iCell]
            for iEdge in edges:
                if iEdge != 0:
                    cell1 = cellsOnEdge[iEdge-1, 0] - 1
                    cell2 = cellsOnEdge[iEdge-1, 1] - 1

                    if (cell1 != -1 and cell2 != -1):
                        if (LTSRegion[cell1] == 5 and
                            LTSRegion[cell2] == 2):
                            
                            LTSRegion[cell2] = 5
                            changedCells[(i+1)%2].append(cell2)

                        elif (LTSRegion[cell1] == 2 and
                              LTSRegion[cell2] == 5):
                            
                            LTSRegion[cell1] = 5
                            changedCells[(i+1)%2].append(cell1)
                        # END if
                    # END if
                # END if
            # END for
        # END for
    # END for

    # n layers of interface region with label 4
    print('Labeling interface cells...')
    for i in range(0, num_interface):
        changedCells[(i+1)%2] = []
        
        for iCell in changedCells[i%2]:
            edges = edgesOnCell[iCell]
            for iEdge in edges:
                if iEdge != 0:
                    cell1 = cellsOnEdge[iEdge-1, 0] - 1
                    cell2 = cellsOnEdge[iEdge-1, 1] - 1
                    
                    if (cell1 != -1 and cell2 != -1):
                        # for the first layer, need to check neighbors are
                        # 5 and 2
                        # for further layers, need to check neighbors are
                        # 3 and 2
                        if (i == 0):
                            if (LTSRegion[cell1] == 5 and 
                                LTSRegion[cell2] == 2):
                                
                                LTSRegion[cell2] = 3
                                changedCells[(i+1)%2].append(cell2)

                            elif (LTSRegion[cell1] == 2 and
                                  LTSRegion[cell2] == 5):
                                
                                LTSRegion[cell1] = 3
                                changedCells[(i+1)%2].append(cell1)
                            # END if
                        else:
                            if (LTSRegion[cell1] == 3 and 
                                LTSRegion[cell2] == 2):
                                
                                LTSRegion[cell2] = 3
                                changedCells[(i+1)%2].append(cell2)

                            elif (LTSRegion[cell1] == 2 and
                                  LTSRegion[cell2] == 3):
                                
                                LTSRegion[cell1] = 3
                                changedCells[(i+1)%2].append(cell1)
                            # END if
                        # END if 
                    # END if
                # END if
            # END for
        # END for
    # END for

    changedCells[0] = changedCells[num_interface%2]
    
    # n layers of interface region with label 3
    for i in range(0, num_interface):
        changedCells[(i+1)%2] = []
        
        for iCell in changedCells[i%2]:
            edges = edgesOnCell[iCell]
            for iEdge in edges:
                if iEdge != 0:
                    cell1 = cellsOnEdge[iEdge-1, 0] - 1
                    cell2 = cellsOnEdge[iEdge-1, 1] - 1
                    
                    if (cell1 != -1 and cell2 != -1):
                        # for the first layer, need to check neighbors are
                        # 5 and 2
                        # for further layers, need to check neighbors are
                        # 3 and 2
                        if (i == 0):
                            if (LTSRegion[cell1] == 3 and 
                                LTSRegion[cell2] == 2):
                                
                                LTSRegion[cell2] = 4
                                changedCells[(i+1)%2].append(cell2)

                            elif (LTSRegion[cell1] == 2 and
                                  LTSRegion[cell2] == 3):
                                
                                LTSRegion[cell1] = 4
                                changedCells[(i+1)%2].append(cell1)
                            # END if
                        else:
                            if (LTSRegion[cell1] == 4 and 
                                LTSRegion[cell2] == 2):
                                
                                LTSRegion[cell2] = 4
                                changedCells[(i+1)%2].append(cell2)

                            elif (LTSRegion[cell1] == 2 and
                                  LTSRegion[cell2] == 4):
                                
                                LTSRegion[cell1] = 4
                                changedCells[(i+1)%2].append(cell1)
                            # END if
                        # END if 
                    # END if
                # END if
            # END for
        # END for
    # END for


    #####
    ### create lts_mesh.nc
    #####

    print('Creating lts_mesh.nc...')

    LTSDataArray = xr.DataArray(LTSRegion)
    LTSDict = {'LTSRegionLocal': (['nCells'], LTSDataArray.data)}
    LTSds = xr.Dataset(LTSDict)

    encodingDict = {'LTSRegionLocal': {'_FillValue': -1}}
    
    combined = xr.merge([ds, LTSds], combine_attrs='drop_conflicts')
    combined.to_netcdf(path='lts_mesh.nc', 
                       encoding=encodingDict,
                       format='NETCDF3_64BIT')

    shCommand = 'paraview_vtk_field_extractor.py -f lts_mesh.nc \
                 -v allOnCells,allOnEdges -d TWO= maxEdges=  maxEdges2= \
                 -o lts_mesh_vtk'
    sp.call(shCommand.split())


    #####
    # label cells in graph.info
    #####

    print('Weighting ' + graph_info + '...')

    fineCells = 0
    coarseCells = 0
    
    newf = ""
    with open(graph_info,'r') as f:
        lines = f.readlines()
        newf += lines[0].strip() + " 010 3 \n"
        for iCell in range(1,len(lines)):
            if (LTSRegion[iCell-1] == 1 or
                LTSRegion[iCell-1] == 5):  # fine 
                
                newf+= "0 1 0 " + lines[iCell].strip() + "\n"
                fineCells = fineCells + 1
            
            elif (LTSRegion[iCell-1] == 2):  # coarse
                newf+= "1 0 0 " + lines[iCell].strip() + "\n"
                coarseCells = coarseCells + 1
            
            else:  # interface 1 and 2
                newf+= "0 0 1 " + lines[iCell].strip() + "\n"
                coarseCells = coarseCells + 1
            # END if
        # END for
    # END with
        
    with open(graph_info, 'w') as f:
        f.write(newf)

    areaRatio = max(areaCell) / min(areaCell)
    numberRatio = coarseCells / fineCells

    print('Number of fine cells = ' + str(fineCells))
    print('Number of coarse cells = ' + str(coarseCells))
    print('Ratio of largest cell area to smallest cell area = ' 
          + str(areaRatio))
    print('Ratio of number of coarse cells to number of fine cells = '
          + str(numberRatio))

    return fineCells, coarseCells, areaRatio, numberRatio


# END of main()


if __name__ == '__main__':
    # If called as a primary module, run main

    parser = ap.ArgumentParser(description='Python script to label the cells \
                               from a given MPAS mesh for LTS.')

    parser.add_argument('-m', '--mesh', dest='mesh',
                        default='culled_mesh.nc',
                        help='File containing the mesh to label. Default is \
                        `culled_mesh.nc`.')

    parser.add_argument('-g', '--graph-info', dest='graph_info',
                        default='culled_graph.info',
                        help='The graph.info file corresponding to the \
                        mesh. Default is `culled_graph.info`.')

    parser.add_argument('-l', '--num-interface-layers', 
                        dest='num_interface', default=2, type=int,
                        help='Number of extra interface layers to add for load \
                        balancing. Default is 2 (the minimum for an ocean \
                        simulation).')

    args = parser.parse_args()


    main(args.mesh, args.graph_info, args.num_interface)

