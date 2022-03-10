#!/usr/bin/env python

import os
import subprocess as sp
import shutil
import numpy as np
import xarray as xr
import argparse
import math


def main(mesh, graph_info, num_blocks):

    shCommand = 'gpmetis ' + graph_info + ' ' + str(num_blocks)
    sp.call(shCommand.split())
    
    ds = xr.open_dataset(mesh)

    nCells = ds['nCells'].size
    LTSRegionLocal = np.zeros([nCells])

    with open(graph_info, 'r') as gi:
        lines = gi.readlines()
        for iCell in range(1, len(lines)):
            if lines[iCell][0:5] == '0 1 0':
                LTSRegionLocal[iCell-1] = int(1)  # fine
            elif lines[iCell][0:5] == '1 0 0':
                LTSRegionLocal[iCell-1] = int(2)  # coarse
            else:
                LTSRegionLocal[iCell-1] = int(3)  # interface

    numBlocks = num_blocks  # usually 3 * NUM_PROCS

    newf=""
    with open(graph_info + '.part.' + str(int(numBlocks)), 'r') as f:
        lines = f.readlines()
        procFoundCell = [0] * nCells

        for iCell in range(0, len(lines)):
            procCell = int(numBlocks / 3 + 1)  # just some initialization
            for iProc in range(0, int(numBlocks/3)) :
                # NOTE: the if statement below is like that because the 
                # blocks assigned to a given processor go like 
                # iProc*2 + iProc + k, where k=0,1,2
                for k in range(0, 3) :
                    if int(lines[iCell][0:10]) == (iProc*2 + iProc + k) :
                        procCell = int(iProc)
                        procFoundCell[iCell] = 1
                        break
                if procFoundCell[iCell] == 1 :
                    break

            blockID = str(numBlocks + 1)  # just some initialization
            if LTSRegionLocal[iCell] == 1:
                # blocks 0,3,6,9 etc are FINE
                tmp = procCell * 2 + procCell
                blockID = str(tmp)
            elif LTSRegionLocal[iCell] == 2:   
                # blocks 1,4,7,10 etc are COARSE
                tmp = procCell * 2 + procCell + 1
                blockID = str(tmp)
            else:
                # blocks 2,5,8,11 etc are INTERFACE
                tmp = procCell * 2 + procCell + 2
                #tmp = 2 #this gives all the interface cells to block 2 
                # (case C in the paper)
                blockID = str(tmp)
                #if procCell % 2 == 0 :
                #   tmp = 2
                #   blockID = str(tmp)
                #else :
                #   tmp = 5
                #   blockID = str(tmp)
            
            newf+= blockID + "\n"
       
        print('If all procCells have been found, these two numbers are equal:',
              sum(procFoundCell[:]), nCells)
    
    with open(graph_info + '.part.lts.' + str(int(numBlocks)), 'w') as f:
        f.write(newf)


if __name__ == '__main__':
    # If called as a primary module, run main

    parser = argparse.ArgumentParser(description='Python script to repartition \
                                     an existing graph.info.part. file so that \
                                     each block contains only one type of cell \
                                     from fine, coarse, and interface. Before \
                                     running this script, one must generate a \
                                     graph.info.part.NUM_BLOCKS file from an \
                                     existing graph.info file by running \
                                     `gpmetis graph.info NUM_BLOCKS`.')

    parser.add_argument('-m', '--mesh', dest='mesh',
                        default='input.nc',
                        help='File containing the mesh. Default is \
                        `input.nc`.')

    parser.add_argument('-g', '--graph-info', dest='graph_info',
                        default='graph.info',
                        help='graph.info file corresponding to the mesh. \
                        Default is `graph.info.lts`.')

    parser.add_argument('-b', '--num-blocks', dest='num_blocks',
                        default=12,
                        help='Number of blocks to partition for. Usually this \
                        is 3 * NUM_PROCS, where NUM_PROCS is the number of MPI \
                        processors that will be used. Default is 12.')

    args = parser.parse_args()


    main(args.mesh, args.graph_info, args.num_blocks)

