#!/usr/bin/env python

import argparse as ap
import os
import subprocess as sp
import xml.etree.ElementTree as et
from build_base_mesh_test5 import main as build_base_mesh
from build_weighted_graph_info_with_lts_regions_test5 import main as build_graph
from build_graph_info_part_for_multi_block_run import main as partition_graph


def main(outDir, modelRepo, baseMesh, graphInfo, model, numProcs, multiBlocks,
         numInterface, coarseRegionDist, disableOutput, doLTS2):
    
    startingDir = os.getcwd()
    
    if outDir[-1] != '/':
        outDir += '/'

    if modelRepo[-1] != '/':
        modelRepo += '/'

    if not os.path.exists(outDir):
        os.makedirs(outDir)


    # -- BEGIN --
    # mesh generation and MPI partitioning
    
    os.chdir(outDir)

    print('\n\n\n--- Building base mesh...\n\n\n')
    build_base_mesh(baseMesh, True)
    print('\n\n\n--- Done\n\n\n')

    print('\n\n\n--- Converting ' + baseMesh + ' to a MPAS mesh...\n\n\n')
    shCommand = 'MpasMeshConverter.x ' + baseMesh
    sp.call(shCommand.split())

    shCommand = 'mv graph.info ' + graphInfo
    sp.call(shCommand.split())
    print('\n\n\n--- Done\n\n\n')

    print('\n\n\n--- Weighting ' + graphInfo + ' for LTS regions...\n\n\n')
    build_graph(baseMesh, graphInfo, numInterface, coarseRegionDist, doLTS2)
    print('\n\n\n--- Done\n\n\n')

    print('\n\n\n--- Partitioning cells across MPI blocks with gpmetis...\n\n\n')
    if multiBlocks:
        numBlocks = 3 * numProcs
    else:
        numBlocks = numProcs

    shCommand = 'gpmetis ' + graphInfo + ' ' + str(numBlocks)
    sp.call(shCommand.split())
    print('\n\n\n--- Done\n\n\n')

    print('\n\n\n--- Resorting cells so that each MPI block only has one \n \
          type of cell--fine, interface, or coarse...\n\n\n')
    partition_graph(baseMesh, graphInfo, numBlocks)
    print('\n\n\n--- Done\n\n\n')

    # mesh generation and MPI partitioning
    # -- END --


    # -- BEGIN --
    # model config and build

    os.chdir(modelRepo)

    if doLTS2:
        nLTSHalos = 1
        LTSX = 'LTS2'
    else:
        nLTSHalos = 3
        LTSX = 'LTS3'

    if numInterface != 1:
        nLTSHalosCopy = (nLTSHalos - 1) + numInterface
        moreCellsOnInterface = 1
    else:
        nLTSHalosCopy = nLTSHalos
        moreCellsOnInterface = 0

    if disableOutput:
        modelOutput = 'none'
    else:
        modelOutput = 'output'


    print('\n\n\n--- Editing Registry.xml...\n\n\n')
    registryXML = modelRepo + 'src/core_sw/Registry.xml'

    registryTree = et.parse(registryXML)
    registryRoot = registryTree.getroot()

    registryDims = registryRoot.find('dims').findall('dim')
    for dim in registryDims:
        name = dim.attrib['name']
        if name == 'nLTSHalos':
            dim.set('definition', str(nLTSHalos))
        elif name == 'moreCellsOnInterface':
            dim.set('definition', str(moreCellsOnInterface))
        elif name == 'nLTSHalosPExtra':
            dim.set('definition', str(nLTSHalosCopy))
        elif name == 'nVertLevels':
            dim.set('definition', str(100))  # TODO

    registryVarStructs = registryRoot.findall('var_struct')
    for vs in registryVarStructs:
        name = vs.attrib['name']
        if name == 'LTS':
            variables = vs.findall('var')
            for var in variables:
                name = var.attrib['name']
                if name == 'coarseRegionDist':
                    var.set('default_value', str(coarseRegionDist))

    registryTree.write(registryXML) 
    print('\n\n\n--- Done\n\n\n')

    print('\n\n\n--- Building model...')
    shCommand = 'make clean CORE=sw'
    sp.call(shCommand.split())

    shCommand = 'make gfortran CORE=sw DEBUG=false'
    sp.call(shCommand.split())

    shCommand = 'cp sw_model ' + outDir + model
    sp.call(shCommand.split())
    print('\n\n\n--- Done\n\n\n')

    print('\n\n\n--- Configuring namelist.sw and streams.sw...\n\n\n')
    shCommand = 'cp namelist.sw streams.sw ' + outDir
    sp.call(shCommand.split())

    os.chdir(outDir)

    # config namelist
    newTxt = ''
    with open('namelist.sw', 'r') as namelistFile:
        namelistTxt = namelistFile.read()

        for line in namelistTxt.split('\n'):
            words = line.split()
            if 'config_time_integration' in words:
                words[-1] = "'" + LTSX + "'"
                newTxt += '    ' + ' '.join(words) + '\n'
            elif 'config_dt' in words:
                words[-1] = str(200)  # TODO
                newTxt += '    ' + ' '.join(words) + '\n'
            elif 'config_run_duration' in words:
                words[-1] = "'" + '00:30:00' + "'"  # TODO
                newTxt += '    ' + ' '.join(words) + '\n'
            elif 'config_block_decomp_file_prefix' in words:
                words[-1] = "'" + graphInfo + '.part.' + "'"
                newTxt += '    ' + ' '.join(words) + '\n'
            elif 'config_number_of_blocks' in words:
                words[-1] = str(numBlocks)
                newTxt += '    ' + ' '.join(words) + '\n'
            elif 'config_proc_decomp_file_prefix' in words:
                words[-1] = "'" + graphInfo + '.part.' + "'"
                newTxt += '    ' + ' '.join(words) + '\n'
            elif 'config_use_local_time_stepping' in words:
                words[-1] = 'true'
                newTxt += '    ' + ' '.join(words) + '\n'
            elif 'config_dt_scaling_LTS' in words:
                words[-1] = str(25)  # TODO this is our M
                newTxt += '    ' + ' '.join(words) + '\n'
            else:
                newTxt += line + '\n'

    with open('namelist.sw', 'w') as namelistFile:
        namelistFile.write(newTxt)

    # config streams.sw
    streamsXML = 'streams.sw'

    streamsTree = et.parse(streamsXML)
    streamsRoot = streamsTree.getroot()

    for element in streamsRoot:
        name = element.attrib['name']
        if name == 'input':
            element.set('filename_template', baseMesh)
        elif name == 'output':
            element.set('type', modelOutput)
            element.set('clobber_mode', 'truncate')
            element.set('output_interval', '00:30:00')  # TODO
    
    streamsTree.write(streamsXML) 
    print('\n\n\n--- Done\n\n\n')

    # model config and build
    # -- END --


    # Write chosen parameters to a text file for later reference
    paraTxt = 'Parameter list for test case 5 from Williamson et al. for LTS.'
    paraTxt += '\n\n'

    paraTxt += 'outDir = ' + outDir + '\n'
    paraTxt += 'modelRepo = ' + modelRepo + '\n'
    paraTxt += 'baseMesh = ' + baseMesh + '\n'
    paraTxt += 'graphInfo = ' + graphInfo + '\n'
    paraTxt += 'model = ' + model + '\n'
    paraTxt += 'numProcs = ' + str(numProcs) + '\n'
    paraTxt += 'multiBlocks = ' + str(multiBlocks) + '\n'
    paraTxt += 'numBlocks = ' + str(numBlocks) + '\n'
    paraTxt += 'numInterface = ' + str(numInterface) + '\n'
    paraTxt += 'coarseRegionDist = ' + str(0.55) + '\n'
    paraTxt += 'disableOutput = ' + str(disableOutput) + '\n'
    paraTxt += 'doLTS2 = ' + str(doLTS2) + '\n'

    paraTxt += '\nThis test case should be run from ' + outDir + ' with:\n'
    paraTxt += ('    mpirun -n ' + str(numProcs) + ' ' + str(model) + ' '
                + 'namelist.sw streams.sw\n')

    with open('parameterList.txt', 'w') as f:
        f.write(paraTxt)


# END main


if __name__ == '__main__':
    # If called as the primary module, run main

    parser = ap.ArgumentParser(description='Python script to generate the base \
                               mesh, graph.info file, block partition, and \
                               executable for running a variation of test case \
                               5 from Williamson et al.')

    parser.add_argument('-o', '--out-dir', dest='out_dir', required=True,
                        type=str,
                        help='Directory to store all output files in.')

    parser.add_argument('-r', '--model-repo', dest='model_repo', required=True,
                        type=str,
                        help='Path to the `MPAS-Model` repo to build the \
                        model from.')

    parser.add_argument('-b', '--base-mesh', dest='base_mesh',
                        default='base_mesh.nc', type=str,
                        help='Filename for the base mesh. Default is \
                        `base_mesh.nc`.')

    parser.add_argument('-g', '--graph-info', dest='graph_info',
                        default='graph.info', type=str,
                        help='Filename for graph.info file. Default is \
                        `graph.info`.')

    parser.add_argument('-m', '--model', default='sw_model', type=str,
                        help='Filename for model executable. Default is \
                        `sw_model`.')

    parser.add_argument('-n', '--num-procs', dest='num_procs',
                        default=128, type=int,
                        help='Number of MPI processors to config for. Default \
                        is 128.')

    parser.add_argument('-k', '--multi-blocks-per-proc', dest='multi_blocks',
                        action='store_true',
                        help='Set to partition the graph so that each MPI \
                        processor has three blocks--one for each type of LTS \
                        region i.e. NUM_BLOCKS = 3 * NUM_PROCS.')

    parser.add_argument('-l', '--num-interface-layers', 
                        dest='num_interface', default=1, type=int,
                        help='Number of interface layers to add for load \
                        balancing. Default is 1. NOTE: LTS3 requires that 2 is \
                        added to this number--THIS IS HANDELED AUTOMATICALLY. \
                        For example, running LTS3 and setting this flag to 38 \
                        results in 40 total interface layers for each region.')

    parser.add_argument('-c', '--coarse-region-dist', 
                        dest='coarse_region_dist', default=0.55, type=float,
                        help='Cells more than this distance way from the \
                        mountain at the north pole will be part of the coarse \
                        region. Default is 0.55.')

    parser.add_argument('-d', '--disable-output', dest='disable_output',
                        action='store_true',
                        help='Configure `streams.sw` to forgo producing \
                        visualization output.')

    parser.add_argument('--lts2', dest='do_lts2', action='store_true',
                        help='Prepare the mesh for LTS2 rather than LTS3 which \
                        is the default.')

    args = parser.parse_args()


    main(args.out_dir, args.model_repo, args.base_mesh, args.graph_info, 
         args.model, args.num_procs, args.multi_blocks, 
         args.num_interface, args.coarse_region_dist, 
         args.disable_output, args.do_lts2)

