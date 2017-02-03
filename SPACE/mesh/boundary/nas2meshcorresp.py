#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np
from optparse import OptionParser

import operator

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-b", "--boundary", dest="boundary",
                      help="read config from BOUNDARY", metavar="BOUNDARY")
    parser.add_option("-s", "--surface", dest="surface",
                      help="read config from SURFACE", metavar="SURFACE")

    (options, args)=parser.parse_args()

    if options.boundary == None:
        raise Exception("No config file provided. Use -b flag")
    if options.surface == None:
        raise Exception("No config file provided. Use -s flag")

    out_file = options.boundary.split('.')[0] + '.mesh'
    corresp_file = 'correspondance.dat'

    oml_ref = 1
    farfield_ref = 2
    symmetry_ref = 3

    # READ

    inp = open(options.boundary)
    data = inp.readline().split(',')
    while 'GRID' not in data[0]:
        data = inp.readline().split(',')
        if len(data) == 0: break
    node_boundary = []
    while 'GRID' in data[0]:
        node_boundary.append([float(data[3].strip()), float(data[4].strip()), float(data[5].strip())])
        data = inp.readline().split(',')
        if len(data) == 0: break
    nPoint_boundary = len(node_boundary)
    while 'CTRIA3' not in data[0]:
        data = inp.readline().split(',')
        if len(data) == 0: break
    elem_boundary = []
    while 'CTRIA3' in data[0]:
        elem_boundary.append([int(data[3].strip()), int(data[4].strip()), int(data[5].strip()), int(data[2].strip())])
        data = inp.readline().split(',')
        if len(data) == 0: break
    nElem_boundary = len(elem_boundary)
    inp.close()

    inp = open(options.surface)
    data = inp.readline().split(',')
    while 'GRID' not in data[0]:
        data = inp.readline().split(',')
        if len(data) == 0: break
    node_surface = []
    index = 0
    while 'GRID' in data[0]:
        index += 1
        node_surface.append([round(float(data[3].strip()),5), round(float(data[4].strip()),5), round(float(data[5].strip()),5), index])
        data = inp.readline().split(',')
        if len(data) == 0: break
    nPoint_surface = len(node_surface)
    inp.close()

    # CORRESPONDANCE

    ids_boundary_oml = []
    for iElem in range(nElem_boundary):
        elem_i = elem_boundary[iElem]
        if (elem_i[3] == oml_ref):
            ids_boundary_oml.append(elem_i[0])
            ids_boundary_oml.append(elem_i[1])
            ids_boundary_oml.append(elem_i[2])
    ids_boundary_oml = list(set(ids_boundary_oml))

    node_boundary_oml = []
    for id_boundary_oml in ids_boundary_oml:
        node_boundary_oml.append([round(node_boundary[id_boundary_oml-1][0],5),round(node_boundary[id_boundary_oml-1][1],5),round(node_boundary[id_boundary_oml-1][2],5),id_boundary_oml])
    
    node_surface.sort(key = operator.itemgetter(0, 1, 2))
    node_boundary_oml.sort(key = operator.itemgetter(0, 1, 2))

    correspondance = []
    skip = 0
    for index in range(len(node_boundary_oml)):
        if (abs(node_surface[index-skip][0]-node_boundary_oml[index][0]) > 1e-3 or abs(node_surface[index-skip][1]-node_boundary_oml[index][1]) > 1e-3 or abs(node_surface[index-skip][2]-node_boundary_oml[index][2]) > 1e-3):
            #print node_surface[index], node_boundary_oml[index+skip]
            skip = skip+1
        else:
            correspondance.append([node_surface[index-skip][3],node_boundary_oml[index][3]])
    correspondance.sort(key = operator.itemgetter(0))

    # WRITE

    out = open(out_file,'w')
    out.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nVertices\n' + str(nPoint_boundary) + '\n\n')
    for iPoint in range(nPoint_boundary):
        out.write("%f %f %f 0\n" % (node_boundary[iPoint][0], node_boundary[iPoint][1], node_boundary[iPoint][2]))
    out.write('\nTriangles\n' + str(nElem_boundary) + '\n\n')
    for iElem in range(nElem_boundary):
        out.write("%d %d %d %d\n" % (elem_boundary[iElem][0], elem_boundary[iElem][1], elem_boundary[iElem][2], elem_boundary[iElem][3]))
    out.close()

    corresp = open(corresp_file,'w')
    for iPoint in range(len(correspondance)):
        corresp.write("%d\n" % (correspondance[iPoint][1]))
    corresp.close()

#: def main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
