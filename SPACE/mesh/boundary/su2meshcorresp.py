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
    parser.add_option("-v", "--volume", dest="volume",
                      help="read config from VOLUME", metavar="VOLUME")
    parser.add_option("-s", "--surface", dest="surface",
                      help="read config from SURFACE", metavar="SURFACE")

    (options, args)=parser.parse_args()

    if options.volume == None:
        raise Exception("No volume file provided. Use -v flag")
    if options.surface == None:
        raise Exception("No surface file provided. Use -s flag")

    corresp_file = 'correspondance.dat'

    # READ

    inp = open(options.volume)

    data = inp.readline().split('=')
    while 'NPOIN' not in data[0]:
        data = inp.readline().split('=')
    nPoint_volume = int(data[1])
    #print nPoint_volume
    node_volume = []
    for index in range(nPoint_volume):
        data = inp.readline().split()
        node_volume.append([round(float(data[0].strip()),5), round(float(data[1].strip()),5), round(float(data[2].strip()),5), int(data[3].strip())])

    data = inp.readline()
    while 'BASELINE' not in data:
        data = inp.readline()
    data = inp.readline().split('=')
    nElem_baseline = int(data[1])
    #print nElem_baseline
    elem_baseline = []
    for index in range(nElem_baseline):
        data = inp.readline().split()
        elem_baseline.append([int(data[1].strip()), int(data[2].strip()), int(data[3].strip())])

    inp.close()

    inp = open(options.surface)
    line = inp.readline()
    while line[0] != 'V':
        line = inp.readline()
    line = inp.readline()
    nPoint_surface = int(line.split()[0])
    print nPoint_surface
    #line = inp.readline()
    node_surface = []
    for index in range(nPoint_surface):
        data = inp.readline().split()
        node_surface.append([round(float(data[0].strip()),5), round(float(data[1].strip()),5), round(float(data[2].strip()),5), index])
    inp.close()

    # CORRESPONDANCE

    ids_volume_oml = []
    for iElem in range(nElem_baseline):
        elem_i = elem_baseline[iElem]
        ids_volume_oml.append(elem_i[0])
        ids_volume_oml.append(elem_i[1])
        ids_volume_oml.append(elem_i[2])
    ids_volume_oml = list(set(ids_volume_oml))

    print len(ids_volume_oml)

    node_volume_oml = []
    for id_volume_oml in ids_volume_oml:
        node_volume_oml.append(node_volume[id_volume_oml])
    
    node_surface.sort(key = operator.itemgetter(0, 1, 2))
    node_volume_oml.sort(key = operator.itemgetter(0, 1, 2))

    correspondance = []
    for index in range(len(node_volume_oml)):
        if (abs(node_surface[index][0]-node_volume_oml[index][0]) > 1e-3 or abs(node_surface[index][1]-node_volume_oml[index][1]) > 1e-3 or abs(node_surface[index][2]-node_volume_oml[index][2]) > 1e-3):
            print node_surface[index], node_volume_oml[index]
        else:
            correspondance.append([node_surface[index][3],node_volume_oml[index][3]])
    correspondance.sort(key = operator.itemgetter(0))

    # WRITE

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
