#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def surf2pressuresol( config ):

    surface_flow = open(config.FLUID_SURFACE_FLOW + '.dat')
    sol = open(config.FLUID_SURFACE + '.sol','w')

    # Read Write

    line = surface_flow.readline()
    line = surface_flow.readline()

    data = surface_flow.readline().split(',')
    nVertex = int(data[0].split('=')[-1])

    sol.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nSolAtVertices\n' + str(nVertex) + '\n1 1\n')
    for iVertex in range(nVertex):
       data = surface_flow.readline().split()
       sol.write(data[10] + "\n")

    sol.write("\nEnd\n")

    surface_flow.close()
    sol.close()

#: def surf2pressuresol()
