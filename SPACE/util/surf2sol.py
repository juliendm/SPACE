#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def surf2sol( config ):

    surface_flow = open(config.FLUID_SURFACE_FLOW + '.dat')

    mesh = open(config.FLUID_SURFACE + '.mesh','w')    
    sol = open(config.FLUID_SURFACE + '.sol','w')

    # Read Write

    line = surface_flow.readline()
    line = surface_flow.readline()

    data = surface_flow.readline().split(',')
    nVertex = int(data[0].split('=')[-1])
    nElem = int(data[1].split('=')[-1])

    mesh.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nVertices\n' + str(nVertex) + '\n')
    sol.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nSolAtVertices\n' + str(nVertex) + '\n4 1 1 1 1\n')

    for iVertex in range(nVertex):
       data = surface_flow.readline().split()
       mesh.write(data[0] + " " + data[1] + " " + data[2] + " 0 \n")
       sol.write(data[12] + " " + data[15] + " " + data[16] + " " + data[17] + "\n")

    mesh.write('\nTriangles\n' + str(nElem) + '\n')

    for iElem in range(nElem):
       data = surface_flow.readline().split()
       mesh.write(data[0] + " " + data[1] + " " + data[2] + " 1 \n")

    mesh.write("\nEnd\n")
    sol.write("\nEnd\n")

    surface_flow.close()
    mesh.close()
    sol.close()

#: def surf2pressuresol()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':

    sys.path.append(os.environ['SPACE_RUN'])
    import SPACE
    from SPACE import io   as spaceio

    config = spaceio.Config('config_DSN.cfg')
    surf2pressuresol( config )
