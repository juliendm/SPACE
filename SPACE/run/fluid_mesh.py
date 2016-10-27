#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, time, sys, shutil, copy
import numpy as np

import operator

from .. import io  as spaceio
from interface import GHS       as SPACE_GHS

# -------------------------------------------------------------------
#  Fluid Mesh Simulation
# -------------------------------------------------------------------

def fluid_mesh ( config ): 

    # local copy
    konfig = copy.deepcopy(config)

    # Read Correspondance
    corresp = np.loadtxt(konfig.CORRESPONDANCE_FILENAME, dtype='int')

    # Read Fluid Surface and Fluid Boundary
    node_fluid_surface = []
    fluid_surface = open(konfig.FLUID_SURFACE + '.mesh')
    line = fluid_surface.readline()
    while line[0] != 'V':
        line = fluid_surface.readline()
    nPoint = int(fluid_surface.readline().split()[0])
    fluid_surface.readline()
    for i in range(nPoint):
        data = fluid_surface.readline().split()
        node_fluid_surface.append([float(data[0]),float(data[1]),float(data[2])])
    fluid_surface.close()

    node_fluid_boundary = []
    fluid_boundary = open(konfig.FLUID_BOUNDARY_FILENAME)
    line = fluid_boundary.readline()
    while line[0] != 'V':
        line = fluid_boundary.readline()
    nPoint = int(fluid_boundary.readline().split()[0])
    fluid_boundary.readline()
    for i in range(nPoint):
        data = fluid_boundary.readline().split()
        node_fluid_boundary.append([float(data[0]),float(data[1]),float(data[2])])
    remaining = ''
    line = fluid_boundary.readline()
    while line:
        remaining += line
        line = fluid_boundary.readline()
    fluid_boundary.close()
    nDim = 3
    for index in range(len(corresp)):
        for iDim in range(nDim):
            node_fluid_boundary[corresp[index]-1][iDim] = node_fluid_surface[index][iDim]

    # Write FLuid Boundary Updated
    fluid_boundary_updated_filename = konfig.FLUID_BOUNDARY_FILENAME.split('.')
    fluid_boundary_updated_filename = fluid_boundary_updated_filename[0] + '_updated.' + fluid_boundary_updated_filename[1]
    fluid_boundary_updated = open(fluid_boundary_updated_filename,'w')
    nPoint = len(node_fluid_boundary)
    fluid_boundary_updated.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nVertices\n' + str(nPoint) + '\n\n')
    for iPoint in range(nPoint):
        fluid_boundary_updated.write("%f %f %f 0\n" % (node_fluid_boundary[iPoint][0], node_fluid_boundary[iPoint][1], node_fluid_boundary[iPoint][2]))
    fluid_boundary_updated.write(remaining) 
    fluid_boundary_updated.close()

    # Run Solution
    SPACE_GHS(konfig)

    # info out
    info = spaceio.State()
    info.FILES.FLUID_BOUNDARY_UPDATED = fluid_boundary_updated_filename
    info.FILES.FLUID_VOLUME_MESH = konfig.FLUID_VOLUME + '.meshb'
    info.FILES.FLUID_VOLUME_SU2 = konfig.FLUID_VOLUME + '.su2'

    return info

#: def boundary()
