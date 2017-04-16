#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, time, sys, shutil, copy
import numpy as np

import operator

from .. import io  as spaceio
from interface import SYM       as SPACE_SYM
from interface import VOL       as SPACE_VOL
from interface import BLG       as SPACE_BLG

# -------------------------------------------------------------------
#  Fluid Mesh Simulation
# -------------------------------------------------------------------

def fluid_mesh ( config ): 

    # local copy
    konfig = copy.deepcopy(config)

    # Volume Mesh
    
    konfig.SYMMETRY_FILENAME = 'symmetry.mesh'
    konfig.BOUNDARY_FILENAME = 'boundary.mesh'
    konfig.FLUID_SURFACE_CUR = konfig.FLUID_SURFACE

    SPACE_SYM(konfig)

    for index in range(3):

        SPACE_VOL(konfig)
        SPACE_BLG(konfig)

        if os.path.exists(konfig.FLUID_VOLUME + '.meshb'): break

    os.system('meshutils -O 3 -in ' + konfig.FLUID_VOLUME + '.meshb -out ' + konfig.FLUID_VOLUME + ' > log_meshutil.out')

    # Back Mesh

    konfig.SYMMETRY_FILENAME = 'symmetry_back.mesh'
    konfig.BOUNDARY_FILENAME = 'boundary_back.mesh'
    konfig.FLUID_SURFACE_CUR = konfig.FLUID_SURFACE + '_back'
    SPACE_SYM(konfig)

    # info out
    info = spaceio.State()
    info.FILES.FLUID_VOLUME_MESH = konfig.FLUID_VOLUME + '.su2'
    info.FILES.BOUNDARY_BACK_MESH = 'boundary_back.mesh'

    return info

