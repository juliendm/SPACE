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

# -------------------------------------------------------------------
#  Fluid Mesh Simulation
# -------------------------------------------------------------------

def fluid_mesh ( config ): 

    # local copy
    konfig = copy.deepcopy(config)

    konfig.SYMMETRY_FILENAME = 'symmetry.mesh'
    konfig.BOUNDARY_FILENAME = 'boundary.mesh'
    SPACE_SYM(konfig)

    SPACE_VOL(konfig)

    # info out
    info = spaceio.State()
    info.FILES.FLUID_VOLUME_MESH = konfig.FLUID_VOLUME + '.meshb'
    info.FILES.FLUID_VOLUME_SU2 = konfig.FLUID_VOLUME + '.su2'

    return info

