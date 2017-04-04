#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy

from .. import io  as spaceio
from merge     import merge     as spacemerge
from interface import CFD       as SPACE_CFD

# ----------------------------------------------------------------------
#  Direct Simulation
# ----------------------------------------------------------------------

def direct ( config ): 
    
    # local copy
    konfig = copy.deepcopy(config)

    # setup direct problem
    konfig['MATH_PROBLEM']  = 'DIRECT'
    konfig['CONV_FILENAME'] = konfig['CONV_FILENAME'] + '_direct'    

    # Run Solution
    SPACE_CFD(konfig)

    # merge
    konfig['SOLUTION_FLOW_FILENAME'] = konfig['RESTART_FLOW_FILENAME']
    spacemerge(konfig)
    
    # filenames
    plot_format      = konfig['OUTPUT_FORMAT']
    plot_extension   = spaceio.get_extension(plot_format)
    history_filename = konfig['CONV_FILENAME'] + plot_extension

    # get history and objectives
    #history      = spaceio.read_history( history_filename )
    aerodynamics = spaceio.read_aerodynamics( history_filename )

    # info out
    info = spaceio.State()
    info.FUNCTIONS.update( aerodynamics )
    info.FILES.DIRECT = konfig['RESTART_FLOW_FILENAME']
    #info.HISTORY.DIRECT = history

    info.FILES.FLUID_SURFACE_FLOW = konfig.SURFACE_FLOW_FILENAME + '.dat'
    
    return info
