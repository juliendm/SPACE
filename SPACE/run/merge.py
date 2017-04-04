#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy

from .. import io  as spaceio
from interface import SOL as SPACE_SOL

# ----------------------------------------------------------------------
#  Merge Mesh
# ----------------------------------------------------------------------

def merge( config ):
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # check if needed
    partitions = konfig['NUMBER_PART']
    if partitions <= 1:
        return spaceio.State()
    
    # # MERGING # #
    merge_solution(konfig)
        
    # info out (empty)
    info = spaceio.State()
    
    return info

#: merge

def merge_solution( config ):
    """ SU2.io.merge.merge_solution(config)
        general volume surface merging with SPACE_SOL
    """
    
    SPACE_SOL( config )
    
    return

#: merge_solution( config )
