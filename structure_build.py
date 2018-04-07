#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, glob, re
# from .. import io   as spaceio
# from .  import func as spacefunc
# from ..io import redirect_folder, save_data

import subprocess
import numpy

SPACE_RUN = os.environ['SPACE_RUN']

from optparse import OptionParser

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE import io   as spaceio
from SPACE.eval import func as spacefunc
from SPACE.eval import design as spacedesign
from SPACE.io import redirect_folder, save_data

from SPACE.util import LHC_unif, DesignVariables
from SPACE.surfpack import Surfpack

# ----------------------------------------------------------------------
#  Main Function Interface
# ----------------------------------------------------------------------

def build(regime = 'ON'):

    regime = 'OFF'

    models_folder = 'MODELS'
    build_points_folder = 'BUILD_POINTS'
    desvar = DesignVariables()

    if regime == 'ON':
        ndim_struct = desvar.ndim_struct_on
    elif regime == 'OFF':
        ndim_struct = desvar.ndim_struct_off

    n_models = 132
    for index in range(0,n_models):

        mass_model = Surfpack('MASS_%05d' % (index+1), ndim_struct)
        mass_model.load_data(os.path.join(build_points_folder,'enriched_points_mass_%05d.dat' % (index+1)))
        mass_model.build('kriging')
        mass_model.save_model(os.path.join(models_folder,'model_mass_%05d.sps' % (index+1)))

    structure_mass_model = Surfpack('STRUCTURE_MASS', ndim_struct)
    structure_mass_model.load_data(os.path.join(build_points_folder,'enriched_points_STRUCTURE_MASS.dat'))
    structure_mass_model.build('kriging')
    structure_mass_model.save_model(os.path.join(models_folder,'model_structure_mass.sps'))

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-r", "--regime", dest="regime", default="ON",
                      help="regime", metavar="REGIME")


    (options, args)=parser.parse_args()

    build( options.regime )


