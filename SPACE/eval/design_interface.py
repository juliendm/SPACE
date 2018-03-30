#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, glob, re
# from .. import io   as spaceio
# from .  import func as spacefunc
# from ..io import redirect_folder, save_data

import subprocess

SPACE_RUN = os.environ['SPACE_RUN']

# from optparse import OptionParser
sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE import io   as spaceio
from SPACE.eval import func as spacefunc
from SPACE.eval import design as spacedesign
from SPACE.io import redirect_folder, save_data

from SPACE.surfpack import Surfpack
from SPACE.util import LHC_unif, DesignVariables

# ----------------------------------------------------------------------
#  Main Function Interface
# ----------------------------------------------------------------------

def eval_design(func_name, config):
    """ interface for async jobs
    """

    design_filename = 'design.pkl'
    if os.path.exists(design_filename):
        design = spaceio.load_data(design_filename)
        design.config = config
        save_data(design.filename,design)
    else:
        design = init_design(config)
        save_data(design.filename,design)

    # ##################
    # config = SPACE.io.Config('config_DSN.cfg')
    # design.config = config
    # ##################







    load_models = False

    if load_models:

        desvar = DesignVariables()

        n_models = 2180

        models_folder = '../../../RESPONSE_SURFACE_DV_SUP/DESIGNS/MODELS'

        for index in range(n_models):
            cp_model = Surfpack('CP_%05d' % (index+1), desvar.ndim)
            cp_model.load_model(os.path.join(models_folder,'model_cp_%05d.sps' % (index+1)))

        for index in range(n_models):
            cfx_model = Surfpack('CFX_%05d' % (index+1), desvar.ndim)
            cfx_model.load_model(os.path.join(models_folder,'model_cfx_%05d.sps' % (index+1)))

        for index in range(n_models):
            cfy_model = Surfpack('CFY_%05d' % (index+1), desvar.ndim)
            cfy_model.load_model(os.path.join(models_folder,'model_cfy_%05d.sps' % (index+1)))

        for index in range(n_models):
            cfz_model = Surfpack('CFZ_%05d' % (index+1), desvar.ndim)
            cfz_model.load_model(os.path.join(models_folder,'model_cfz_%05d.sps' % (index+1)))




    if func_name == 'AERODYNAMICS':
#        design.func('GEOMETRY')   # this way, will save intermadiate design
#        design.func('FLUID_MESH') # this way, will save intermadiate design
        design.func(func_name)
    else:
        design.func(func_name)

def init_design(config): # At this point, it is already determined that the design don't already exists
    """ starts a new design
        works in design folder
    """
    
    konfig = copy.deepcopy(config)
    #ztate  = copy.deepcopy(self.state)

    # start new design
    design = spacedesign.Design(konfig)
    
    # update local state filenames ( ??? why not in Design() )
    for key in design.files:
        name = design.files[key]
        name = os.path.split(name)[-1]
        design.files[key] = name
       
    return design

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    
    # Command Line Options
    # parser=OptionParser()
    # parser.add_option("-f", "--func",       dest="func_name",
    #                   help="read function", metavar="FUNCTION_NAME")
                      
    # (options, args)=parser.parse_args()

    config = spaceio.Config('config_DSN.cfg')
    eval_design(sys.argv[1], config)



