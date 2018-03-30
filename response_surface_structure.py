#!/usr/bin/env python2.7

import time, os, gc, sys, shutil, copy, math
import numpy as np
import pwd
import pickle

from optparse import OptionParser

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE.surfpack import Surfpack
from SPACE.util import LHC_unif, DesignVariables

sys.path.append(os.environ['SU2_RUN'])
import SU2

def main():

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--project", dest="project_folder",
                      help="project folder", metavar="PROJECT_FOLDER")
    parser.add_option("-i", "--initiate", dest="initiate", default="False",
                      help="initiate", metavar="INITIATE")
    parser.add_option("-n", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")


    (options, args)=parser.parse_args()
    options.initiate = options.initiate.upper() == 'TRUE'
    options.partitions = int( options.partitions )

    response_surface( options.filename       ,
                      options.project_folder ,
                      options.initiate       ,
                      options.partitions  )

#: main()

def response_surface( filename          ,
                      project_folder    ,
                      initiate = False  ,
                      partitions  = 0  ):

    # Project

    if os.path.exists(project_folder):
        project = SPACE.io.load_data(os.path.join(project_folder,'project.pkl'))
        project.compile_designs()
        config = project.config
    else:
        config = SPACE.io.Config(filename)
        state  = SPACE.io.State()
        project = SPACE.project.Project(config, state, folder=project_folder)

    print '%d design(s) so far' % len(project.designs)

    konfig = copy.deepcopy(config)
    konfig.NUMBER_PART = 0

    # Design Variables

    desvar = DesignVariables()


    XB = desvar.XB_STRUCT

    nd = 300
    dvs_struct_filename = 'dvs_struct.dat'

    # X = LHC_unif(XB,nd)
    # np.savetxt(dvs_struct_filename,X)

    X = np.loadtxt(dvs_struct_filename)

    # # Load All Models

    # n_models = 2180

    # models_folder = 'RESPONSE_SURFACE_DV_SUP/DESIGNS/MODELS'

    # for index in range(n_models):
    #     cp_model = Surfpack('CP_%05d' % (index+1), desvar.ndim)
    #     cp_model.load_model(os.path.join(models_folder,'model_cp_%05d.sps' % (index+1)))

    # for index in range(n_models):
    #     cfx_model = Surfpack('CFX_%05d' % (index+1), desvar.ndim)
    #     cfx_model.load_model(os.path.join(models_folder,'model_cfx_%05d.sps' % (index+1)))

    # for index in range(n_models):
    #     cfy_model = Surfpack('CFY_%05d' % (index+1), desvar.ndim)
    #     cfy_model.load_model(os.path.join(models_folder,'model_cfy_%05d.sps' % (index+1)))

    # for index in range(n_models):
    #     cfz_model = Surfpack('CFZ_%05d' % (index+1), desvar.ndim)
    #     cfz_model.load_model(os.path.join(models_folder,'model_cfz_%05d.sps' % (index+1)))

    # Compute Structure

    number_dsn = 15

#    for index in range(number_dsn*10-1,number_dsn*10+9):
    for index in [partitions-1]:

        dvs = X[index]

        # dvs = [1.1,0.2450089127,0.866,0.806527,-1.804901,338617.115272,13415.473225,22250,-0.5,0.5,0.5,0.0,0.0,0.0] # ~ DSN_003_A1

        desvar.unpack_structure(konfig, dvs)

        proc = project.func('AERODYNAMICS', konfig)

        # force_redo_dsn_folder = None #'DSN_001'
        # proc = project.func('STRUCTURE', konfig, force_redo_dsn_folder)

        #proc.wait()



#: response_surface()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
