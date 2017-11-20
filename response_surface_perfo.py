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
from SPACE.eval import model as spacemodel

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
                      initiate = True  ,
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


    # Start new model

    konfig = copy.deepcopy(config)

    aero = spacemodel.AeroModel(konfig,project_folder)

    # Design Variables

    XB = np.array([
        [0.3,  8.0],       # dv_mach
        [-1.0, 1.0],       # dv_rey
        [-1.0, 3.0],       # dv_aoa            # CHEAT UPPER BOUND AOA

        [-0.5, 0.5],       # dv_geo1
        [-0.5, 0.5],       # dv_geo2
        [-0.5, 0.5],       # dv_geo3
        [-0.2, 0.5],       # dv_geo4
    ])

    ndim = 7

    flag = 'TRIM_AFT'

    if initiate:

        nd = 10*ndim

        X = LHC_unif(XB,nd)

        model = Surfpack(flag,ndim)

        for index in range(len(X)):

            perfo_dvs = X[index]
            aero_dvs = perfo_to_aero(perfo_dvs)
            val = eval_function(flag,aero,aero_dvs)
            model.add(perfo_dvs,val)

        model.save_data(os.path.join(project_folder,'build_points_' + flag + '.dat'))

    else:

        na = 500

        model = Surfpack(flag,ndim)
        model.load_data(os.path.join(project_folder,'build_points_' + flag + '.dat'))
        model.build('kriging')
        model.save_model(os.path.join(project_folder,'model_' + flag + '.sps'))

        for ite in range(na):

            print 'Ite:', ite
            model.build('kriging')
            print 'Model built'

            perfo_dvs = model.max_variance(XB)
            print 'dvs: ', perfo_dvs
            aero_dvs = perfo_to_aero(perfo_dvs)
            val = eval_function(flag,aero,aero_dvs)
            model.add(perfo_dvs,val)

            model.save_data(os.path.join(project_folder,'enriched_points_' + flag + '.dat'))

        # Save Model

        model.build('kriging')
        model.save_model(os.path.join(project_folder,'model_' + flag + '.sps'))

#: response_surface()

def eval_function(flag,aero,aero_dvs):

    cog_x_aft = 0.63*17
    cog_x_fwd = 0.60*17
    cog_z = -0.3

    if flag == 'TRIM_FWD':
        val = aero.raw_trim(aero_dvs,cog_x_fwd,cog_z)
    elif flag == 'TRIM_AFT':
        val = aero.raw_trim(aero_dvs,cog_x_aft,cog_z)
    elif (flag == 'TRIMMED_STATIC_MARGIN_FWD'):
        trimmed_dvs = aero.trim(aero_dvs,cog_x_fwd,cog_z)
        val = aero.static_margin(trimmed_dvs,cog_x_fwd,cog_z)
    elif (flag == 'TRIMMED_STATIC_MARGIN_AFT'):
        trimmed_dvs = aero.trim(aero_dvs,cog_x_aft,cog_z)
        val = aero.static_margin(trimmed_dvs,cog_x_aft,cog_z)
    elif (flag == 'K_ALPHA_FWD'):
        trimmed_dvs = aero.trim(aero_dvs,cog_x_fwd,cog_z)
        val = aero.k_alpha(trimmed_dvs,cog_x_fwd,cog_z)
    elif (flag == 'K_ALPHA_AFT'):
        trimmed_dvs = aero.trim(aero_dvs,cog_x_aft,cog_z)
        val = aero.k_alpha(trimmed_dvs,cog_x_aft,cog_z)

    return val

def perfo_to_aero(dvs):

    return [dvs[0],dvs[1],dvs[2], 0.0,0.0, dvs[3],dvs[4],dvs[5],dvs[6], 0.5,0.5]

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
