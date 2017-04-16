#!/usr/bin/env python2.7

import time, os, gc, sys, shutil, copy, math
import numpy as np
import pwd
import pickle

from optparse import OptionParser

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE.surfpack import Surfpack

#from pySBO.pyGPR import LHC_unif

sys.path.append(os.environ['SU2_RUN'])
import SU2

from SU2.io import Config, State
from SU2.eval import func, grad
from SU2.opt.project import Project

def main():

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--precalc", dest="precalulate", default="False",
                      help="precalulate", metavar="PRECALCULATE")
    parser.add_option("-n", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")


    (options, args)=parser.parse_args()
    options.precalulate = options.precalulate.upper() == 'TRUE'
    options.partitions = int( options.partitions )

    response_surface( options.filename    ,
                        options.precalulate ,
                        options.partitions  )

#: main()

def response_surface( filename           ,
                        precalulate      ,
                        partitions  = 0  ):



    # Config
    config = Config(filename)
    config.NUMBER_PART = partitions
    config.CONSOLE = 'CONCISE'
    config.EXT_ITER = 350
#    config.EXT_ITER = 600

    # State
    state  = State()

    # Project
    project = Project(config, state, None, 'SUP_ADD')


    XB = np.array([[1.2, 9.0],
           [0.0, 20.0],
           [-0.5, 0.5],
           [-0.5, 0.5],
           [-0.5, 0.5]])

    ndim = len(XB)

    if precalulate:

        nd = 10*ndim

        X = LHC_unif(XB,nd)

        lift_model = Surfpack('LIFT',ndim)
        drag_model = Surfpack('DRAG',ndim)

        for i in range(nd):

            dvs = X[i]

            konfig = copy.deepcopy(config)
                    
            konfig.MACH_NUMBER = dvs[0]
            if dvs[0] > 8.0: konfig.CFL_NUMBER = 1.1
            konfig.AoA = dvs[1]
            konfig.unpack_dvs([dvs[2],dvs[3],dvs[4]])

            lift_model.add(dvs,project.func('LIFT',konfig))
            drag_model.add(dvs,project.func('DRAG',konfig))

        lift_model.save_data('build_points_lift.dat')
        drag_model.save_data('build_points_drag.dat')

    else:

        na = 30

        lift_model = Surfpack('LIFT',ndim)
        lift_model.load_data('build_points_lift.dat')

        drag_model = Surfpack('DRAG',ndim)
        drag_model.load_data('build_points_drag.dat')

        # Objective: LIFT
        lift_model.build('kriging')

        for ite in range(na):

            dvs = lift_model.max_variance(XB)

            print '-------------------------------'
            print dvs
            print '-------------------------------'

            konfig = copy.deepcopy(config)
            
            konfig.MACH_NUMBER = dvs[0]
            if dvs[0] > 8.0: konfig.CFL_NUMBER = 1.1
            konfig.AoA = dvs[1]
            konfig.unpack_dvs([dvs[2],dvs[3],dvs[4]])

            lift_model.add(dvs,project.func('LIFT',konfig))
            drag_model.add(dvs,project.func('DRAG',konfig))

            # Objective: LIFT
            lift_model.build('kriging')

        lift_model.save_data('enriched_points_lift.dat')
        drag_model.save_data('enriched_points_drag.dat')

        lift_model.save_model('lift.sps')
        drag_model.save_model('drag.sps')

#: response_surface()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
