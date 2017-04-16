#!/usr/bin/env python 

import time, os, gc, sys, shutil, copy, math
import numpy as np
import pwd
import pickle

from optparse import OptionParser

from Surfpack.surfpack import Surfpack

from pySBO.pyGPR import LHC_unif

sys.path.append(os.environ['MISSION_ANALYSIS_RUN'])
import MAP

def main():

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--precalc", dest="precalulate", default="False",
                      help="precalulate", metavar="PRECALCULATE")

    (options, args)=parser.parse_args()
    options.precalulate = options.precalulate.upper() == 'TRUE'

    response_surface( options.filename    ,
                        options.precalulate )

#: main()

def response_surface( filename           ,
                        precalulate       ):

    # Project
    project = MAP.project.Project(filename)

    XB = np.array([[-0.5, 0.5],
                   [-0.5, 0.5]])
    ndim = len(XB)

    if precalulate:

        nd = 10*ndim

        X = LHC_unif(XB,nd)
        PERFO = np.zeros(nd)

        for i in range(nd):
            dvs = X[i]
            dummy, PERFO[i] = project.eval(dvs[0],0.0,dvs[1])

        out_perfo = open("build_points/build_points_perfo.dat", "w")
        for i in range(nd):
            string = "%.11f %.11f %.11f\n" % (X[i][0], X[i][1], PERFO[i])
            out_lift.write(string)
        out_perfo.close()

    else:

        na = 10

        X = None
        PERFO = None

        # Kriging
        os.system('surfpack spk/sp_build_surrogate_perfo.spk')

        surfpack = Surfpack(ndim)

        for ite in range(na):
        
            surfpack.load('model_'+str(ite), 'sps/sp_gp_model.perfo.sps')

            dvs = surfpack.max_variance('model_'+str(ite), XB)

            print '-------------------------------'
            print dvs
            print '-------------------------------'

            if X != None: X = np.append(X,[dvs], axis=0)
            else: X = np.array([dvs])

            dummy, perfo = project.eval(dvs[0],0.0,dvs[1])

            if PERFO != None: PERFO = np.append(PERFO, [perfo], axis=0)
            else: PERFO = np.array([perfo])

            inp_perfo = open("build_points/build_points_perfo.dat")
            out_perfo = open("build_points/temp_perfo.dat", "w")
            out_perfo.writelines([l for l in inp_perfo.readlines()])
            for i in range(len(PERFO)):
                string = "%.11f %.11f %.11f\n" % (X[i][0], X[i][1], PERFO[i])
                out_perfo.write(string)
            inp_perfo.close()
            out_perfo.close()
            os.rename("build_points/temp_perfo.dat", "build_points/enrich_points_perfo.dat")

            os.system('surfpack spk/sp_enrich_surrogate_perfo.spk')


#: response_surface()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
