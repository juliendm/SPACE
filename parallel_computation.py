#!/usr/bin/env python2.7

import os, time, sys, shutil, copy
import numpy as np
from optparse import OptionParser
import random

sys.path.append(os.environ['SPACE_RUN'])
import SPACE

from SPACE.util import LHC_unif, DesignVariables

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    # Command Line Options
    parser=OptionParser()
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--project", dest="project_folder",
                      help="project folder", metavar="PROJECT_FOLDER")
    parser.add_option("-r", "--regime", dest="regime", default="SUP",
                      help="regime", metavar="REGIME")
    parser.add_option("-n", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")
                      
    (options, args)=parser.parse_args()
    options.partitions  = int( options.partitions )

    parallel_computation( options.filename        ,
                          options.project_folder  ,
                          options.regime          ,
                          options.partitions  )

#: main()

def parallel_computation( filename           ,
                          project_folder     ,
                          regime = 'SUP'     ,
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
    konfig.NUMBER_PART = partitions

    procs = []

    desvar = DesignVariables(regime)

    if False:

        nd = 10*desvar.ndim

        X = LHC_unif(desvar.XB,nd)

        np.savetxt(os.path.join(project_folder,'dvs_new.dat'), X, fmt='%.18e', delimiter=', ', newline='\n', header='', footer='', comments='# ')

    else:

        X = np.loadtxt(os.path.join(project_folder,'dvs.dat'), delimiter=', ', comments='# ')

        for index in [79]: #range(0,len(X)):

            dvs = X[index]
            desvar.unpack(konfig, dvs)
            
            proc = project.func('AERODYNAMICS', konfig)
            proc.wait()


    for proc in procs:
        #if not proc is None: proc.Disconnect()
        if not proc is None: proc.wait()

    project.compile_designs()



# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()
