#!/usr/bin/env python2.7

import time, os, gc, sys, shutil, copy, math
import numpy as np

sys.path.append(os.environ['SPACE_RUN'])
import SPACE

# from pySBO.pyGPR import LHC_unif    

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    # Config
    config = SPACE.io.Config('config.cfg')
    # State
    state  = SPACE.io.State()
    #state.find_files(config)

    # Project
    project = SPACE.project.Project(config, state, None, 'MAX_VEL_DVS_SPACE')

    konfig = copy.deepcopy(config)

    XB = np.array([[-0.5, 0.5]])
    dvs = np.linspace(-0.5,0.1,num=3)   #LHC_unif(XB,20)
#    dvs = np.append(dvs, [[0.5,0.5],[0.5,-0.5],[-0.5,0.5],[-0.5,-0.5]], axis=0)
#    dvs = np.array([[0.5,0.5],[0.5,-0.5],[-0.5,0.5],[-0.5,-0.5]])
#    dvs = np.array([[0.0, 0.0]])

    procs = []
    for index in range(0,len(dvs)):
        konfig.DV1 = '0.0'; konfig.DV2 = '0.0'; konfig.DV3 = str(dvs[index])
        procs.append(project.func('MISSION', konfig)) # , v_i, m_o, m_i, max_pdyn, max_mach, max_alpha 
#        print mass_dry, dv1, v_o, m_o
        # print dv1, v_o

    for proc in procs:
        proc.wait()

#: main()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
