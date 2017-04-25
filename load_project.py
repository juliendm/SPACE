#!/usr/bin/env python2.7

import os, time, sys, shutil, copy
import numpy as np
from optparse import OptionParser
import random

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE.surfpack import Surfpack

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from SPACE.util import DesignVariables

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    # Command Line Options
    parser=OptionParser()
    parser.add_option("-p", "--project", dest="project_folder",
                      help="project folder", metavar="PROJECT_FOLDER")
    parser.add_option("-r", "--regime", dest="regime", default="SUP",
                      help="regime", metavar="REGIME")

    (options, args)=parser.parse_args()

    load_project( options.project_folder ,
                  options.regime )

#: main()

def load_project( project_folder   ,
                  regime = 'SUP' ):

    desvar = DesignVariables(regime)
    lift_model = Surfpack('LIFT',desvar.ndim)

    project = SPACE.io.load_data(os.path.join(project_folder,'project.pkl'))
    project.compile_designs()

    for index in range(len(project.designs)):
 
        design = project.designs[index].design
    
        funcs = design.funcs
        config = design.config
        dvs = desvar.pack(config)

        if hasattr(funcs,'LIFT'):
            lift_model.add(dvs,funcs.LIFT)


    lift_model.save_data('build_points_lift_sup.dat')
    lift_model.build('kriging')
    lift_model.save_model('lift_sup.sps')

    # lift_model.load_model('lift.sps')

    nx = 50

    fig = plt.figure()
    ax = Axes3D(fig)


    model = lift_model
    model_range = [-0.1,0.6]

    mach_range = [[ 1.1, 1.5],[ 1.5, 2.0],[ 2.0, 2.5],[ 2.5, 3.0],[ 3.0, 4.0],[ 4.0, 5.0],[ 5.0, 8.0]]
    aoa_range  = [[ 0.0,15.0],[ 2.5,20.0],[ 5.0,25.0],[ 7.5,27.5],[10.0,30.0],[12.5,40.0],[15.0,50.0]]

    for index in range(len(mach_range)):

        MACH = np.linspace(mach_range[index][0], mach_range[index][1], nx)
        AOA = np.linspace(aoa_range[index][0], aoa_range[index][1], nx)
        MACH,AOA = np.meshgrid(MACH,AOA)
        COEFF = MACH*0.

        for ix in range(nx):
            for iy in range(nx):

                dv_mach = MACH[ix,iy]
                dv_re, dv_aoa = desvar.reverse_variable_change(dv_mach, 1E7, AOA[ix,iy])
                COEFF[ix,iy] = model.eval([dv_mach,0.0,dv_aoa,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

        ax.plot_surface(MACH,AOA,COEFF, norm=colors.Normalize(vmin=model_range[0],vmax=model_range[1]), cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)

    plt.xlabel('Mach', fontsize=18)
    plt.ylabel('AoA', fontsize=18)

    fig.savefig('fig.png')   # save the figure to file
    plt.close(fig)    # close the figure

    # plt.draw()
    # plt.show()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()
