#!/usr/bin/env python2.7

import time, os, gc, sys, shutil, copy, math
import numpy as np

import matplotlib.pyplot as plt

sys.path.append(os.environ['SPACE_RUN'])
import SPACE

def main():

  project_folder = 'RESPONSE_SURFACE_DV_SUP'
  project = SPACE.io.load_data(project_folder + '/project.pkl')
  project.compile_designs()

  for design_container in project.designs:
    design = design_container.design
    config = design.config

    print design.folder

    if hasattr(design.funcs,'LIFT'):

      history_filename = os.path.join(project.folder,'DESIGNS',design.folder,'DIRECT','history_direct.dat')

      data = np.loadtxt(history_filename, skiprows=3, delimiter=',')

      fig = plt.figure()

      plt.plot(data[:,0],data[:,1],data[:,0],data[:,2],data[:,0],data[:,5])

      plt.xlabel('Iterations', fontsize=18)
      plt.ylabel('Coeff', fontsize=18)

      plt.legend(['LIFT','DRAG','MOMENT_Y'], loc='upper right')

      plt.figtext(0.2, 0.75, design.folder + '\nMach = ' + config.MACH_NUMBER + '\nReynods = ' + config.REYNOLDS_NUMBER + '\nAoA = ' + config.AoA)



      figure_filename = os.path.join('figures','history_direct_' + design.folder + '.png')
      fig.savefig(figure_filename)

      plt.close(fig)

#: main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
