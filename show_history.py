#!/usr/bin/env python2.7

import time, os, gc, sys, shutil, copy, math
import numpy as np
from optparse import OptionParser

import matplotlib.pyplot as plt

sys.path.append(os.environ['SPACE_RUN'])
import SPACE

def main():

  # Command Line Options
  parser=OptionParser()
  parser.add_option("-p", "--project", dest="project_folder",
                    help="project folder", metavar="PROJECT_FOLDER")
  parser.add_option("-f", "--figures", dest="figures_folder",
                    help="figures folder", metavar="FIGURES_FOLDER")
                    
  (options, args)=parser.parse_args()

  show_history( options.project_folder ,
                options.figures_folder )

#: main()

def show_history( project_folder ,
                  figures_folder ):

  project = SPACE.io.load_data(project_folder + '/project.pkl')
  project.compile_designs()

  for design_container in project.designs:
    design = design_container.design
    
    if design is not None:

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

        plt.figtext(0.2, 0.75, design.folder + '\nMach = ' + config.MACH_NUMBER + '\nReynods = ' + config.REYNOLDS_NUMBER + '\nAoA = ' + config.AoA + '\nBF = ' + config.BODY_FLAP_DEF + '\nEL = ' + config.ELEVON_DEF)



        figure_filename = os.path.join(figures_folder,'history_direct_' + design.folder + '.png')
        fig.savefig(figure_filename)

        plt.close(fig)

#: main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__=='__main__':
    main()
