#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np
from optparse import OptionParser

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    corresp = np.loadtxt('correspondance.dat', dtype='int')
    inp_original = open('fluid_surface_original.mesh')
    inp_deformed = open('fluid_surface.mesh')
    out = open('mesh_motion.dat',"w")

    # Read inp Write out
    #---------------------

    # Points    

    line = inp_original.readline()
    while line[0] != 'V':
        line = inp_original.readline()
    line = inp_original.readline()
    nPoint = int(line.split()[0])
    original = []
    for index in range(nPoint):
        line = inp_original.readline()
        data = line.split()
        original.append([float(data[0]), float(data[1]), float(data[2])])

    line = inp_deformed.readline()
    while line[0] != 'V':
        line = inp_deformed.readline()
    line = inp_deformed.readline()
    deformed = []
    for index in range(nPoint):
        line = inp_deformed.readline()
        data = line.split()
        deformed.append([float(data[0]), float(data[1]), float(data[2])])

    for index in range(nPoint):
        out.write("%d %f %f %f\n" % (corresp[index], deformed[index][0], deformed[index][1], deformed[index][2]))

    inp_original.close()
    inp_deformed.close()
    out.close()

 
#: def main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
