#!/usr/bin/env python2.7

import os, time, sys, shutil, copy
import numpy as np
from optparse import OptionParser
import random

sys.path.append(os.environ['SPACE_RUN'])
import SPACE

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    data = np.loadtxt('RESPONSE_SURFACE_DV_SUP/dvs_sup_new.dat', delimiter=', ', comments='# ')

    dvs_geo2 = []

    for index in range(len(data)):

        # data[index][6] = data[index][6] - data[index][8]	
        dvs_geo2.append(data[index][6])

    print min(dvs_geo2), max(dvs_geo2)
    print dvs_geo2.index(min(dvs_geo2)), dvs_geo2.index(max(dvs_geo2))

    print [i for i,v in enumerate(dvs_geo2) if v <= -0.8]

    # np.savetxt('RESPONSE_SURFACE_DV_SUB/dvs_sub_new.dat', data, delimiter=', ')

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()

