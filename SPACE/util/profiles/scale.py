#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    scale = 1.09

    data = np.loadtxt('saved/edge_1.dat')
    plt.plot(data[:,0]*scale,data[:,1]*scale, '.', color='b', lw=2)

    np.savetxt('struct/edge_1_new.dat',data*scale)

    data = np.loadtxt('fluid/edge_3.dat')
    plt.plot(data[:,0],data[:,1], '.', color='r', lw=2)

    data = np.loadtxt('fluid/edge_23.dat')
    plt.plot(data[:,0],data[:,1], '.', color='g', lw=2)

    data = np.loadtxt('fluid/edge_45.dat')
    plt.plot(data[:,0],data[:,1], '.', color='k', lw=2)




    plt.axis('equal')
    plt.show()
						
		
#: def main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
