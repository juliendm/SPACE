#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    
    plt.figure(1)
	
    profile = np.loadtxt('../GeoMACH-master/GeoMACH/PGM/airfoils/rae2822.dat')
    x_u = profile[0:65,0]
    x_l = profile[0:65,0]

    y_u = 0.05*np.ones(len(x_u))

    n_smooth = 18
    y_target = 0
    a_u_l = np.array([[x_u[n_smooth-1]*x_u[n_smooth-1],x_u[n_smooth-1],1],[x_u[0]*x_u[0],x_u[0],1],[2.0*x_u[n_smooth-1],1,0]])
    b_u_l = np.array([y_u[n_smooth-1],y_target,(y_u[n_smooth-1]-y_u[n_smooth])/(x_u[n_smooth-1]-x_u[n_smooth])])
    a_u_r = np.array([[x_u[-n_smooth]*x_u[-n_smooth],x_u[-n_smooth],1],[x_u[-1]*x_u[-1],x_u[-1],1],[2.0*x_u[-n_smooth],1,0]])
    b_u_r = np.array([y_u[-n_smooth],y_target,(y_u[-n_smooth]-y_u[-n_smooth-1])/(x_u[-n_smooth]-x_u[-n_smooth-1])])
    coeffs_u_l = np.linalg.solve(a_u_l,b_u_l)
    coeffs_u_r = np.linalg.solve(a_u_r,b_u_r)
    for ind in range(n_smooth):
        index_l = ind
        y_u[index_l] = coeffs_u_l[0]*x_u[index_l]*x_u[index_l] + coeffs_u_l[1]*x_u[index_l] + coeffs_u_l[2] 
        index_r = -(ind+1)
        y_u[index_r] = coeffs_u_r[0]*x_u[index_r]*x_u[index_r] + coeffs_u_r[1]*x_u[index_r] + coeffs_u_r[2]

    y_l = -y_u 

    out_file = 'profile_cs.dat'
    out = open(out_file,"w")
    for k in range(len(x_u)):
        out.write(str(x_u[len(x_u)-1-k]) + ' ' + str(y_u[len(x_u)-1-k]) + '\n')
    for k in range(1,len(x_l)):
        out.write(str(x_l[k]) + ' ' + str(y_l[k]) + '\n')
    out.close()


    plt.plot(x_u,y_u, '.', color='r', lw=2)
    plt.plot(x_l,y_l, '.', color='b', lw=2)

    plt.axis('equal')
    plt.show()
				
		
#: def main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
