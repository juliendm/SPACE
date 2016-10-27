#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    
    in_file = 'slices.dat'
    inp = open(in_file)

    slice = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
    n_slice = 45 #len(slice)
	
    #---------------------

    inp.readline()
    inp.readline()
    inp.readline()
    inp.readline()

    x_l_previous = None
    x_u_previous = None
    y_l_previous = None
    y_u_previous = None

    plt.figure(1)
	
    for i in range(n_slice):
	
        inp.readline()
        inp.readline()
        line = inp.readline()
        nNode = int(line.split('Nodes=')[-1].split(',')[0])
        nElem = int(line.split('Elements=')[-1].split(',')[0])
        print nNode
        inp.readline()
        inp.readline()

        node = [[0 for z in range(3)] for k in range(nNode)]
        for k in range(nNode):
            data = inp.readline().split()
            node[k][0] = float(data[0])
            node[k][1] = float(data[1])
            if (node[k][1] < 0.001): node[k][1] = 0.0
            node[k][2] = float(data[2])
			
        elem = [[0 for z in range(2)] for k in range(nElem)]
        for k in range(nElem):
            data = inp.readline().split()
            elem[k][0] = int(data[0])
            elem[k][1] = int(data[1])
			
        ordered_x = []
        ordered_z = []
		
        current_pos = 0
        index_x = 0
        index_z = 2
		
        initial_index = elem[current_pos][0]-1

        current_index = elem[current_pos][0]-1
        ordered_x.append(node[current_index][index_x])
        ordered_z.append(node[current_index][index_z])
		
        current_index = elem[current_pos][1]-1
        ordered_x.append(node[current_index][index_x])
        ordered_z.append(node[current_index][index_z])

        for n in range(nElem-2):
            for k in range(nElem):
                if ((elem[k][0]-1 == current_index) and (k != current_pos)):
                    current_pos = k
                    current_index = elem[current_pos][1]-1
                    break
                elif ((elem[k][1]-1 == current_index) and (k != current_pos)):
                    current_pos = k
                    current_index = elem[current_pos][0]-1
                    break
            if (current_index == initial_index): break
            ordered_x.append(node[current_index][index_x])
            ordered_z.append(node[current_index][index_z])

        length = len(ordered_x)
				
        index_min_x = ordered_x.index(min(ordered_x))
        index_max_x = ordered_x.index(max(ordered_x))
				
        min_index = min(index_min_x,index_max_x)
        max_index = max(index_min_x,index_max_x)

	X_u = np.zeros(max_index+1-min_index)
        Y_u = np.zeros(max_index+1-min_index)

        X_l = np.zeros(length+2-len(X_u))
        Y_l = np.zeros(length+2-len(Y_u))		

        for k in range(min_index,max_index+1):
            X_u[k-min_index] = ordered_x[k]
            Y_u[k-min_index] = ordered_z[k]
        for k in range(max_index,length):
            X_l[k-max_index] = ordered_x[k]
            Y_l[k-max_index] = ordered_z[k]
        for k in range(min_index+1):
            X_l[k+length-max_index] = ordered_x[k]
            Y_l[k+length-max_index] = ordered_z[k]

        if i == 0:
            min_x = min([min(X_u),min(X_l)])
            min_y = 0.5*(Y_u[np.argmin(X_u)]+Y_l[np.argmin(X_l)])

        X_u = X_u - min_x
        X_l = X_l - min_x
        Y_u = Y_u - min_y
        Y_l = Y_l - min_y

        if X_l[0] < X_l[1]:
            X_u = X_u[::-1]
            Y_u = Y_u[::-1]
        else:
            X_l = X_l[::-1]
            Y_l = Y_l[::-1]

#        if i in [42,43]:
#            temp = Y_u
#            Y_u = Y_l
#            Y_l = temp

        wing_length = 14.5

        if i < 0:

            cutting_index_u = 0
            for k in range(len(X_u)):
                if (X_u[k] + min_x > wing_length) :
                    cutting_index_u = k
                    break

            cutting_index_l = 0
            for k in range(len(X_l)):
                if (X_l[k] + min_x > wing_length) :
                    cutting_index_l = k
                    break

            if cutting_index_u > 0 :
                a_u = (Y_u[cutting_index_u-1]-Y_u[cutting_index_u])/(X_u[cutting_index_u-1]-X_u[cutting_index_u])
                b_u = Y_u[cutting_index_u] - a_u*X_u[cutting_index_u]
                X_u = X_u[0:cutting_index_u+1]
                X_u[-1] = wing_length-min_x
                Y_u = Y_u[0:cutting_index_u+1]
                Y_u[-1] = a_u*X_u[-1] + b_u 

            if cutting_index_l > 0 :
                a_l = (Y_l[cutting_index_l-1]-Y_l[cutting_index_l])/(X_l[cutting_index_l-1]-X_l[cutting_index_l])
                b_l = Y_l[cutting_index_l] - a_l*X_l[cutting_index_l]
                X_l = X_l[0:cutting_index_l+1]
                X_l[-1] = wing_length-min_x
                Y_l = Y_l[0:cutting_index_l+1]
                Y_l[-1] = a_l*X_l[-1] + b_l

#        angle = 0.0 * np.pi / 180.0
#        X_u_p = X_u; Y_u_p = Y_u
#        X_u = X_u_p * np.cos(angle) - Y_u_p * np.sin(angle)
#        Y_u = X_u_p * np.sin(angle) + Y_u_p * np.cos(angle)
#        X_l_p = X_l; Y_l_p = Y_l
#        X_l = X_l_p * np.cos(angle) - Y_l_p * np.sin(angle)
#        Y_l = X_l_p * np.sin(angle) + Y_l_p * np.cos(angle)

        max_x = max([max(X_u),max(X_l)])
#        X_u = X_u / max_x
#        Y_u = Y_u / max_x
#        X_l = X_l / max_x
#        Y_l = Y_l / max_x
#        X_u[-1] = 1.0
#        X_l[-1] = 1.0




        profile = np.loadtxt('../GeoMACH-master/GeoMACH/PGM/airfoils/rae2822.dat')
        x_u = profile[0:65,0]
        x_l = profile[0:65,0]
        x_u = x_u*(np.max(X_u)-np.min(X_u))+np.min(X_u)
        x_l = x_l*(np.max(X_l)-np.min(X_l))+np.min(X_l)
        x_u[0] = X_u[0]
        x_l[0] = X_l[0]
        x_u[-1] = X_u[-1]       
        x_l[-1] = X_l[-1]

        f_u = interpolate.interp1d(X_u,Y_u)
        f_l = interpolate.interp1d(X_l,Y_l)
##        s_u = interpolate.UnivariateSpline(x_u,f_u(x_u),s=1)
##        s_l = interpolate.UnivariateSpline(x_l,f_l(x_l),s=1)
        y_u = f_u(x_u)
        y_l = f_l(x_l)

#        x_u = X_u
#        x_l = X_l
#        y_u = Y_u
#        y_l = Y_l

#        x_l = -x_l+max_x
#        x_l = x_l[::-1]
#        x_u = -x_u+max_x
#        x_u = x_u[::-1]
#        y_l = y_l[::-1]
#        y_u = y_u[::-1]

        avg_begin = 0.5*(y_u[0]+y_l[0])
        y_u[0] = avg_begin
        y_l[0] = avg_begin

        n_smooth = 12
        y_target = y_l[-1]
#        y_target = -0.05
        a_u_l = np.array([[x_u[n_smooth-1]*x_u[n_smooth-1],x_u[n_smooth-1],1],[x_u[0]*x_u[0],x_u[0],1],[2.0*x_u[n_smooth-1],1,0]])
        b_u_l = np.array([y_u[n_smooth-1],y_target,(y_u[n_smooth-1]-y_u[n_smooth])/(x_u[n_smooth-1]-x_u[n_smooth])])
        a_u_r = np.array([[x_u[-n_smooth]*x_u[-n_smooth],x_u[-n_smooth],1],[x_u[-1]*x_u[-1],x_u[-1],1],[2.0*x_u[-n_smooth],1,0]])
        b_u_r = np.array([y_u[-n_smooth],y_target,(y_u[-n_smooth]-y_u[-n_smooth-1])/(x_u[-n_smooth]-x_u[-n_smooth-1])])
        coeffs_u_l = np.linalg.solve(a_u_l,b_u_l)
        coeffs_u_r = np.linalg.solve(a_u_r,b_u_r)
        a_l_l = np.array([[x_l[n_smooth-1]*x_l[n_smooth-1],x_l[n_smooth-1],1],[x_l[0]*x_l[0],x_l[0],1],[2.0*x_l[n_smooth-1],1,0]])
        b_l_l = np.array([y_l[n_smooth-1],y_target,(y_l[n_smooth-1]-y_l[n_smooth])/(x_l[n_smooth-1]-x_l[n_smooth])])
        a_l_r = np.array([[x_l[-n_smooth]*x_l[-n_smooth],x_l[-n_smooth],1],[x_l[-1]*x_l[-1],x_l[-1],1],[2.0*x_l[-n_smooth],1,0]])
        b_l_r = np.array([y_l[-n_smooth],y_target,(y_l[-n_smooth]-y_l[-n_smooth-1])/(x_l[-n_smooth]-x_l[-n_smooth-1])])
        coeffs_l_l = np.linalg.solve(a_l_l,b_l_l)
        coeffs_l_r = np.linalg.solve(a_l_r,b_l_r)
        for ind in range(n_smooth):
            index_l = ind
            index_r = -(ind+1)
#            y_u[index_l] = coeffs_u_l[0]*x_u[index_l]*x_u[index_l] + coeffs_u_l[1]*x_u[index_l] + coeffs_u_l[2] 
            y_u[index_r] = coeffs_u_r[0]*x_u[index_r]*x_u[index_r] + coeffs_u_r[1]*x_u[index_r] + coeffs_u_r[2]
#            y_l[index_l] = coeffs_l_l[0]*x_l[index_l]*x_l[index_l] + coeffs_l_l[1]*x_l[index_l] + coeffs_l_l[2]            
            y_l[index_r] = coeffs_l_r[0]*x_l[index_r]*x_l[index_r] + coeffs_l_r[1]*x_l[index_r] + coeffs_l_r[2]

        out_file = 'profile_' + str(i+1) + '.dat'
        out = open(out_file,"w")
        for k in range(len(x_u)):
            out.write(str(x_u[len(x_u)-1-k]) + ' ' + str(y_u[len(x_u)-1-k]) + '\n')
        for k in range(1,len(x_l)):
            out.write(str(x_l[k]) + ' ' + str(y_l[k]) + '\n')
        out.close()



        if i in [0,2,22,44]:

            plt.plot(x_u,y_u, '.', color='r', lw=2)
            plt.plot(x_l,y_l, '.', color='b', lw=2)

            if not x_l_previous is None:
                for m in range(len(x_l)) :
                    if m%10 == 0:
                        plt.plot([x_l_previous[m],x_l[m]],[y_l_previous[m],y_l[m]],color = 'black')
                        plt.plot([x_u_previous[m],x_u[m]],[y_u_previous[m],y_u[m]],color = 'black')

            out_file = 'edge_' + str(i+1) + '.dat'
            out = open(out_file,"w")
            for k in range(len(x_u)):
                out.write(str(x_u[k]) + ' ' + str(y_u[k]) + '\n')
            out.write('\n')
            for k in range(len(x_l)):
                out.write(str(x_l[k]) + ' ' + str(y_l[k]) + '\n')
            out.close()

            x_l_previous = x_l
            x_u_previous = x_u
            y_l_previous = y_l
            y_u_previous = y_u

    plt.axis('equal')
    plt.show()
				
		
#: def main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
