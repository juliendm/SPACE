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

    # Command Line Options
    parser=OptionParser()
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")
                      
    (options, args)=parser.parse_args()
    options.partitions  = int( options.partitions )

    # Project
    project_folder = 'RESPONSE_SURFACE_DV_SUP'
    if os.path.exists(project_folder):
        project = SPACE.io.load_data(project_folder + '/project.pkl')
        project.compile_designs()
        config = project.config
    else:
        config = SPACE.io.Config(options.filename)
        state  = SPACE.io.State()
        project = SPACE.project.Project(config, state, folder=project_folder)

    print '%d design(s) so far' % len(project.designs)

    konfig = copy.deepcopy(config)
    konfig.NUMBER_PART = options.partitions

    procs = []


    # dv_geo1  # [-0.5;0.5]
    # dv_geo2  # [-0.5;0.5]
    # dv_geo3  # [-0.5;0.5]
    # dv_geo4  # [-0.3;0.5]
    # dv_geo5  # [-0.5;0.5]
    # dv_geo6  # [-0.3;0.5]

    # dv_bf    # [-0.25;0.5]
    # dv_el    # [-0.35;0.5]

    # dv_mach  # [0.3;8.0]
    # dv_re    # [-1.0;1.0]
    # dv_aoa   # [-1.0;1.0]


    if False:

        XB = np.array([

            [0.3, 0.95],        # dv_mach
            [-1.0, 1.0],       # dv_re
            [-1.0, 1.0],       # dv_aoa

            [-0.25, 0.5],      # dv_bf
            [-0.35, 0.5],      # dv_el

            [-0.5, 0.5],       # dv_geo1
            [-0.5, 0.5],       # dv_geo2
            [-0.5, 0.5],       # dv_geo3
            [-0.2, 0.5],       # dv_geo4
            [-0.5, 0.5],       # dv_geo5
            [-0.2, 0.5]        # dv_geo6

        ])

        ndim = len(XB)
        nd = 10*ndim

        X = LHC_unif(XB,nd)

        np.savetxt('dvs.dat', X, fmt='%.18e', delimiter=', ', newline='\n', header='', footer='', comments='# ')

    else:

        X = np.loadtxt('dvs_sup.dat', delimiter=', ', comments='# ')

        for index in range(len(X)):

            dvs = X[index]

            dv_mach = dvs[0]
            dv_re = dvs[1]
            dv_aoa = dvs[2]
            dv_bf = dvs[3]
            dv_el = dvs[4]
            dv_geo1 = dvs[5]
            dv_geo2 = dvs[6]
            dv_geo3 = dvs[7]
            dv_geo4 = dvs[8]
            dv_geo5 = dvs[9]
            dv_geo6 = dvs[10]

            Reynolds = 10.0**(-3.0/8.0*dv_mach+7.0+dv_re)

            if dv_mach >= 1.0:
                AoA = 3.92857*dv_mach+3.57143+dv_aoa*(1.78571*dv_mach+5.71429) # Supersonic
            else:
                AoA = (dv_aoa+1.0)*7.5 # Subsonic

            konfig.DV1 = str(dv_geo1); konfig.DV2 = str(dv_geo2); konfig.DV3 = str(dv_geo3); konfig.DV4 = str(dv_geo4); konfig.DV5 = str(dv_geo5); konfig.DV6 = str(dv_geo6)
            konfig.ELEVON_DEF = str(dv_el*180.0/np.pi); konfig.BODY_FLAP_DEF = str(dv_bf*180.0/np.pi)
            konfig.MACH_NUMBER = str(dv_mach); konfig.REYNOLDS_NUMBER = str(Reynolds); konfig.AoA= str(AoA)
            proc = project.func('AERODYNAMICS', konfig)
            proc.wait()




    for proc in procs:
        #if not proc is None: proc.Disconnect()
        if not proc is None: proc.wait()

    project.compile_designs()

def LHC_unif(XB,NI,XI=None,maxits=10):
    ''' Latin Hypercube Sampling with uniform density
        iterates to maximize minimum L2 distance
    '''
    
    # dimension
    ND = XB.shape[0]
    
    # initial points to respect
    if XI is None:
        XI = np.empty([0,ND])
       
    # output points
    XO = []
    
    # initialize
    mindiff = 0;
    
    # maximize minimum distance
    for it in range(maxits):
        
        # samples
        S = np.zeros([NI,ND])
        
        # populate samples
        for i_d in range(ND):
            S[:,i_d] = ( np.random.random([1,NI]) + np.random.permutation(NI) ) / NI
        XS = S*(XB[:,1]-XB[:,0]) + XB[:,0]        
        
        # add initial points
        XX = np.vstack([ XI , XS ])
        
        # calc distances
        vecdiff = VecDist(XX)[0]
        
        # update
        if vecdiff > mindiff:
            mindiff = vecdiff
            XO = XX
        
    #: for iterate
    
    return XO

def VecDist(X,P=None):
    ''' calculates distance between points in matrix X 
        with each other, or optionally to given point P
        returns min, max and matrix/vector of distances
    '''
    
    # distance matrix among X
    if P is None:
        
        nK,nD = X.shape
        
        d = np.zeros([nK,nK,nD])
        for iD in range(nD):
            d[:,:,iD] = np.array([X[:,iD]])-np.array([X[:,iD]]).T
        D = np.sqrt( np.sum( d**2 , 2 ) )
        
        diag_inf = np.diag( np.ones([nK])*np.inf )
        dmin = np.min(np.min( D + diag_inf ))
        dmax = np.max(np.max( D ))
        
    # distance vector to P
    else:
        assert P.shape[0] == 1 , 'P must be a horizontal vector'
        D = np.array([ np.sqrt( np.sum( (X-P)**2 , 1 ) ) ]).T
        dmin = D.min()
        dmax = D.max()
        
    return (dmin,dmax,D)

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
	main()
