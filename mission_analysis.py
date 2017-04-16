#!/usr/bin/env python2.7

import time, os, gc, sys, shutil, copy, math
import numpy as np

from optparse import OptionParser

sys.path.append(os.environ['SPACE_RUN'])
import SPACE

# from pySBO.pyGPR import LHC_unif    

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main(compute):

    # Project
    project_folder = 'MAX_VEL_DV3'
    if os.path.exists(project_folder):
        project = SPACE.io.load_data(project_folder + '/project.pkl')
        project.compile_designs(force=True)
        config = project.config
    else:
        config = SPACE.io.Config('config.cfg')
        state  = SPACE.io.State()
        project = SPACE.project.Project(config, state, folder=project_folder)

    print '%d design(s) so far' % len(project.designs)

    if compute:

        konfig = copy.deepcopy(config)

        #XB = np.array([[-0.5, 0.5],[-0.5, 0.5],[-0.5, 0.5]]) # LHC_unif(XB,60)
        dvs = np.linspace(-0.1,0.5,num=12)

        comms = []
        for index in range(len(dvs)):
            konfig.DV1 = str(0.0); konfig.DV2 = str(dvs[index]); konfig.DV3 = str(0.0) #dvs[index][2]
            comms.append(project.func('MISSION', konfig))
        for comm in comms:
            if not comm is None: comm.Disconnect()

        project.compile_designs()

    else:

        for index in range(len(project.designs)):
            design = project.designs[index].design
            if "SPEED_MECO" in design.state.FUNCTIONS.keys():
                print design.config.DV1, design.config.DV2, design.config.DV3, design.state.FUNCTIONS.SPEED_MECO

#: main()

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
if __name__=='__main__':

    # Command Line Options
    parser=OptionParser()
    parser.add_option("-c", "--compute",    dest="compute",    default="True",
                      help="COMPUTE direct and adjoint problem", metavar="COMPUTE")
                      
    (options, args)=parser.parse_args()
    options.compute     = options.compute.upper() == 'TRUE'

    main(options.compute)
