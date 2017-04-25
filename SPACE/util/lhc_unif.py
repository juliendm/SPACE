#!/usr/bin/env python 

import os, time, sys, shutil
import numpy as np

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

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
    