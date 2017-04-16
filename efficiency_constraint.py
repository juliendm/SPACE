#!/usr/bin/env python2.7

import time, os, gc, sys, shutil, copy, math
import numpy as np

from optparse import OptionParser

sys.path.append(os.environ['SPACE_RUN'])
import SPACE

from SPACE import io   as spaceio
from SPACE import util as spaceutil
from SPACE.eval import model as spacemodel

import matplotlib.pyplot as plt

# from pySBO.pyGPR import LHC_unif    

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    config = spaceio.Config('config.cfg')

    data = spaceutil.ordered_bunch()
    data.lift_sub = spaceutil.ordered_bunch()
    data.lift_sub.sps = config.DATABASE_LIFT_SUB
    data.lift_sub.ranges = np.array([[0.5,0.95],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
    data.lift_sup = spaceutil.ordered_bunch()
    data.lift_sup.sps = config.DATABASE_LIFT_SUP
    data.lift_sup.ranges = np.array([[1.1,9.0],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
    lift_model = spacemodel.RangedModel(data)

    data = spaceutil.ordered_bunch()
    data.drag_sub = spaceutil.ordered_bunch()
    data.drag_sub.sps = config.DATABASE_DRAG_SUB
    data.drag_sub.ranges = np.array([[0.5,0.95],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
    data.drag_sup = spaceutil.ordered_bunch()
    data.drag_sup.sps = config.DATABASE_DRAG_SUP
    data.drag_sup.ranges = np.array([[1.1,9.0],[-5.0,20.0],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])
    drag_model = spacemodel.RangedModel(data)

    mach = 0.5
    aoa = 10.0

    # XB = np.array([[-0.5, 0.5],[-0.5, 0.5]]) 
    # dvs = LHC_unif(XB,600)

    x = np.arange(-0.5, 0.5, 0.01)
    y = np.arange(-0.5, 0.5, 0.01)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros(X.shape)

    for i in range(X.shape[0]):
        for j in range(X.shape[1]):

            model_l = lift_model.find_model([mach, aoa, 0.0, X[i][j], Y[i][j]])
            model_d = drag_model.find_model([mach, aoa, 0.0, X[i][j], Y[i][j]])
            Z[i][j] = model_l.eval([mach, aoa, 0.0, X[i][j], Y[i][j]]) / model_d.eval([mach, aoa, 0.0, X[i][j], Y[i][j]])

    fig = plt.figure()
    plt.contour(X,Y,Z)
    fig.savefig('image.png', dpi=fig.dpi)

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

    main()
