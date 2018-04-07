import time, os, gc, sys, shutil, copy, math
import numpy as np
from ctypes import cdll, c_double

SPACE_RUN = os.environ['SPACE_RUN']
surfpack_lib = cdll.LoadLibrary(SPACE_RUN + '/SPACE/surfpack/lib/surfpack.so')

import pyOpt

class Surfpack(object):

    def __init__(self,name,ndim):

        self._eval = surfpack_lib.eval
        self._eval.restype = c_double

        self._gradient = surfpack_lib.gradient
        self._gradient.restype = c_double * ndim

        self._variance = surfpack_lib.variance
        self._variance.restype = c_double

        self.name = name
        self.ndim = ndim
        self.array = c_double * ndim

    def load_data(self,data):
        
        surfpack_lib.load_data(self.name,data,self.ndim,1,0)

    def save_data(self,data):

        surfpack_lib.save_data(self.name,data)

    def load_model(self,model):
        
        surfpack_lib.load_model(self.name,model)

    def save_model(self,model):

        surfpack_lib.save_model(self.name,model)

    def add(self,dvs,f):

        if isinstance(dvs, np.ndarray): dvs = dvs.tolist()
        assert len(dvs) == self.ndim, 'wrong dimension'

        c_dvs = self.array()
        for ix in range(self.ndim):
            c_dvs[ix] = dvs[ix];

        surfpack_lib.add(self.name,c_dvs,self.ndim,c_double(f))

    def build(self,type):

        surfpack_lib.build(self.name,type)

    def eval(self,dvs):

        if isinstance(dvs, np.ndarray): dvs = dvs.tolist()
        assert len(dvs) == self.ndim, 'wrong dimension'

        c_dvs = self.array()
        for ix in range(self.ndim):
            c_dvs[ix] = dvs[ix];

        return self._eval(self.name,c_dvs,self.ndim)

    def gradient(self,dvs):

        if isinstance(dvs, np.ndarray): dvs = dvs.tolist()
        assert len(dvs) == self.ndim, 'wrong dimension'

        c_dvs = self.array()
        c_grad = self.array()

        for ix in range(self.ndim):
            c_dvs[ix] = dvs[ix];

        self._gradient(self.name,c_dvs,c_grad,self.ndim)

        grad = np.zeros(self.ndim)
        for ix in range(self.ndim):
            grad[ix] = c_grad[ix];

        return grad

    def variance(self,dvs):

        if isinstance(dvs, np.ndarray): dvs = dvs.tolist()
        assert len(dvs) == self.ndim, 'wrong dimension'

        c_dvs = self.array()
        for ix in range(self.ndim):
            c_dvs[ix] = dvs[ix];

        return self._variance(self.name,c_dvs,self.ndim)

    def max_variance(self,XB, number = 1, exclude = []):

        assert XB.shape[0] == self.ndim, 'wrong dimension'
        
        X_min = [0.0]*len(XB)
        X_min[0] = XB[0][0]-1.0

        while not is_in(X_min,XB):

            obj = self.objective
            prob = pyOpt.Optimization('Variance Maximization',obj)
                    
            for ix in range(XB.shape[0]):
                prob.addVar('X%i'%ix,'c',lower=XB[ix,0],upper=XB[ix,1],value=0.)
                
            prob.addObj('Estimated Variance')

            opt_ALPSO = pyOpt.ALPSO(pll_type=None)
            #opt_ALPSO = pyOpt.ALPSO(pll_type='MP',args=[1.0])
            opt_ALPSO.setOption('fileout',0)
            opt_ALPSO.setOption('maxOuterIter',10)
            opt_ALPSO.setOption('stopCriteria',1)       
#            opt_ALPSO.setOption('SwarmSize',self.ndim*100)
            opt_ALPSO.setOption('SwarmSize',self.ndim*20)

            opt_SLSQP = pyOpt.SLSQP()
            opt_SLSQP.setOption('IPRINT',-1)
            opt_SLSQP.setOption('ACC',1e-5)

            vec = []
            for index in range(number*10):
                print index+1, ' so far: ', len(vec)
                [YI_min,X_min,Info] = opt_ALPSO(prob)
                [YI_min,X_min,Info] = opt_SLSQP(prob.solution(index),sens_type='FD')
                if not is_already_in(X_min,vec+exclude) and is_in(X_min,XB): vec.append(X_min.tolist())
                if len(vec) >= number: break

        if len(vec) == 1: return vec[0]
        return vec

    def objective(self,XC):
        F = [-self.variance(XC)]
        G = []
        if len(F)==1: F=F[0]
        fail = 0
        return F,G,fail 


def is_already_in(X,vec):

    for in_X in vec:
        if len([i for i, j in zip(X, in_X) if i == j]) == len(X):
            return True
    return False

def is_in(X,XB):

    eps = 1e-6

    for index,val in enumerate(X):
        if (val < XB[index][0]-eps) or (val > XB[index][1]+eps):
            return False
    return True
