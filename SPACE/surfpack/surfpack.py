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

    def max_variance(self,XB):

        assert XB.shape[0] == self.ndim, 'wrong dimension'
        
        obj = self.objective
        prob = pyOpt.Optimization('Variance Maximization',obj)
                
        for ix in range(XB.shape[0]):
            prob.addVar('X%i'%ix,'c',lower=XB[ix,0],upper=XB[ix,1],value=0.)
            
        prob.addObj('Estimated Variance')

        opt = pyOpt.ALPSO(pll_type=None)
        #opt = pyOpt.ALPSO(pll_type='MP',args=[1.0])
        opt.setOption('fileout',0)
        opt.setOption('maxOuterIter',10)
        opt.setOption('stopCriteria',1)       
        opt.setOption('SwarmSize',self.ndim*100)
        opt(prob)

        opt = pyOpt.SLSQP()
        opt.setOption('IPRINT',-1)
        opt.setOption('ACC',1e-5)
        [YI_min,X_min,Info] = opt(prob.solution(0),sens_type='FD')

        return X_min

    def objective(self,XC):
        F = [-self.variance(XC)]
        G = []
        if len(F)==1: F=F[0]
        fail = 0
        return F,G,fail 

