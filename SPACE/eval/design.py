#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, glob, re
from .. import io   as spaceio
from .  import func as spacefunc
from ..io import redirect_folder, save_data

# ----------------------------------------------------------------------
#  Design Class
# ----------------------------------------------------------------------

class Design(object):

    def __init__(self, config, state=None):
        """ Initializes a SPACE Design """
        
        config = copy.deepcopy(config)
        state  = copy.deepcopy(state)
        state  = spaceio.State(state)
        state.find_files(config)
        
        self.config = config
        self.state  = state
        self.files  = state.FILES
        self.funcs  = state.FUNCTIONS
        self.grads  = state.GRADIENTS
        self.folder = os.getcwd().split('/')[-1]
        
        self.filename = 'design.pkl'

        # save design, config
        save_data(self.filename,self)
        
    def _eval(self,eval_func,*args):
        """ Evaluates a SPACE Design 
            always adds config and state to the inputs list
        """
        
        config = self.config
        state  = self.state
        files  = self.files
        folder = self.folder
        
        filename = self.filename
            
        # get timestamp
        timestamp = state.tic()
        
        # run 
        inputs = args + (config,state)
        vals = eval_func(*inputs)

        # save design
        if state.toc(timestamp):
            save_data(filename,self)
        
        # update files
        files.update(state['FILES'])
        
        return vals

    def func(self,func_name):
        """ Evaluates SU2 Design Functions by Name """
        return self._eval(spacefunc,func_name)
    
    def touch(self):
        return self._eval(touch)
    
    def skip(self,*args,**kwarg):
        return self._eval(skip)
        
    
    def __repr__(self):
        return '<Design> %s' % self.folder
    def __str__(self):
        output = self.__repr__()
        output += '\n%s' % self.state
        return output
    
#: class Design()

def touch(config,state):
    """ SU2.eval.touch(config,state)
        resets state timestamp 
    """
    state.set_timestamp()
    
def skip(config,state):
    """ SU2.eval.skip(config,state)
        does nothing
    """
    pass
