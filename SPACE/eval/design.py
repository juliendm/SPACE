#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, glob, re
from .. import io   as spaceio
from .  import func as spacefunc
from ..io import redirect_folder, save_data

import subprocess

SPACE_RUN = os.environ['SPACE_RUN']

# ----------------------------------------------------------------------
#  Design Class
# ----------------------------------------------------------------------

class Design(object):

    def __init__(self, config, state=None, folder='DESIGNS/DSN_*'):
        """ Initializes a SPACE Design """
        
        if '*' in folder: folder = spaceio.next_folder(folder)
        
        config = copy.deepcopy(config)
        state  = copy.deepcopy(state)
        state  = spaceio.State(state)
        state.find_files(config)
        
        self.config = config
        self.state  = state
        self.files  = state.FILES
        self.funcs  = state.FUNCTIONS
        self.grads  = state.GRADIENTS
        self.folder = folder
        
        self.filename = 'design.pkl'
            
        # initialize folder with files
        pull,link = state.pullnlink(config)
        with redirect_folder(folder,pull,link,force=True):
            # save design, config
            save_data(self.filename,self)
            config.dump('config_DSN.cfg')
        
    def _eval(self,eval_func,*args):
        """ Evaluates a SPACE Design 
            always adds config and state to the inputs list
        """
        
        config = self.config
        state  = self.state
        files  = self.files
        folder = self.folder
        
        filename = self.filename

        # check folder
        assert os.path.exists(folder) , 'cannot find design folder %s' % folder
        
        # list files to pull and link
        pull,link = state.pullnlink(config)
        
        # output redirection, don't re-pull files
        with redirect_folder(folder,pull,link,force=False) as push:
            
            # get timestamp
            timestamp = state.tic()
            
            # run 
            inputs = args + (config,state)
            #vals = eval_func(*inputs)

            # PUT EVAL_FUNC IN CONFIG AND MODIFY FUNCTIONS.PY
            Command = 'python2.7 ' + SPACE_RUN + '/SPACE/eval/functions.py'
            proc = subprocess.Popen( Command, shell=True    ,
                             stdout=sys.stdout      , 
                             stderr=subprocess.PIPE  )

            
            # save design
            if state.toc(timestamp):
                save_data(filename,self)
            
        #: with redirect folder
        
        # update files
        files.update(state['FILES'])
        
        return proc

    def func(self,func_name):
        """ Evaluates SU2 Design Functions by Name """
        return self._eval(su2func,func_name)
    
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
