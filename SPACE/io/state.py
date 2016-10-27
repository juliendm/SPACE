#!/usr/bin/env python

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, time
from ..io   import expand_part, expand_time, get_adjointSuffix, add_suffix, \
                   get_specialCases, Config
from ..util import bunch
from ..util import ordered_bunch


# ----------------------------------------------------------------------
#  State Factory
# ----------------------------------------------------------------------

def State_Factory(state=None,config=None):
    
    if isinstance(state,Config) and not config:
        config = state
        state = None
    
    if not state is None:
        assert isinstance(state,State) , 'input is must be a state instance'
        return state
    
    NewClass = State()
    
    for key in ['FUNCTIONS','GRADIENTS','VARIABLES','FILES','HISTORY']:
        NewClass[key] = ordered_bunch()
            
    if config:
        NewClass.find_files(config)
            
    return NewClass


# ----------------------------------------------------------------------
#  State Class
# ----------------------------------------------------------------------

class State(ordered_bunch):
    
    _timestamp = 0
    
    def update(self,ztate):
        """ Updates self given another state
        """
        
        if not ztate: return
        assert isinstance(ztate,State) , 'must update with another State-type'
        for key in self.keys():
            if isinstance(ztate[key],dict):
                self[key].update( ztate[key] )
            elif ztate[key]:
                self[key] = ztate[key]
        
        self.set_timestamp()
                
        
    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        output = 'STATE:'
        for k1,v1 in self.iteritems():
            output += '\n    %s:' % k1
            if isinstance(v1,dict):
                for k2,v2 in v1.iteritems():
                    output += '\n        %s: %s' % (k2,v2)
            else:
                output += '\n        %s' % v1
        return output
    
    def pullnlink(self,config):
        """ pull,link = SU2.io.State.pullnlink(config)
            returns lists pull and link of files for folder
            redirection, based on a given config
        """
        
        pull = []; link = []
        
        # choose files to pull and link
        for key,value in self.FILES.iteritems():
            pull.append(value)
        #: for each filename
        
        return pull,link
    
    def design_vector(self):
        """ vectorizes State.VARIABLES
        """
        vector = []
        for value in self.VARIABLES.values():
            if isinstance(value,dict):
                for v in value.values():
                    vector.append(v)
            elif not isinstance(value,list):
                value = [value]
            vector.extend(value)
        return vector
    
    def find_files(self,config):
        """ SU2.io.State.find_files(config)
            finds mesh and solution files for a given config.
            updates state.FILES with filenames.
            files already logged in state are not overridden.
            will ignore solutions if config.RESTART_SOL == 'NO'.
        """
        
        files = self.FILES
        
        fluid_boundary_name     = config.FLUID_BOUNDARY_FILENAME
        correspondance_name     = config.CORRESPONDANCE_FILENAME
        config_aero_name        = config.CONFIG_AERO_FILENAME
        
        def register_file(label,filename):
            if not files.has_key(label):
                if os.path.exists(filename):
                    files[label] = filename
                    print 'Found: %s' % filename
            else:
                assert os.path.exists(files[label]) , 'state expected file: %s' % filename
        #: register_file()                

        # mesh
        register_file('FLUID_BOUNDARY',fluid_boundary_name)
        register_file('CORRESPONDANCE',correspondance_name)
        register_file('CONFIG_AERO',config_aero_name)
        
        return
    
    def __setitem__(self,k,v):
        if self._initialized:
            self.set_timestamp()
        super(State,self).__setitem__(k,v)
    
    def set_timestamp(self):
        self._timestamp = time.time()
    
    def tic(self):
        """ timestamp = State.tic() 
            returns the time that this state was last modified
        """
        return self._timestamp
    
    def toc(self,timestamp):
        """ updated = State.toc(timestamp)
            returns True if state was modified since last timestamp
        """
        return self._timestamp > timestamp
        
    
#: def State
