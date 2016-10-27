#!/usr/bin/env python 

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, sys, shutil, copy, glob, time
import numpy as np
from .. import io   as spaceio
from .. import eval as spaceeval
from .. import util as spaceutil
from ..io import redirect_folder

from warnings import warn, simplefilter
#simplefilter(Warning,'ignore')

inf = 1.0e20


# -------------------------------------------------------------------
#  Project Class
# -------------------------------------------------------------------

class Project(object):
    
    _design_folder = 'DESIGNS/DSN_*'
    _design_number = '%03d'
    
    def __init__( self, config, state=None , 
                  designs=None, folder='.' ,
                  warn = True                ):
        
        folder = folder.rstrip('/')+'/'
        if '*' in folder: folder = spaceio.next_folder(folder)        
        if designs is None: designs = []
        
        print 'New Project: %s' % (folder)
        
        # setup config
        config = copy.deepcopy(config)
        
        # setup state
        if state is None:
            state = spaceio.State()
        else:
            state  = copy.deepcopy(state)
            state  = spaceio.State(state)

        state.find_files(config)

        if 'FLUID_BOUNDARY' not in state.FILES:
            raise Exception , 'Could not find mesh file: %s' % config.FLUID_BOUNDARY_FILENAME
        if 'CORRESPONDANCE' not in state.FILES:
            raise Exception , 'Could not find dat file: %s' % config.CORRESPONDANCE_FILENAME
        if 'CONFIG_AERO' not in state.FILES:
            raise Exception , 'Could not find config file: %s' % config.CONFIG_AERO_FILENAME
        
        self.config  = config      # base config
        self.state   = state       # base state
        self.files   = state.FILES # base files
        self.designs = designs     # design list
        self.folder  = folder      # project folder
        self.results = spaceutil.ordered_bunch() # project design results
        
        # output filenames
        self.filename = 'project.pkl' 
        self.results_filename = 'results.pkl' 
        
        # initialize folder with files
        pull,link = state.pullnlink(config)
        with redirect_folder(folder,pull,link,force=True):
        
            # look for existing designs
            folders = glob.glob(self._design_folder)
            if len(folders)>0:
                sys.stdout.write('Removing old designs in 10s.')
                sys.stdout.flush()
                if warn: time.sleep(10)
                sys.stdout.write(' Done!\n\n')
                for f in folders: shutil.rmtree(f)
            #: if existing designs
            
            # save project
            spaceio.save_data(self.filename,self)
            
        return
    
    def _eval(self,config,func,*args):
        """ evalautes a config, checking for existing designs
        """
        
        konfig = copy.deepcopy(config) # design config
        config = self.config           # project config
        state  = self.state            # project state
        folder = self.folder           # project folder
        
        filename = self.filename
        
        # check folder
        assert os.path.exists(folder) , 'cannot find project folder %s' % folder        
        
        # list project files to pull and link
        pull,link = state.pullnlink(config)
        
        # project folder redirection, don't overwrite files
        with redirect_folder(folder,pull,link,force=False) as push:        
        
            # start design
            design = self.new_design(konfig)
            
            if config.get('CONSOLE','VERBOSE') == 'VERBOSE':
                print os.path.join(self.folder,design.folder)
            timestamp = design.state.tic()
            
            # run design+
            vals = design._eval(func,*args)
            
            # check for update
            if design.state.toc(timestamp):
                
                # recompile design results
                self.compile_results()

                # save data
                spaceio.save_data(filename,self)
                
            #: if updated
            
        #: with redirect folder
        
        # done, return output
        return vals
    
    def unpack_dvs(self,dvs):
        dvs = copy.deepcopy(dvs)
        konfig = copy.deepcopy( self.config )
        if isinstance(dvs, np.ndarray): dvs = dvs.tolist()
        konfig.unpack_dvs(dvs)
        return konfig, dvs
    
    def func(self,func_name,config):
        func = spaceeval.func
        konfig = copy.deepcopy(config)
        return self._eval(konfig, func, func_name)
    
    def grad(self,func_name,method,config):
        func = spaceeval.grad
        konfig = copy.deepcopy(config)
        return self._eval(konfig, func, func_name,method)
    
    def user(self,user_func,config,*args):
        raise NotImplementedError
        #return self._eval(config, user_func,*args) 
    
    def add_design(self,config):
        #func = spaceeval.touch # hack - TWL
        func = spaceeval.skip 
        konfig = copy.deepcopy(config)
        return self._eval(konfig, func)
        
    def new_design(self,config):
        """ finds an existing design for given config
            or starts a new design with a closest design 
            used for restart data
        """
         # local konfig
        konfig = copy.deepcopy(config)
        
        # find closest design
        closest,delta = self.closest_design(konfig)
        # found existing design
        if delta == 0.0 and closest:
            design = closest
        # start new design
        else:
            design = self.init_design(konfig,closest)
        #: if new design    
        
        return design
    
    def get_design(self,config):
        konfig = copy.deepcopy(config)
        closest,delta = self.closest_design(konfig)
        if delta == 0.0 and closest:
            design = closest
        else:
            raise Exception, 'design not found for this config'
        return design
        
    def closest_design(self,config):
        """ looks for an existing or closest design 
            given a config
        """        
                
        designs = self.designs
        
        keys_check = ['MACH_NUMBER','AoA','DV_VALUE_NEW']
        
        if not designs: 
            return [] , inf
        
        diffs = []
        for this_design in designs:
            this_config = this_design.config
            distance = config.dist(this_config,keys_check)
            diffs.append(distance) 
                        
        #: for each design 
        
        # pick closest design
        i_min = np.argmin(diffs)
        delta  = diffs[i_min]
        closest = designs[i_min]
        
        return closest, delta 
    
    def init_design(self,config,closest=None):
        """ starts a new design
            works in project folder
        """
        
        konfig = copy.deepcopy(config)
        ztate  = copy.deepcopy(self.state)
        if closest is None: closest = []
        
        # use closest design as seed
        if closest:
            # copy useful state info
            seed_folder = closest.folder
            seed_files  = closest.files
            for key in seed_files.keys():
                # ignore mesh
                if key == 'MESH': continue 
                # build file path
                name = seed_files[key]
                name = os.path.join(seed_folder,name)
                # update pull files
                ztate.FILES[key] = name
            
        # name new folder
        folder = self._design_folder.replace('*',self._design_number)
        folder = folder % (len(self.designs) + 1)

        # start new design (pulls files to folder)
        design = spaceeval.Design(konfig,ztate,folder)
        
        # update local state filenames ( ??? why not in Design() )
        for key in design.files:
            name = design.files[key]
            name = os.path.split(name)[-1]
            design.files[key] = name
        
        # add design to project 
        self.designs.append(design)        
        
        return design
    
    def compile_results(self,default=np.nan):
        """ results = SU2.opt.Project.compile_results(default=np.nan)
            builds a Bunch() of design results
            
            Inputs:
                default - value for missing values
                
            Outputs:
                results - state with items filled with list of
                values ordered by each design iteration.
                
                results.VARIABLES
                results.FUNCTIONS
                results.GRADIENTS
                results.HISTORY.DIRECT
                results.HISTORY.ADJOINT_*
                
        """
        
        results = spaceio.State()
        results.VARIABLES = []
        del results.FILES
        filename = self.results_filename
        
        n_dv = 0
        
        # populate fields
        for i,design in enumerate(self.designs):
            for key in design.state.FUNCTIONS.keys():
                results.FUNCTIONS[key] = []
            for key in design.state.GRADIENTS.keys():
                results.GRADIENTS[key] = []
            for TYPE in design.state.HISTORY.keys():
                if not results.HISTORY.has_key(TYPE):
                    results.HISTORY[TYPE] = spaceutil.ordered_bunch()
                for key in design.state.HISTORY[TYPE].keys():
                    results.HISTORY[TYPE][key] = []
            this_ndv = len( design.state.design_vector() )
            
            # check design vectors are of same length
            if i == 0:
                n_dv = this_ndv
            else:
                if n_dv != this_ndv:
                    warn('different dv vector length during compile_results()')
        #: for each design
            
        # populate results
        for design in self.designs:
            this_designvector = design.state.design_vector()
            results.VARIABLES.append( this_designvector )
            for key in results.FUNCTIONS.keys():
                if design.state.FUNCTIONS.has_key(key):
                    new_func = design.state.FUNCTIONS[key]
                else:
                    new_func = default
                results.FUNCTIONS[key].append(new_func)
            for key in results.GRADIENTS.keys():
                if design.state.GRADIENTS.has_key(key):
                    new_grad = design.state.GRADIENTS[key]
                else:
                    new_grad = [default] * len( this_designvector )
                results.GRADIENTS[key].append(new_grad)
            for TYPE in results.HISTORY.keys():
                for key in results.HISTORY[TYPE].keys():
                    if key in results.FUNCTIONS.keys():
                        new_func = results.FUNCTIONS[key][-1]
                    elif ( TYPE in design.state.HISTORY.keys() and
                            key in design.state.HISTORY[TYPE].keys() ):
                        new_func = design.state.HISTORY[TYPE][key][-1]
                    else:
                        new_func = default
                    results.HISTORY[TYPE][key].append(new_func)
        #: for each design
        
        # save
        self.results = results
        spaceio.save_data(filename,results)
            
        return self.results
    
    def deep_compile(self):
        """ Project.deep_compile()
            recompiles project using design files saved in each design folder
            useful if designs were run outside of project class
        """
        
        project_folder = self.folder
        designs = self.designs
        
        with spaceio.redirect_folder(project_folder):
            for i_dsn,design in enumerate(designs):
                design_filename = os.path.join(design.folder,design.filename)
                self.designs[i_dsn] = spaceio.load_data(design_filename)
            
            self.compile_results()
            spaceio.save_data(self.filename,self)
            
        return
        
    def save(self):
        with spaceio.redirect_folder(self.folder):
            spaceio.save_data(self.filename,self)
        
    def __repr__(self):
        return '<Project> with %i <Design>' % len(self.designs)
    def __str__(self):
        output = self.__repr__()
        return output    
