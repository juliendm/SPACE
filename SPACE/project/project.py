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


import subprocess
#from multiprocessing import Process
#from mpi4py import MPI    ############## MUST NOT USE SHARED MEMORY: ALL EXE MUST BE COMPILED IN SERIAL
SPACE_RUN = os.environ['SPACE_RUN']



# -------------------------------------------------------------------
#  Project Class
# -------------------------------------------------------------------

class Project(object):
    
    _design_folder = 'DESIGNS/DSN_*'
    _design_number = '%03d'
    
    def __init__( self, config, state=None , 
                  designs=None, models=None , folder='.' ,
                  warn = True                ):
        
        folder = folder.rstrip('/')+'/'
        if '*' in folder: folder = spaceio.next_folder(folder)        
        
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

        # check has the needed files
        if ('MODEL_LIFT_SUP' in config.keys()):
            if 'MODEL_LIFT_SUP' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.MODEL_LIFT_SUP
        if ('MODEL_DRAG_SUP' in config.keys()):
            if 'MODEL_DRAG_SUP' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.MODEL_DRAG_SUP
        if ('MODEL_FORCE_X_SUP' in config.keys()):
            if 'MODEL_FORCE_X_SUP' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.MODEL_FORCE_X_SUP
        if ('MODEL_FORCE_Z_SUP' in config.keys()):
            if 'MODEL_FORCE_Z_SUP' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.MODEL_FORCE_Z_SUP
        if ('MODEL_MOMENT_Y_SUP' in config.keys()):
            if 'MODEL_MOMENT_Y_SUP' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.MODEL_MOMENT_Y_SUP

        if ('MODEL_LIFT_SUB' in config.keys()):
            if 'MODEL_LIFT_SUB' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.MODEL_LIFT_SUB
        if ('MODEL_DRAG_SUB' in config.keys()):
            if 'MODEL_DRAG_SUB' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.MODEL_DRAG_SUB
        if ('MODEL_FORCE_X_SUB' in config.keys()):
            if 'MODEL_FORCE_X_SUB' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.MODEL_FORCE_X_SUB
        if ('MODEL_FORCE_Z_SUB' in config.keys()):
            if 'MODEL_FORCE_Z_SUB' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.MODEL_FORCE_Z_SUB
        if ('MODEL_MOMENT_Y_SUB' in config.keys()):
            if 'MODEL_MOMENT_Y_SUB' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.MODEL_MOMENT_Y_SUB

        if ('CONFIG_AERO_FILENAME' in config.keys()):
            if 'FARFIELD' not in state.FILES:
                raise Exception , 'Could not find mesh file: %s' % config.FARFIELD_FILENAME
        if ('CONFIG_AERO_FILENAME' in config.keys()):
            if 'CONFIG_AERO' not in state.FILES:
                raise Exception , 'Could not find config file: %s' % config.CONFIG_AERO_FILENAME
        

        self.config  = config      # base config
        self.state   = state       # base state
        self.files   = state.FILES # base files


        if designs is None: designs = []
        self.designs = designs     # design list
        if models is None: models = spaceutil.ordered_bunch()
        self.models  = models      # model

        self.folder  = folder      # project folder
        self.results = spaceutil.ordered_bunch() # project design results



        
        # output filenames
        self.filename = 'project.pkl' 
        self.results_filename = 'results.pkl' 
        



        # initialize folder with files
        pull,link = state.pullnlink(True,config)
        with redirect_folder(folder,pull,link,force=True):
            
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
        pull,link = state.pullnlink(False,config)
        
        # project folder redirection, don't overwrite files
        with redirect_folder(folder,pull,link,force=False) as push:        
        
            if not args[1] is None:
                for design_container in self.designs:
                    if (design_container.folder == args[1].strip()): break
            else:
                # start design
                design_container = self.new_design_container(konfig)
            
            if config.get('CONSOLE','VERBOSE') == 'VERBOSE':
                print os.path.join(self.folder,design_container.folder)

            if (design_container.design is None) and (args[1] is None):
                #timestamp = design.state.tic()

                # run design: initialize folder with files
                pull,link = state.pullnlink(False,config)
                with redirect_folder(design_container.folder,pull,link,force=True):
                    konfig.dump('config_DSN.cfg')

                    Command = 'python2.7 ' + SPACE_RUN + '/SPACE/eval/design_interface.py ' + args[0]
                    proc = subprocess.Popen(Command, shell=True, stdout=sys.stdout, stderr=subprocess.PIPE)

                    # proc = Process(target=spaceeval.eval_design, args=(args[0], config))
                    # proc.start()

                    #proc = MPI.COMM_SELF.Spawn(sys.executable, args=[SPACE_RUN + '/SPACE/eval/design_interface.py',args[0]], maxprocs=1)
                
                # # check for update
                # if design.state.toc(timestamp):
                #     # recompile design results
                #     self.compile_results()

                # update data with new design_container (if subprocess is stopped, project will still know this design was ongoing)
                spaceio.save_data(filename,self)

            else:
                with redirect_folder(design_container.folder,force=False):
                    Command = 'python2.7 ' + SPACE_RUN + '/SPACE/eval/design_interface.py ' + args[0]
                    proc = subprocess.Popen(Command, shell=True, stdout=sys.stdout, stderr=subprocess.PIPE)

        # done, return output
        return proc
    
    def unpack_dvs(self,dvs):
        dvs = copy.deepcopy(dvs)
        konfig = copy.deepcopy( self.config )
        if isinstance(dvs, np.ndarray): dvs = dvs.tolist()
        konfig.unpack_dvs(dvs)
        return konfig, dvs
    




    def func(self,func_name,config,dsn_folder=None):
        func = spaceeval.func
        konfig = copy.deepcopy(config)
        return self._eval(konfig, func, func_name, dsn_folder)
    










        
    def new_design_container(self,config):
        """ finds an existing design for given config
            or starts a new design with a closest design 
            used for restart data
        """
         # local konfig
        konfig = copy.deepcopy(config)
        
        # find closest design
        # closest,delta = self.closest_design(konfig)
        # found existing design

        # IMPROVE
        closest = None
        delta = 1e8

        if delta == 0.0 and closest:
            design_container = closest
        else:
            # name new folder
            folder = self._design_folder.replace('*',self._design_number)
            folder = folder % (len(self.designs) + 1)
            design_container = spaceutil.ordered_bunch()
            design_container.folder = folder
            design_container.design = None
            self.designs.append(design_container)      
        
        return design_container
        
    def closest_design(self,config):
        """ looks for an existing or closest design 
            given a config
        """        
                
        designs = self.designs
        
        keys_check = ['MACH_NUMBER','AoA','DV1','DV2','DV3']
        
        if not designs: 
            return [] , inf
        
        diffs = []
        for this_design_container in designs:
            if not this_design_container.design is None:
                this_config = this_design_container.design.config
                distance = config.dist(this_config,keys_check)
                diffs.append(distance)
        
        # pick closest design
        i_min = np.argmin(diffs)
        delta  = diffs[i_min]
        closest = designs[i_min]
        
        return closest, delta 


    def compile_designs(self, force=False, update_name=False):
        """
            recompiles project using design files saved in each design folder
        """
        
        project_folder = self.folder
        designs = self.designs

        with spaceio.redirect_folder(project_folder):
            for design_container in designs:
                if ((design_container.design is None) or force):
                    design_filename = os.path.join(design_container.folder,'design.pkl')
                    if os.path.exists(design_filename):
                        design = spaceio.load_data(design_filename)
                        if update_name:
                            design.folder = design_container.folder.split('/')[-1]
                            spaceio.save_data(design_filename,design)
                        design_container.design = design
            
            #self.compile_results()
            spaceio.save_data(self.filename,self)
            
        return

    def fresh_compile_designs(self, project_folder=None):
        """
            recompiles project using design files saved in each design folder
        """
        
        if not project_folder is None:
            self.folder = project_folder
        else:
            project_folder = self.folder

        self.designs = []

        ls_dsn = glob.glob(project_folder + '/DESIGNS/*')
        ls_dsn.sort()

        for dsn in ls_dsn:
            print dsn
            design_container = spaceutil.ordered_bunch()
            design_container.folder = 'DESIGNS/' + dsn.split('/')[-1]
            design_container.design = None
            self.designs.append(design_container)      

        self.compile_designs(update_name=True)

        return



    def compile_models(self, force=False):
        """
            compile models
        """
            
        return

    









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
    

        
    def save(self):
        with spaceio.redirect_folder(self.folder):
            spaceio.save_data(self.filename,self)
        
    def __repr__(self):
        return '<Project> with %i <Design>' % len(self.designs)
    def __str__(self):
        output = self.__repr__()
        return output    
