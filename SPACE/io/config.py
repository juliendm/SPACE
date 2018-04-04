#!/usr/bin/env python

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import numpy as np
from ..util import bunch, ordered_bunch, switch
from .tools import *

try:
    from collections import OrderedDict
except ImportError:
    from ..util.ordered_dict import OrderedDict

inf = 1.0e20


# ----------------------------------------------------------------------
#  Configuration Class
# ----------------------------------------------------------------------

class Config(ordered_bunch):
    """ config = SU2.io.Config(filename="")
        
        Starts a config class, an extension of 
        ordered_bunch()
        
        use 1: initialize by reading config file
            config = SU2.io.Config('filename')
        use 2: initialize from dictionary or bunch
            config = SU2.io.Config(param_dict)
        use 3: initialize empty
            config = SU2.io.Config()
        
        Parameters can be accessed by item or attribute
        ie: config['MESH_FILENAME'] or config.MESH_FILENAME
        
        Methods:
            read()       - read from a config file
            write()      - write to a config file (requires existing file)
            dump()       - dump a raw config file
            unpack_dvs() - unpack a design vector 
            diff()       - returns the difference from another config
            dist()       - computes the distance from another config
    """    

    _filename = 'config.cfg'
    
    def __init__(self,*args,**kwarg):
        
        # look for filename in inputs
        if args and isinstance(args[0],str):
            filename = args[0]
            args = args[1:]
        elif kwarg.has_key('filename'):
            filename = kwarg['filename']
            del kwarg['filename']
        else:
            filename = ''
        
        # initialize ordered bunch
        super(Config,self).__init__(*args,**kwarg)
        
        # read config if it exists
        if filename:
            try:
                self.read(filename)
            except IOError:
                print 'Could not find config file: %s' % filename
	    except:
		print 'Unexpected error: ',sys.exc_info()[0]
		raise
        
        self._filename = filename
    
    def read(self,filename):
        """ reads from a config file """
        konfig = read_config(filename)
        self.update(konfig)
        
    def write(self,filename=''):
        """ updates an existing config file """
        if not filename: filename = self._filename
        assert os.path.exists(filename) , 'must write over an existing config file'
        write_config(filename,self)
        
    def dump(self,filename=''):
        """ dumps all items in the config bunch, without comments """
        if not filename: filename = self._filename
        dump_config(filename,self)
    
    def __getattr__(self,k):
        try:
            return super(Config,self).__getattr__(k)
        except AttributeError:
            raise AttributeError , 'Config parameter not found'
        
    def __getitem__(self,k):
        try:
            return super(Config,self).__getitem__(k)
        except KeyError:
            raise KeyError , 'Config parameter not found: %s' % k

    def unpack_dvs(self,dv_new,dv_old=None):
        """ updates config with design variable vectors
            will scale according to each DEFINITION_DV scale parameter
            
            Modifies:
                DV_KIND
                DV_MARKER
                DV_PARAM
                DV_VALUE_OLD
                DV_VALUE_NEW
            
            Inputs:
                dv_new - list or array of new dv values
                dv_old - optional, list or array of old dv values, defaults to zeros
                         
        """
        
        dv_new = copy.deepcopy(dv_new)
        dv_old = copy.deepcopy(dv_old)
        
        # handle unpacking cases
        def_dv = self['DEFINITION_DV']

        n_dv   = sum(def_dv['SIZE'])

        if not dv_old: dv_old = [0.0]*n_dv
        assert len(dv_new) == len(dv_old) , 'unexpected design vector length'
        
        # handle param
        param_dv = self['DV_PARAM']

        # apply scale
        dv_scales = def_dv['SCALE']

        k = 0
        for i, dv_scl in enumerate(dv_scales):
            for j in range(def_dv['SIZE'][i]):
                dv_new[k] = dv_new[k]*dv_scl;
                dv_old[k] = dv_old[k]*dv_scl;
                k = k + 1
        
        # Change the parameters of the design variables

        self['DV_KIND'] = def_dv['KIND']
        param_dv['PARAM'] = def_dv['PARAM']
        param_dv['FFDTAG'] = def_dv['FFDTAG']
        param_dv['SIZE']   = def_dv['SIZE']

        self.update({ 'DV_MARKER'        : def_dv['MARKER'][0] ,
                      'DV_VALUE_OLD'     : dv_old              ,
                      'DV_VALUE_NEW'     : dv_new              })
        
    def __eq__(self,konfig):
        return super(Config,self).__eq__(konfig)
    def __ne__(self,konfig):
        return super(Config,self).__ne__(konfig)
    
    
    def local_files(self):
        """ removes path prefix from all *_FILENAME params
        """
        for key,value in self.iteritems():
            if key.split('_')[-1] == 'FILENAME':
                self[key] = os.path.basename(value)    
    
    def diff(self,konfig):
        """ compares self to another config
            
            Inputs: 
                konfig - a second config
            
            Outputs:
                config_diff - a config containing only the differing 
                    keys, each with values of a list of the different 
                    config values.
                for example: 
                config_diff.MATH_PROBLEM = ['DIRECT','CONTINUOUS_ADJOINT']
                
        """
        
        keys = set([])
        keys.update( self.keys() )
        keys.update( konfig.keys() )
        
        konfig_diff = Config()
        
        for key in keys:
            value1 = self.get(key,None)
            value2 = konfig.get(key,None)
            if not value1 == value2:
                konfig_diff[key] = [value1,value2]
        
        return konfig_diff
    
    def dist(self,konfig,keys_check='ALL'):
        """ calculates a distance to another config
            
            Inputs: 
                konfig     - a second config
                keys_check - optional, a list of keys to check
            
            Outputs:
                distance   - a float
                
            Currently only works for DV_VALUE_NEW and DV_VALUE_OLD
            Returns a large value otherwise
                
        """        

        konfig_diff = self.diff(konfig)
        
        if keys_check == 'ALL':
            keys_check = konfig_diff.keys()
    
        distance = 0.0
        
        for key in keys_check:
            if konfig_diff.has_key(key):
                
                val1 = konfig_diff[key][0]
                val2 = konfig_diff[key][1]
                
                if key in ['DV_VALUE_NEW','DV_VALUE_OLD']:

                    val1 = np.array( val1 )
                    val2 = np.array( val2 )


                    this_diff = np.sqrt( np.sum( (val1-val2)**2 ) )
                
                elif key in ['MACH_NUMBER','REYNOLDS_NUMBER','AoA','BODY_FLAP_DEF','ELEVON_DEF','DV1','DV2','DV3','DV4','DV5','DV6','DRY_MASS','FUEL_MASS','THRUST','P_DYN_INF']:

                    val1 = float( val1 )
                    val2 = float( val2 )


                    this_diff = np.sqrt( np.sum( (val1-val2)**2 ) )

                else:
                    print 'Warning, unexpected config difference'
                    this_diff = inf
                    
                distance += this_diff
            
            #: if key different
        #: for each keys_check
        
        return distance
    
    def __repr__(self):
        #return '<Config> %s' % self._filename
        return self.__str__()
    
    def __str__(self):
        output = 'Config: %s' % self._filename
        for k,v in self.iteritems():
            output +=  '\n    %s= %s' % (k,v)
        return output
#: class Config







# -------------------------------------------------------------------
#  Get SU2 Configuration Parameters
# -------------------------------------------------------------------

def read_config(filename):
    """ reads a config file """
      
    # initialize output dictionary
    data_dict = OrderedDict()
    
    input_file = open(filename)
    
    # process each line
    while 1:
        # read the line
        line = input_file.readline()
        if not line:
            break
        
        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):
            continue
        # split across equals sign
        line = line.split("=",1)
        this_param = line[0].strip()
        this_value = line[1].strip()
        
        assert not data_dict.has_key(this_param) , ('Config file has multiple specifications of %s' % this_param )
        for case in switch(this_param):
            
            # otherwise
            # string parameters
            if case():
                data_dict[this_param] = this_value
                break              
            
            #: if case DEFINITION_DV
                        
        #: for case
        
    #: for line
 
    return data_dict
    
#: def read_config()



# -------------------------------------------------------------------
#  Set SU2 Configuration Parameters
# -------------------------------------------------------------------

def write_config(filename,param_dict):
    """ updates an existing config file """
    
    temp_filename = "temp.cfg"
    shutil.copy(filename,temp_filename)
    output_file = open(filename,"w")

    # break pointers
    param_dict = copy.deepcopy(param_dict)
    
    for raw_line in open(temp_filename):
        # remove line returns
        line = raw_line.strip('\r\n')
        
        # make sure it has useful data
        if not "=" in line:
            output_file.write(raw_line)
            continue
        
        # split across equals sign
        line = line.split("=")
        this_param = line[0].strip()
        old_value  = line[1].strip()
        
        # skip if parameter unwanted
        if not param_dict.has_key(this_param):
            output_file.write(raw_line)
            continue
        
        # start writing parameter
        new_value = param_dict[this_param] 
        output_file.write(this_param + "= ")
        
        # handle parameter types
        for case in switch(this_param):  
            
            # default, assume string, integer or unformatted float 
            if case():
                output_file.write('%s' % new_value)
                break                         
                
        #: for case
        
        # remove from param dictionary
        del param_dict[this_param]
        
        # next line
        output_file.write("\n")        
        
    #: for each line
    
    # check that all params were used
    for this_param in param_dict.keys():
        if not this_param in ['JOB_NUMBER']:
            print ( 'Warning: Parameter %s not found in config file and was not written' % (this_param) )
        
    output_file.close()
    os.remove( temp_filename )
    
#: def write_config()


def dump_config(filename,config):
    ''' dumps a raw config file with all options in config 
        and no comments
    '''
    
    # HACK - twl
    if config.has_key('DV_VALUE_NEW'):
        config.DV_VALUE = config.DV_VALUE_NEW
        
    config_file = open(filename,'w')
    # write dummy file
    for key in config.keys():
        config_file.write( '%s= 0 \n' % key )
    config_file.close()
    # dump data
    write_config(filename,config)    

