#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import subprocess
from ..io import Config
from ..util import which

# ------------------------------------------------------------
#  Setup
# ------------------------------------------------------------

BIN_PATH = os.environ['BIN_PATH'] 
sys.path.append( BIN_PATH )

# SU2 suite run command template
base_Command = os.path.join(BIN_PATH,'%s')

# check for slurm
slurm_job = os.environ.has_key('SLURM_JOBID')

#check for tacc
tacc_job = os.environ.has_key('TACC_PUBLIC_MACHINE')

# set mpi command
if slurm_job:
    mpi_Command = 'srun -n %i %s'
    if tacc_job:
        mpi_Command = 'ibrun -o 0 -n %i %s'
elif not which('mpirun') is None:
    mpi_Command = 'mpirun -n %i %s'
elif not which('mpiexec') is None:
    mpi_Command = 'mpiexec -n %i %s'
else:
    mpi_Command = ''
    
from .. import EvaluationFailure, DivergenceFailure
return_code_map = {
    1 : EvaluationFailure ,
    2 : DivergenceFailure ,
}
    
# ------------------------------------------------------------
#  SU2 Suite Interface Functions
# ------------------------------------------------------------

def GHS(config):
    """ run GHS
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    tempname = 'config_GHS.cfg'
    konfig.dump(tempname)
    
    # must run with rank 1
    #processes = konfig['NUMBER_PART']
    #processes = min([1,processes])
    processes = 1

    in_file = konfig.FLUID_BOUNDARY_FILENAME.split('.')
    in_file = in_file[0] + '_updated.' + in_file[1]
    
    the_Command = 'ghs3d -O 1 -in ' + in_file + ' -out ' + konfig.FLUID_VOLUME + '.meshb'
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    os.system('meshutils -O 3 -in ' + konfig.FLUID_VOLUME + '.meshb -out ' + konfig.FLUID_VOLUME + ' > log_meshutil.out')

    #os.remove(tempname)
    
    return

def INT(config):
    """ run GHS
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    os.system('mshint ' + konfig.FLUID_SURFACE + '.mesh ' + konfig.STRUCT_SURFACE + '.mesh > log_mshint.out')
    
    return

# ------------------------------------------------------------
#  Helper functions
# ------------------------------------------------------------

def build_command( the_Command , processes=0 ):
    """ builds an mpi command for given number of processes """
    the_Command = base_Command % the_Command
    if processes >= 1:
        if not mpi_Command:
            raise RuntimeError , 'could not find an mpi interface'
        the_Command = mpi_Command % (processes,the_Command)
    return the_Command

def run_command( Command ):
    """ runs os command with subprocess
        checks for errors from command
    """
    
    sys.stdout.flush()
    
    proc = subprocess.Popen( Command, shell=True    ,
                             stdout=sys.stdout      , 
                             stderr=subprocess.PIPE  )
    return_code = proc.wait()
    message = proc.stderr.read()
    
    if return_code < 0:
        message = "SPACE process was terminated by signal '%s'\n%s" % (-return_code,message)
        raise SystemExit , message
    elif return_code > 0:
        message = "Path = %s\nCommand = %s\nSPACE process returned error '%s'\n%s" % (os.path.abspath(','),Command,return_code,message)
        if return_code in return_code_map.keys():
            exception = return_code_map[return_code]
        else:
            exception = RuntimeError
        raise exception , message
    else:
        sys.stdout.write(message)
            
    return return_code

