#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import subprocess
from ..io import Config
from ..util import which, mesh2tri, tri2mesh

# ------------------------------------------------------------
#  Setup
# ------------------------------------------------------------

BIN_PATH = os.environ['BIN_PATH'] 
sys.path.append( BIN_PATH )

SU2_RUN = os.environ['SU2_RUN'] 
sys.path.append( SU2_RUN )

SPACE_RUN = os.environ['SPACE_RUN']

MISSION_ANALYSIS_RUN = os.environ['MISSION_ANALYSIS_RUN'] 
sys.path.append( MISSION_ANALYSIS_RUN )

# SU2 suite run command template
base_Command = os.path.join(SU2_RUN,'%s')

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
#  SPACE Suite Interface Functions
# ------------------------------------------------------------

def MIS(config):
    """ run MIS
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)

    Command = MISSION_ANALYSIS_RUN + '/optimizer/bin/mission_analysis > optimizer.log'
    os.system(Command)
    #subprocess.call(command, shell=True)
    
    return

def CFD(config):
    """ run CFD
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)

    tempname = 'config_CFD.cfg'
    konfig.dump(tempname)

    os.system('SU2_CFD config_CFD.cfg > log_cfd.out') # ONLY WORKS WITH A SERIAL COMPILED SU2
    
    return

def SUR(config):
    """ run SUR
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    yams = open(konfig.FLUID_SURFACE + '.yams',"w")
    yams.write('Absolute\nGradation 1.3\nMinSize 0.01\nMaxsize 0.2\nGeomApp 0.001\n')
    yams.close()

    os.system('yams2 -f -O -1 -in ' + konfig.FLUID_SURFACE + '.mesh -out ' + konfig.FLUID_SURFACE + '.mesh > log_yams.out')

    return

def SYM(config):
    """ run SYM
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    mesh2tri(konfig)
    os.system('triangle -p ' + konfig.FLUID_SURFACE + '.poly > log_tri.out')
    tri2mesh(konfig)
    os.system('spider2 -O 6 -f -f64 -eps 0.0001 -in ' + konfig.FLUID_SURFACE + '.mesh ' + konfig.FARFIELD_FILENAME + ' ' + konfig.SYMMETRY_FILENAME + ' -out ' + konfig.BOUNDARY_FILENAME + ' >> log_spider.out 2>&1')

    return

def VOL(config):
    """ run VOL
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    adap_surf = open('adap.surf',"w")
    adap_surf.write('3    1\n\n')
    adap_surf.close()

    adap_source = open('adap.source',"w")
    adap_source.write('volume boundingbox  -2 20  -1  8  -5 7 h 0.2\n')
    adap_source.close()

    os.system('amg -novol -in ' + konfig.BOUNDARY_FILENAME + ' -hgrad 2.0 -out ' + konfig.BOUNDARY_FILENAME + ' > log_amg.out 2>&1')
    os.system('amg -novol -in ' + konfig.BOUNDARY_FILENAME + ' -hgrad 1.05 -source adap.source -out ' + konfig.BOUNDARY_FILENAME + ' >> log_amg.out 2>&1')
    os.system('amg -novol -in ' + konfig.BOUNDARY_FILENAME + ' -hgrad 1.05 -source adap.source -out ' + konfig.BOUNDARY_FILENAME + ' >> log_amg.out 2>&1')
    os.system('ghs3d -O 1 -in ' + konfig.BOUNDARY_FILENAME + ' -out ' + konfig.FLUID_VOLUME + '.meshb > log_ghs.out 2>&1')
    os.system('amg -in ' + konfig.FLUID_VOLUME + '.meshb -hgrad 1.05 -source adap.source -out ' + konfig.FLUID_VOLUME + '.meshb >> log_amg.out 2>&1')

    os.system('meshutils -O 3 -in ' + konfig.FLUID_VOLUME + '.meshb -out ' + konfig.FLUID_VOLUME + ' > log_meshutil.out')

    return

def MRG(config):
    """ run MRG
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    os.system('spider2 -O 6 -f -f64 -eps 0.0001 -in ' + konfig.FLUID_SURFACE + '.mesh ' + SPACE_RUN + '/SPACE/util/back/back.mesh -out ' + konfig.FLUID_SURFACE + '.mesh > log_mrg.out 2>&1')
    os.system('spider2 -O 1 -in ' + konfig.FLUID_SURFACE + '.mesh -Tri -Ref -from 1:1:10000 -to 1 -out ' + konfig.FLUID_SURFACE + '.mesh -f -f64 > log_spider.out 2>&1')

    return

def INT(config):
    """ run INT
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    os.system('mshint ' + konfig.FLUID_SURFACE + '.mesh ' + konfig.STRUCT + '_surface.mesh > log_mshint.out')
    
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

