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

    #os.system('SU2_CFD config_CFD.cfg > log_cfd.out') # ONLY WORKS WITH A SERIAL COMPILED SU2

    processes = int(konfig['NUMBER_PART'])
    the_Command = 'SU2_CFD ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    return

def SOL(config):
    """ run SOL
      partitions set by config.NUMBER_PART
    """
  
    konfig = copy.deepcopy(config)
    
    tempname = 'config_SOL.cfg'
    konfig.dump(tempname)

    processes = int(konfig['NUMBER_PART'])
    the_Command = 'SU2_SOL ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    return

def AMG(config):
    
    # adap_surf = open('adap.surf',"w")
    # adap_surf.write('\n')
    # adap_surf.close()

    adap_surf = open('adap.surf',"w")
    adap_surf.write('1    -1\n')
    adap_surf.write('2    1\n')
    adap_surf.write('3    1\n\n')
    adap_surf.close()

    adap_source = open('adap.source',"w")
    adap_source.write('\n')
    adap_source.close()
    
    Command = 'amg -in fluid_volume.meshb -sol mach.solb -source adap.source -p 2 -c 160000 -hgrad 1.3 -hmin 0.00001 -hmax 70 -out current.new.meshb -itp restart_flow.solb -back boundary_back.mesh -nordg > log_amg.out'
    
    os.system(Command);

    return

def SUR(config):
    """ run SUR
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    yams = open(konfig.FLUID_SURFACE + '_back.yams',"w")
    #yams.write('Absolute\nGradation 1.1\nMinSize 0.01\nMaxsize 0.2\nGeomApp 0.001\n')
    yams.write('Absolute\nGradation 1.3\nMinSize 0.01\nMaxsize 0.25\nGeomApp 0.005\n')
    yams.close()

    #os.system('yams2 -f -O -1 -in ' + konfig.FLUID_SURFACE + '_back.mesh -out ' + konfig.FLUID_SURFACE + '.mesh > log_yams.out')
    
    #os.system('yams2 -f -O 0 -in ' + konfig.FLUID_SURFACE + '.mesh -out ' + konfig.FLUID_SURFACE + '.mesh >> log_yams.out')
    shutil.copy(konfig.FLUID_SURFACE + '_back.mesh', konfig.FLUID_SURFACE + '.mesh')

    return

def SYM(config):
    """ run SYM
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    mesh2tri(konfig)
    os.system('triangle -p ' + konfig.FLUID_SURFACE_CUR + '.poly > log_tri.out')
    tri2mesh(konfig)
    os.system('spider2 -O 6 -f -f64 -eps 0.0001 -in ' + konfig.FLUID_SURFACE_CUR + '.mesh ' + konfig.FARFIELD_FILENAME + ' ' + konfig.SYMMETRY_FILENAME + ' -out ' + konfig.BOUNDARY_FILENAME + ' >> log_spider.out 2>&1')

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
    #adap_source.write('volume boundingbox  -4 20  -1  8  -6 7 h 0.2\n')
    #adap_source.write('volume boundingbox  -9 25  -1  13  -11 12 h 0.8\n')z
    adap_source.write('volume boundingbox  -2.0 18.5  -1  2.5  -3.2  2.5 h 0.2\n')
    adap_source.write('volume boundingbox   2.0 15.5  -1  4.0  -3.2 -0.3 h 0.2\n')
    adap_source.write('volume boundingbox   8.5 15.5  -1  6.5  -3.2 -0.3 h 0.2\n')
    adap_source.close()

    os.system('amg -novol -in ' + konfig.BOUNDARY_FILENAME + ' -hgrad 2.0 -out ' + konfig.BOUNDARY_FILENAME + ' > log_amg.out 2>&1')
    os.system('amg -novol -in ' + konfig.BOUNDARY_FILENAME + ' -hgrad 1.1 -source adap.source -out ' + konfig.BOUNDARY_FILENAME + ' >> log_amg.out 2>&1')
    os.system('amg -novol -in ' + konfig.BOUNDARY_FILENAME + ' -hgrad 1.1 -source adap.source -out ' + konfig.BOUNDARY_FILENAME + ' >> log_amg.out 2>&1')
    os.system('amg -novol -in ' + konfig.BOUNDARY_FILENAME + ' -hgrad 1.1 -source adap.source -out ' + konfig.BOUNDARY_FILENAME + ' >> log_amg.out 2>&1')
    os.system('ghs3d -O 1 -in ' + konfig.BOUNDARY_FILENAME + ' -out ' + konfig.FLUID_VOLUME + '_euler.meshb > log_ghs.out 2>&1')
    # -grad 1.03 
    os.system('amg -in ' + konfig.FLUID_VOLUME + '_euler.meshb -hgrad 1.1 -source adap.source -out ' + konfig.FLUID_VOLUME + '_euler.meshb >> log_amg.out 2>&1')

    return

def BLG(config):
    """ run BLG
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    bloom = open(konfig.FLUID_VOLUME + '_euler.bloom',"w")

    bloom.write('BLReference\n1\n1\n\n')
    bloom.write('BLSymmetryReference\n1\n3\n\n')
    bloom.write('BLImprintReference\n1\n3\n\n')

    bloom.write('InitialSpacing\n' + konfig.INITIAL_SPACING + '\n\n')
    bloom.write('GrowthRate\n1.3\n\n')

    bloom.write('BLThickness\n0.2\n\n')
    #bloom.write('NumberOfLayers\n50\n\n')

    #bloom.write('TurbulentBoundaryLayer\n\n')
    bloom.write('MixedElementsBoundaryLayer\n\n')

    bloom.write('MeshDeformationReference\n1\n0\n\n')

    bloom.close()

    os.system('bloom -in ' + konfig.FLUID_VOLUME + '_euler.meshb -out ' + konfig.FLUID_VOLUME + ' > log_bloom.out 2>&1')

    return

def MRG(config):
    """ run MRG
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)
    
    #os.system('spider2 -O 6 -f -f64 -eps 0.0001 -in ' + konfig.FLUID_SURFACE + '_back.mesh ' + SPACE_RUN + '/SPACE/util/back/back.mesh -out ' + konfig.FLUID_SURFACE + '_back.mesh > log_mrg.out 2>&1')
    os.system('spider2 -O 1 -in ' + konfig.FLUID_SURFACE + '_back.mesh -Tri -Ref -from 1:1:10000 -to 1 -out ' + konfig.FLUID_SURFACE + '_back.mesh -f -f64 > log_spider.out 2>&1')

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

