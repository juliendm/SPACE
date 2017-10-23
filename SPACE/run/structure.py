#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import sys, os, copy

sys.path.append(os.environ['STRUCTURE_RUN'])

import numpy, argparse
from mpi4py import MPI
from baseclasses import *
from tacs import *

#from pyOpt import *
from pyoptsparse import *

# from .. import io   as spaceio
# from .. import util as spaceutil
# from interface import INT       as SPACE_INT

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE import io   as spaceio
from SPACE import util as spaceutil
from SPACE.run.interface import INT       as SPACE_INT


# ----------------------------------------------------------------------
#  Structure Simulation
# ----------------------------------------------------------------------

def structure(config):

    # local copy
    konfig = copy.deepcopy(config)

    # Safety Factors

    safetyFactor_thrust = 1.0 ##################################
    safetyFactor_inertial = 1.0 #2.0 ##################################
    safetyFactor_non_inertial = 1.0 #2.0 ##################################



    # Fixed by Trajectory

    fuel_percentage = 0.8
    pdyn_inf = 12169.185366 # float(konfig.P_DYN_INF)
    nx = 0.831404 # float(konfig.ACCELERATION_X)
    ny = 0.0 # float(konfig.ACCELERATION_Y)
    nz = -1.795355 # float(konfig.ACCELERATION_Z)
    half_thrust = 162844.38 # float(konfig.HALF_THRUST)

    loadFactor = numpy.sqrt(nx*nx+ny*ny+nz*nz)
    loadAngle = numpy.arccos(-nz/loadFactor)*180.0/numpy.pi # Angle with vector [0,0,-1]


    # Degree Of Freedom such that Sum M = 0 (Sum F = 0 via iteration with the Trajectory code)
    thrust_angle = -16.0


    # Modify inputs
    #loadFactor = ...
    #loadAngle = 30.0 # 24.8482002334

    # Update Values

    # nx = loadFactor*numpy.sin(loadAngle*numpy.pi/180.0)    
    # ny = 0.0
    # nz = -loadFactor*numpy.cos(loadAngle*numpy.pi/180.0)


    # NEED TO CHECK CONVENTIONS HERE !!!!!!!!!!!!!!!!!!!!!
    gravityVector = -9.81 * numpy.array([-nx,ny,-nz])/loadFactor # Change of Frame: to Structure Frame


    # Initial Guess

    ini_half_mass_guess = 20000 # kg


    # Compute initial Load

    spaceutil.surf2sol(konfig)
    SPACE_INT(konfig)
    load = spaceutil.Load(konfig, loadFactor, gravityVector, pdyn_inf, half_thrust, thrust_angle, fuel_percentage, safetyFactor_thrust, safetyFactor_inertial, safetyFactor_non_inertial)
    load.update(ini_half_mass_guess)




    #computeNastran(config, load)

    computeTacs(config, load)



    # # get history and objectives

    # history      = spaceio.read_history( history_filename )
    # outputs      = spaceutil.ordered_bunch()
    # outputs.MASS = history[obj_name][-1]

    # # info out

    # info = spaceio.State()
    # info.FUNCTIONS.update(outputs)

    # return info

    return 0

#: def structure()




def computeNastran(config, load):

    ini_thickness = 0.01
    min_thickness = 0.0016
    max_thickness = 3.0

    # Write bdf Nastran
    # -----------------

    bdf_nastran = open(config.STRUCT + '_nastran.bdf','w')

    bdf_nastran.write('SOL 200\n')
    bdf_nastran.write('TIME 600\n')
    bdf_nastran.write('CEND\n')

    bdf_nastran.write('ANALYSIS = STATICS\n')
    bdf_nastran.write('DISPLACEMENT = ALL\n')
    bdf_nastran.write('FORCE = ALL\n')
    bdf_nastran.write('STRESS = ALL\n')

    bdf_nastran.write('SPC = 1\n')
    bdf_nastran.write('LOAD = 1\n')
    bdf_nastran.write('DESOBJ = 1\n')
    bdf_nastran.write('DESSUB = 1\n')

    bdf_nastran.write('BEGIN BULK\n')

    # Copy Mesh
    bdf = open(config.STRUCT + '.bdf')
    for line in bdf:
        if (line.strip() != 'BEGIN BULK' and line.strip() != 'END BULK'):
            bdf_nastran.write(line)
    bdf.close()

    # Material Properties
    write_line(bdf_nastran,'MAT1',r=8)
    write_line(bdf_nastran,'1',r=8)
    write_line(bdf_nastran,config.MATERIAL_YOUNG_MODULUS,r=16)
    write_line(bdf_nastran,config.MATERIAL_POISSON_RATIO,r=8)
    write_line(bdf_nastran,config.MATERIAL_DENSITY,r=16)
    bdf_nastran.write('\n')

    # Shell
    for desc in load._descriptions.keys():
        write_line(bdf_nastran,'PSHELL',r=8)
        write_line(bdf_nastran,'%d' % (load._descriptions[desc]),r=8)
        write_line(bdf_nastran,'1',r=8)
        #write_line(bdf_nastran,('%.6f' % x_dvs[descriptions[desc]-1])[1:],r=8)
        write_line(bdf_nastran,('%.6f' % ini_thickness),r=8)
        write_line(bdf_nastran,'1',r=16)
        write_line(bdf_nastran,'1',r=8)
        bdf_nastran.write('\n')

    # Loads
    # -----

    write_line(bdf_nastran,'LOAD',r=8)
    write_line(bdf_nastran,'1',r=8)
    write_line(bdf_nastran,'1.0',r=8)
    write_line(bdf_nastran,'1.0',r=8)
    write_line(bdf_nastran,'2',r=8)
    write_line(bdf_nastran,'1.0',r=8)
    write_line(bdf_nastran,'3',r=8)
    bdf_nastran.write('\n')

    # Force

    for iPoint_bdf in range(load._nPoint_bdf):
        if ((load._load_bdf[iPoint_bdf][0] != 0.0) or (load._load_bdf[iPoint_bdf][1] != 0.0) or (load._load_bdf[iPoint_bdf][2] != 0.0)):
            write_line(bdf_nastran,'FORCE*  ')

            write_line(bdf_nastran,'2',l=16)
            write_line(bdf_nastran,'%d' % (iPoint_bdf+1),l=16)
            write_line(bdf_nastran,'0',l=16)
            write_line(bdf_nastran,'1.0',l=16)

            write_line(bdf_nastran,'*F')
            write_line(bdf_nastran,'%d' % (iPoint_bdf+1),l=6)            
            write_line(bdf_nastran,'\n')
            write_line(bdf_nastran,'*F')
            write_line(bdf_nastran,'%d' % (iPoint_bdf+1),l=6)

            write_line(bdf_nastran,'%.8E' % load._load_bdf[iPoint_bdf][0],l=16)
            write_line(bdf_nastran,'%.8E' % load._load_bdf[iPoint_bdf][1],l=16)
            write_line(bdf_nastran,'%.8E' % load._load_bdf[iPoint_bdf][2],l=16)

            write_line(bdf_nastran,'\n')

    # Acceleration

    write_line(bdf_nastran,'GRAV',r=8)
    write_line(bdf_nastran,'3',r=8)
    write_line(bdf_nastran,'0',r=8)
    write_line(bdf_nastran,'9.81',r=8)
    write_line(bdf_nastran,'%.3f' % (load._gravity_vector[0]/9.81*load._loadFactor*load._safetyFactor_inertial),r=8)
    write_line(bdf_nastran,'%.3f' % (load._gravity_vector[1]/9.81*load._loadFactor*load._safetyFactor_inertial),r=8)
    write_line(bdf_nastran,'%.3f' % (load._gravity_vector[2]/9.81*load._loadFactor*load._safetyFactor_inertial),r=8)
    bdf_nastran.write('\n')

    # Design Variables
    # ----------------

    for desc in load._descriptions.keys():
        write_line(bdf_nastran,'DESVAR',r=8)
        write_line(bdf_nastran,'%d' % (load._descriptions[desc]),r=8)
        write_line(bdf_nastran,'V_' + str(load._descriptions[desc]),r=8)
        write_line(bdf_nastran,str(ini_thickness),r=8)
        write_line(bdf_nastran,str(min_thickness),r=8)
        write_line(bdf_nastran,str(max_thickness),r=8)
        write_line(bdf_nastran,'1.0',r=8)
        bdf_nastran.write('\n')

    for desc in load._descriptions.keys():
        write_line(bdf_nastran,'DVPREL1',r=8)
        write_line(bdf_nastran,'%d' % (load._descriptions[desc]),r=8)
        write_line(bdf_nastran,'PSHELL',r=8)
        write_line(bdf_nastran,'%d' % (load._descriptions[desc]),r=8)
        write_line(bdf_nastran,'T',r=8)
        # write_line(bdf_nastran,str(min_thickness),r=8)
        # write_line(bdf_nastran,str(max_thickness),r=8)
        # write_line(bdf_nastran,str(0.0),r=8)
        bdf_nastran.write('\n')
        write_line(bdf_nastran,'',r=8)
        write_line(bdf_nastran,'%d' % (load._descriptions[desc]),r=8)
        write_line(bdf_nastran,str(1.0),r=8)
        bdf_nastran.write('\n')

    # Objective
    # ---------

    dresp_id = 0

    write_line(bdf_nastran,'DRESP1',r=8)
    dresp_id += 1
    write_line(bdf_nastran,str(dresp_id),r=8)
    write_line(bdf_nastran,'W',r=8)
    write_line(bdf_nastran,'WEIGHT',r=8)
    bdf_nastran.write('\n')

    # Constraints
    # -----------

    con_id = 0

    con_id += 1
    write_line(bdf_nastran,'DCONADD',r=8)
    write_line(bdf_nastran,str(con_id),r=8)
    write_line(bdf_nastran,str(con_id+1),r=8)
    # write_line(bdf_nastran,str(con_id+2),r=8)
    bdf_nastran.write('\n')

    # Stresses
    con_id += 1

    var_stress = [9,17]
    for var in var_stress:
        write_line(bdf_nastran,'DRESP1',r=8)
        dresp_id += 1
        write_line(bdf_nastran,str(dresp_id),r=8)
        if (var == 9):
            write_line(bdf_nastran,'VMSTR1',r=8)
        elif (var == 17):
            write_line(bdf_nastran,'VMSTR2',r=8)
        write_line(bdf_nastran,'STRESS',r=8)
        write_line(bdf_nastran,'ELEM',r=16)
        write_line(bdf_nastran,str(var),r=8)
        write_line(bdf_nastran,'',r=8)
        check = 7
        for iElem_bdf in range(load._nElem_bdf):
            if (check == 8):
                check = 0
                bdf_nastran.write('\n')
                write_line(bdf_nastran,'',r=8)
            check +=1
            write_line(bdf_nastran,str(iElem_bdf+1),r=8)
        if (check != 0):
            bdf_nastran.write('\n')

    write_line(bdf_nastran,'DRESP2',r=8)
    dresp_id += 1
    write_line(bdf_nastran,str(dresp_id),r=8)
    write_line(bdf_nastran,'MAXVM',r=8)
    write_line(bdf_nastran,str(con_id),r=8)
    bdf_nastran.write('\n')
    write_line(bdf_nastran,'',r=8)
    write_line(bdf_nastran,'DRESP1',r=8)
    for iDim in range(len(var_stress)):
        write_line(bdf_nastran,str(dresp_id-len(var_stress)+iDim),r=8)
    bdf_nastran.write('\n')

    write_line(bdf_nastran,'DEQATN',r=8)
    write_line(bdf_nastran,str(con_id),r=8)
    bdf_nastran.write('MAXVM(VMSTR1,VMSTR2)=MAX(VMSTR1,VMSTR2)')
    bdf_nastran.write('\n')

    write_line(bdf_nastran,'DCONSTR',r=8)
    write_line(bdf_nastran,str(con_id),r=8)
    write_line(bdf_nastran,str(dresp_id),r=8)
    write_line(bdf_nastran,'1e-09',r=8)

    write_line(bdf_nastran,'200e6',r=8)   # 324e6
#    write_line(bdf_nastran,config.MATERIAL_YIELD_STRENGTH,r=8)

    bdf_nastran.write('\n')

    # Optimization
    # ------------

    write_line(bdf_nastran,'DOPTPRM',r=8) 
    write_line(bdf_nastran,'DESMAX',r=8)  
    write_line(bdf_nastran,'50',r=8)     
    # write_line(bdf_nastran,'CT',r=8) 
    # write_line(bdf_nastran,'-.03',r=8)  
    # write_line(bdf_nastran,'CTMIN',r=8)    
    # write_line(bdf_nastran,'.003',r=8)     
    # write_line(bdf_nastran,'CONV1',r=8)  
    # write_line(bdf_nastran,'1.-5',r=8)    
    # bdf_nastran.write('\n')
    # write_line(bdf_nastran,'',r=8) 
    # write_line(bdf_nastran,'CONV2',r=8) 
    # write_line(bdf_nastran,'1.-20',r=8)    
    # write_line(bdf_nastran,'CONVDV',r=8) 
    # write_line(bdf_nastran,'1.-6',r=8)  
    # write_line(bdf_nastran,'CONVPR',r=8)    
    # write_line(bdf_nastran,'1.-5',r=8)
    bdf_nastran.write('\n')

# Default

# CONV1 0.001       Objective relative change 2 iterations
# CONV2 1.0E-20     Objective absolute change 2 iterations
# CONVDV 0.0001     Design variables
# CONVPR 0.001      Properties
# CT -0.03          Constraint tolerance (active if value greater than CT)
# CTMIN 0.003                            (violated if value greater than CTMIN)

# DELB 0.0001       Relative finite difference move parameter
# DELP 0.2          Properties: Fractional change allowed
# DELX 0.5          Design variables: Fractional change allowed

# DESMAX 5

# DPMAX 0.5
# DPMIN 0.01        Properties: Minimum move limit

# DXMAX 1.0
# DXMIN 0.05        Design variables: Minimum move limit

# PENAL 0.0         Improve perfo if starting design is infeasible when e.g. 2.0

    bdf_nastran.write('ENDDATA\n')

    bdf_nastran.close()

#: def computeNastran()






















def computeTacs(config, load):

    gcomm = comm = MPI.COMM_WORLD

    # Material properties

    material_rho = float(config.MATERIAL_DENSITY)
    material_E = float(config.MATERIAL_YOUNG_MODULUS)
    material_ys = 80e6 #float(config.MATERIAL_YIELD_STRENGTH) / 1.6 #1.5 #2.0 ##################################
    material_nu = float(config.MATERIAL_POISSON_RATIO)
    kcorr = 5.0/6.0

    t = 0.01
    tMin = 0.0016 # 0.0016
    tMax = 3.0 # 0.020

    KSWeight = 80.0
    evalFuncs = ['mass','ks0','mf0']

    boost_factor = 1.0 #12.0 ##################################

    SPs = [StructProblem('lc0', loadFactor=load._loadFactor*boost_factor, loadFile=config.LOAD_FILENAME, evalFuncs=evalFuncs)]
    #SPs = [StructProblem('lc0', loadFactor=loadFactor*safetyFactor_inertial, loadFile=konfig.LOAD_FILENAME, evalFuncs=['mass','ks0','ks1','ks2'])]
    numLoadCases = len(SPs)

    # Create Solver

    structOptions = {'transferSize':0.5, 'transferGaussOrder':3}
    FEASolver = pytacs.pyTACS(config.STRUCT + '.bdf', comm=comm, options=structOptions)

    # Add Design Variables

    ndv, corresp = addDVGroups(FEASolver)

    def conCallBack(dvNum, compDescripts, userDescript, specialDVs, **kargs):
        print dvNum
        con = constitutive.isoFSDTStiffness(material_rho, material_E, material_nu, kcorr, material_ys, t, dvNum, tMin, tMax)
        # if userDescript in ['JUNCTIONS','FRAMES','LONGERONS','SPARS']:
        #     con = constitutive.isoFSDTStiffness(material_rho, material_E, material_nu, kcorr, material_ys, 0.05, dvNum, 0.0016, 0.05)
        scale = [100.0]
        return con, scale

    FEASolver.createTACSAssembler(conCallBack)

    assert ndv == FEASolver.getNumDesignVars()

    # Add Functions

    # Mass Functions
    FEASolver.addFunction('mass', functions.StructuralMass)

    # KS Functions
    ks0 = FEASolver.addFunction('ks0', functions.AverageKSFailure, KSWeight=KSWeight, loadFactor=1.0)
    #ks0 = FEASolver.addFunction('ks0', functions.AverageKSFailure, KSWeight=KSWeight, include=RIBS+SPARS+FRAMES+LONGERONS+WING_BOX, loadFactor=1.0)
    #ks1 = FEASolver.addFunction('ks1', functions.AverageKSFailure, KSWeight=KSWeight, include=SKIN_U+STRINGERS_U, loadFactor=1.0)
    #ks2 = FEASolver.addFunction('ks2', functions.AverageKSFailure, KSWeight=KSWeight, include=SKIN_L+STRINGERS_L, loadFactor=1.0)

    #ksef0 = FEASolver.addFunction('ksef0', functions.KSElementFailure, KSWeight=KSWeight)
    #ksf0 = FEASolver.addFunction('ksf0', functions.KSFailure, KSWeight=KSWeight)
    mf0 = FEASolver.addFunction('mf0', functions.AverageMaxFailure)

    #ad0 = FEASolver.addFunction('ad0', functions.AggregateDisplacement)


    # # # NO LOAD FACTOR VIA TACS ANYMORE, LOAD FACTOR ON STRUCTURE IS DONE VIA LOAD FILE # # #
    # Load Factor

    # FEASolver.setOption('gravityVector',load._gravity_vector.tolist())
    FEASolver.setOption('gravityVector',(load._gravity_vector*load._loadFactor*load._safetyFactor_inertial).tolist())

    for i in range(numLoadCases):
        FEASolver.addInertialLoad(SPs[i])

    history_filename = 'history_structure.dat'
    history_iteration = {'val':0}

    # Process current state

    x_final = numpy.zeros(FEASolver.getNumDesignVars())
    FEASolver.structure.getDesignVars(x_final)
    load.postprocess(x_final, corresp)

    # Objective

    def obj(x):
        '''Evaluate the objective and constraints'''
        funcs = {}
        FEASolver.setDesignVars(x)
        for i in range(numLoadCases):

            #############################################
            # load.postprocess(x['struct'], corresp) # Update load._structure_mass and load._additional_mass
            # load.update(load._half_structure_mass+load._half_additional_mass)
            # SPs[i].loadFile = config.LOAD_FILENAME # Reset loadFile to read it again
            #############################################

            FEASolver(SPs[i])
            FEASolver.evalFunctions(SPs[i], funcs)
        if comm.rank == 0:
            history_file = open(history_filename,'a')
            history_file.write('%d' % history_iteration['val'])
            for key in funcs.keys():
                history_file.write(',%.16f' % funcs[key])
            history_file.write('\n')
            history_file.close()
            history_iteration['val'] += 1
        return funcs, False

    # Sensitivies

    def sens(x, funcs):
        '''Evaluate the objective and constraint sensitivities'''
        funcsSens = {}
        for i in range(numLoadCases):
            FEASolver.evalFunctionsSens(SPs[i], funcsSens)
        return funcsSens, False


    # Set up the optimization problem

    history_file = open(history_filename,'w')
    history_file.write('VARIABLES = "Iteration"')

    optProb = Optimization('Mass min', obj)
    obj_name = 'lc0_mass'
    optProb.addObj(obj_name)
    history_file.write(',"%s"' % obj_name)
    FEASolver.addVariablesPyOpt(optProb)
    for i in range(numLoadCases):
        for j in xrange(1):
            con_name = '%s_ks%d'% (SPs[i].name, j)
            optProb.addCon(con_name, lower=1.0, upper=1.0)
            history_file.write(',"%s"' % con_name)
            con_name = '%s_mf%d'% (SPs[i].name, j)
            #optProb.addCon(con_name, lower=1.0, upper=1.0)
            history_file.write(',"%s"' % con_name)
    history_file.write('\n')
    history_file.close()

    if comm.rank == 0:
        print optProb
    optProb.printSparsity()

    opt = OPT('snopt',options={
        'Major feasibility tolerance':1e-6,
        'Major optimality tolerance':1e-6,
        'Minor feasibility tolerance':1e-6,
        'Iterations limit':100000,
        'Major iterations limit':3000,
        'Minor iterations limit':500,
        'Major step limit':2.0})

    # Solve

    sol = opt(optProb, sens=sens) #NULL result without error in PyObject_Call

    # Write Files

    write_files(config, FEASolver, SPs[0], corresp, load)

#: def computeTacs()



def addDVGroups(FEASolver):

    # SKIN

    SKIN_FUSE_U = ['FUSE:TOP','FUSE:LFT','FUSE_R','CTAIL:LOW','CTAIL_T::1']
    SKIN_FUSE_L = ['FUSE:BOT','FLAP:UPP','FLAP:LOW','FUSE_F'] # 'FLAP_T',
    SKIN_WING_U = ['LWING:UPP','LWING_T::0']
    SKIN_WING_L = ['LWING:LOW','LWING_T::1']
    SKINS = SKIN_FUSE_U + SKIN_FUSE_L + SKIN_WING_U + SKIN_WING_L + ['MSKINC:a','MSKINC:b']

    # JUNCTIONS

    JUNCTIONS = ['FLAP_FUSE','LWING_FUSE','CTAIL_FUSE']

    # MEMBERS

    FRAMES = ['MFRAME:00','MFRAME:01','MFRAME:02','MFRAME:03','MFRAME:04','MFRAME:05','MFRAME:06','MFRAME:07','MFRAME:08','MFRAME:09',
        'MFRAME:10','MFRAME:11','MFRAME:12']
    LONGERONS = ['MLONG:02:2','MLONG:00:3','MLONG:01:3','MLONG:02:3','MLONG:00:4','MLONG:01:4']
    RIBS = ['MRIBF:00','MRIBF:01','MRIBF:02','MRIBF:03','MRIBF:04','MRIBF:05','MRIBF:06','MRIBF:07',
        'MRIBV:00','MRIBV:01','MRIBV:02','MRIBV:03','MRIBV:04','MRIBV:05','MRIBV:06','MRIBV:07','MRIBV:08','MRIBV:09',
        'MRIBW:00','MRIBW:01','MRIBW:02','MRIBW:03','MRIBW:04','MRIBW:05']
    SPARS = ['MSPARF:00','MSPARF:01', # 'MSPARF:02','MSPARF:03',
        'MSPARV:00','MSPARV:01','MSPARV:02',
        'MSPARC:00','MSPARC:06',
        'MSPARW:00','MSPARW:02','MSPARW:08'] # 'MSPARW:09'
    STRINGERS = ['MSTRINGC:01','MSTRINGC:02','MSTRINGC:03','MSTRINGC:04','MSTRINGC:05',
        'MSTRINGW:01','MSTRINGW:03','MSTRINGW:04','MSTRINGW:05','MSTRINGW:06','MSTRINGW:07']
    MEMBERS = FRAMES + LONGERONS + RIBS + SPARS + STRINGERS

    assert len(FEASolver.selectCompIDs(include=SKINS+JUNCTIONS+MEMBERS)[0]) == FEASolver.nComp

    corresp = [-1 for index in range(FEASolver.nComp)]
    ndv = 0;

    # SKIN_IDS = FEASolver.selectCompIDs(include=SKINS)[0]
    # for i in range(len(SKIN_IDS)):
    #     dv_name = "SKIN_" + str(i)
    #     FEASolver.addDVGroup(dv_name, include = SKIN_IDS[i])
    #     ndv = ndv+1;
    #     corresp[SKIN_IDS[i]] = ndv;

    for i in range(len(SKINS)):
        dv_name = SKINS[i]
        FEASolver.addDVGroup(dv_name, include = SKINS[i])
        ndv = ndv+1;
        SKIN_IDS_I = FEASolver.selectCompIDs(include=SKINS[i])[0]
        for k in range(len(SKIN_IDS_I)):
            corresp[SKIN_IDS_I[k]] = ndv;

    # JUNCTION_IDS = FEASolver.selectCompIDs(include=JUNCTIONS)[0]
    # for i in range(len(JUNCTION_IDS)):
    #     dv_name = "JUNCTION_" + str(i)
    #     FEASolver.addDVGroup(dv_name, include = JUNCTION_IDS[i])
    #     ndv = ndv+1;
    #     corresp[JUNCTION_IDS[i]] = ndv;

    for i in range(len(JUNCTIONS)):
        dv_name = JUNCTIONS[i]
        FEASolver.addDVGroup(dv_name, include = JUNCTIONS[i])
        ndv = ndv+1;
        JUNCTIONS_IDS_I = FEASolver.selectCompIDs(include=JUNCTIONS[i])[0]
        for k in range(len(JUNCTIONS_IDS_I)):
            corresp[JUNCTIONS_IDS_I[k]] = ndv;

    for i in range(len(MEMBERS)):
        dv_name = MEMBERS[i]
        FEASolver.addDVGroup(dv_name, include = MEMBERS[i])
        ndv = ndv+1;
        MEMBERS_IDS_I = FEASolver.selectCompIDs(include=MEMBERS[i])[0]
        for k in range(len(MEMBERS_IDS_I)):
            corresp[MEMBERS_IDS_I[k]] = ndv;

    return ndv, corresp

    # ncoms = FEASolver.nComp
    # for i in range(0,ncoms):
    #     dv_name = 'stru_'+str(i)
    #     FEASolver.addDVGroup(dv_name, include = i)


def write_files(config, FEASolver, SP, corresp, load):

    FEASolver.writeBDFForces(SP, "visualize_forces.bdf")
    FEASolver.writeMeshDisplacements(SP, "struct_tacs.sol")
    FEASolver.writeSolution()

    x_final = numpy.zeros(FEASolver.getNumDesignVars())
    FEASolver.structure.getDesignVars(x_final)

    thickness_point = [0.0 for iPoint_bdf in range(load._nPoint_bdf)]
    thickness_point_count = [0 for iPoint_bdf in range(load._nPoint_bdf)]

    for iElem_bdf in range(load._nElem_bdf):
        for iNode in range(load._nNode):
            thickness_point[load._elem_bdf[iElem_bdf][iNode]-1] += x_final[corresp[load._elem_tag_bdf[iElem_bdf]-1]-1]
            thickness_point_count[load._elem_bdf[iElem_bdf][iNode]-1] += 1
    for iPoint_bdf in range(load._nPoint_bdf):
        thickness_point[iPoint_bdf] = thickness_point[iPoint_bdf]/thickness_point_count[iPoint_bdf]

    write_sol_1('thicknesses.sol',thickness_point)
    write_tecplot('thicknesses.dat',load._coord_bdf,load._elem_bdf,thickness_point)

    dvs_file = 'x_final.dat'
    dvs = open(dvs_file,'w')
    for i in range(len(x_final)):
        dvs.write('%f\n' % x_final[i])
    dvs.close()

    load.postprocess(x_final, corresp)

def write_sol_1(sol_file,solution):
 
    sol = open(sol_file,'w')
    sol.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nSolAtVertices\n' + str(len(solution)) + '\n1 1\n')
    for i_point in range(len(solution)):
        sol.write(str(solution[i_point]) + "\n")
    sol.write('\nEnd\n')
    sol.close()

#: def write_sol_1()

def write_sol_3(sol_file,solution):
 
    sol = open(sol_file,'w')
    sol.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nSolAtVertices\n' + str(len(solution)) + '\n3 1 1 1\n')
    for i_point in range(len(solution)):
        sol.write(str(solution[i_point][0]) + " " + str(solution[i_point][1]) + " " + str(solution[i_point][2]) + "\n")
    sol.write('\nEnd\n')
    sol.close()

#: def write_sol_3()

def write_tecplot(sol_file,coord,elem,thickness):
 
    nPoint = len(coord)
    nElem = len(elem)

    file = open(sol_file,'w')
    file.write('TITLE = "Visualization of the surface solution"\n')
    file.write('VARIABLES = "x""y""z""thickness"\n')
    file.write('ZONE NODES= %d, ELEMENTS= %d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n' % (nPoint, nElem))
    for iPoint in range(nPoint):
        file.write(str(coord[iPoint][0]) + " " + str(coord[iPoint][1]) + " " + str(coord[iPoint][2]) + " ")
        # file.write(str(force[iPoint][0]) + " " + str(force[iPoint][1]) + " " + str(force[iPoint][2]) + " ")
        # file.write(str(disp[iPoint][0]) + " " + str(disp[iPoint][1]) + " " + str(disp[iPoint][2]) + " ")
        # file.write(str(rotation[iPoint][0]) + " " + str(rotation[iPoint][1]) + " " + str(rotation[iPoint][2]) + " ")
        file.write(str(thickness[iPoint]) + "\n")
    for iElem in range(nElem):
        file.write(str(elem[iElem][0]) + " " + str(elem[iElem][1])  + " " + str(elem[iElem][2]) + " " + str(elem[iElem][2]) + "\n")
    file.close()

#: def write_tecplot()


def write_line(file,line,l=0,r=0):
    if l is not 0:
        n = l - len(line)
        for i in range(n):
            line = ' ' + line
    if r is not 0:
        n = r - len(line)
        for i in range(n):
            line = line + ' '
    file.write(line)


def isFloat(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    
    config = spaceio.Config('../config_DSN.cfg')
    structure(config)


