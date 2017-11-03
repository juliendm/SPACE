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

    safetyFactor_thrust = 1.0
    safetyFactor_inertial = 1.0
    safetyFactor_non_inertial = 1.0

    # Initial Guess

    ini_half_mass_guess = 20000 # kg


    loads = []

    # LOAD CASE 1

    fuel_percentage = 0.8
    pdyn_inf = 12169.185366 # float(konfig.P_DYN_INF)
    nx = 0.831404 # float(konfig.ACCELERATION_X)
    ny = 0.0 # float(konfig.ACCELERATION_Y)
    nz = -1.795355 # float(konfig.ACCELERATION_Z)
    half_thrust = 162844.38 # float(konfig.HALF_THRUST)

    loadFactor = numpy.sqrt(nx*nx+ny*ny+nz*nz)
    thrust_angle = 0.0 # Degree Of Freedom such that Sum M = 0 (Sum F = 0 via iteration with the Trajectory code)
    gravity_vector = -9.81 * numpy.array([-nx,ny,-nz])/loadFactor # Change of Frame: to Structure Frame


    load_filename = 'load_1.dat'

    spaceutil.surf2sol(konfig)
    SPACE_INT(konfig)
    load = spaceutil.Load(konfig, load_filename, loadFactor, gravity_vector, pdyn_inf, half_thrust, thrust_angle, fuel_percentage, safetyFactor_thrust, safetyFactor_inertial, safetyFactor_non_inertial)
    load.update(ini_half_mass_guess)

    loads.append(load)


    # LOAD CASE 2

    fuel_percentage = 1.0
    pdyn_inf = 8458.985536
    nx = 5.0
    ny = 0.0
    nz = 0.0
    half_thrust = 162844.38
    mach = 4.072702
    reynolds = 873499.095851
    aoa = 22.956569

    load_filename = 'load_2.dat'

    spaceutil.surf2sol(konfig)
    SPACE_INT(konfig)
    load = spaceutil.Load(konfig, load_filename, loadFactor, gravity_vector, pdyn_inf, half_thrust, thrust_angle, fuel_percentage, safetyFactor_thrust, safetyFactor_inertial, safetyFactor_non_inertial)
    load.update(ini_half_mass_guess)

    loads.append(load)


    # Run

    computeTacs(konfig, loads)




    # # Modify inputs
    # loadAngle = numpy.arccos(-nz/loadFactor)*180.0/numpy.pi # Angle with vector [0,0,-1]
    # loadFactor = ...
    # loadAngle = 30.0 # 24.8482002334
    # nx = loadFactor*numpy.sin(loadAngle*numpy.pi/180.0)    
    # ny = 0.0
    # nz = -loadFactor*numpy.cos(loadAngle*numpy.pi/180.0)



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























def computeTacs(config, loads):

    gcomm = comm = MPI.COMM_WORLD

    # Material properties

    material_rho = float(config.MATERIAL_DENSITY)
    material_E = float(config.MATERIAL_YOUNG_MODULUS)
    material_ys = float(config.MATERIAL_YIELD_STRENGTH) * 0.3333 #2.7 #80e6 #80e6 ##################################
    material_nu = float(config.MATERIAL_POISSON_RATIO)
    kcorr = 5.0/6.0

    t = 0.01
    t_skin = 0.002

    tMin = 0.0016

    tMax = 0.1
    tMax_skin = 0.05

    KSWeight = 80.0

    #evalFuncs = ['mass','ks0','ks1','ks2']
    evalFuncs = ['mass','ksf1','ksf2','ksf3']
    #evalFuncs = ['mass','ks1','ks2','ad0']
    #evalFuncs = ['mass','ad0']





    # Create Solver
    structOptions = {'transferSize':0.5, 'transferGaussOrder':3}
    FEASolver = pytacs.pyTACS(config.STRUCT + '.bdf', comm=comm, options=structOptions)

    # Load Cases
    SPs = []
    numLoadCases = len(loads)

    # Load File
    for i in range(numLoadCases):
        SPs.append(StructProblem('lc%d' % i, loadFactor=loads[i]._loadFactor*loads[i]._safetyFactor_inertial, loadFile=loads[i]._load_filename, evalFuncs=evalFuncs))

    # Add Design Variables
    ndv, corresp, SKINS, JUNCTIONS, MEMBERS = addDVGroups(FEASolver, loads[0])
    def conCallBack(dvNum, compDescripts, userDescript, specialDVs, **kargs):
        if 'SKIN' in userDescript:
            con = constitutive.isoFSDTStiffness(material_rho, material_E, material_nu, kcorr, material_ys, t_skin, dvNum, tMin, tMax_skin)
        else:
            con = constitutive.isoFSDTStiffness(material_rho, material_E, material_nu, kcorr, material_ys, t, dvNum, tMin, tMax)
        scale = [100.0]
        return con, scale
    FEASolver.createTACSAssembler(conCallBack)
    assert ndv == FEASolver.getNumDesignVars()

    # Inertial Forces
    for i in range(numLoadCases):
        FEASolver.setOption('gravityVector',loads[i]._gravity_vector.tolist())
        FEASolver.addInertialLoad(SPs[i])

    # Add Functions

    # Mass Functions
    FEASolver.addFunction('mass', functions.StructuralMass)

    # Contraint Functions
    if 'ksf0' in evalFuncs: ksf0 = FEASolver.addFunction('ksf0', functions.AverageKSFailure, KSWeight=KSWeight, loadFactor=1.0)
    if 'ksf1' in evalFuncs: ksf1 = FEASolver.addFunction('ksf1', functions.AverageKSFailure, KSWeight=KSWeight, include=SKINS, loadFactor=1.0)
    if 'ksf2' in evalFuncs: ksf2 = FEASolver.addFunction('ksf2', functions.AverageKSFailure, KSWeight=KSWeight, include=JUNCTIONS, loadFactor=1.0)
    if 'ksf3' in evalFuncs: ksf3 = FEASolver.addFunction('ksf3', functions.AverageKSFailure, KSWeight=KSWeight, include=MEMBERS, loadFactor=1.0)

    if 'ksb0' in evalFuncs: ksb0 = FEASolver.addFunction('ksb0', functions.AverageKSBuckling, KSWeight=KSWeight, loadFactor=1.0)
    if 'ksb1' in evalFuncs: ksb1 = FEASolver.addFunction('ksb1', functions.AverageKSBuckling, KSWeight=KSWeight, include=SKINS, loadFactor=1.0)
    if 'ksb2' in evalFuncs: ksb2 = FEASolver.addFunction('ksb2', functions.AverageKSBuckling, KSWeight=KSWeight, include=JUNCTIONS, loadFactor=1.0)
    if 'ksb3' in evalFuncs: ksb3 = FEASolver.addFunction('ksb3', functions.AverageKSBuckling, KSWeight=KSWeight, include=MEMBERS, loadFactor=1.0)

    if 'ad0' in evalFuncs: ad0 = FEASolver.addFunction('ad0', functions.AggregateDisplacement, KSWeight=100, loadFactor=1.0)
    if 'ad1' in evalFuncs: ad1 = FEASolver.addFunction('ad1', functions.AggregateDisplacement, KSWeight=100, include=SKINS, loadFactor=1.0)
    if 'ad2' in evalFuncs: ad2 = FEASolver.addFunction('ad2', functions.AggregateDisplacement, KSWeight=100, include=JUNCTIONS, loadFactor=1.0)
    if 'ad3' in evalFuncs: ad3 = FEASolver.addFunction('ad3', functions.AggregateDisplacement, KSWeight=100, include=MEMBERS, loadFactor=1.0)


    # File output
    history_filename = 'history_structure.dat'
    history_iteration = {'val':0}

    # Process current state
    x_final = numpy.zeros(FEASolver.getNumDesignVars())
    FEASolver.structure.getDesignVars(x_final)
    for i in range(numLoadCases):
        loads[i].postprocess(x_final, corresp)

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
            for i in range(numLoadCases):
                for key in evalFuncs:
                    history_file.write(',%.16f' % funcs['%s_%s'% (SPs[i].name, key)])
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

    optProb = Optimization('Mass min', obj)

    history_file = open(history_filename,'w')
    history_file.write('VARIABLES = "Iteration"')


    for i in range(numLoadCases):

        for name in evalFuncs:

            if 'mass' in name:
                obj_name = '%s_%s' % (SPs[i].name, name)
                if obj_name == 'lc0_mass':
                    optProb.addObj(obj_name)
                history_file.write(',"%s"' % obj_name)

            if 'ksf' in name:
                con_name = '%s_%s' % (SPs[i].name, name)
                optProb.addCon(con_name, upper=1.0)
                history_file.write(',"%s"' % con_name)

            if 'ksb' in name:
                con_name = '%s_%s' % (SPs[i].name, name)
                optProb.addCon(con_name, upper=1.0)
                history_file.write(',"%s"' % con_name)

            if 'ad' in name:
                con_name = '%s_%s' % (SPs[i].name, name)
                optProb.addCon(con_name, upper=0.03) # no more than 3 cm disp
                history_file.write(',"%s"' % con_name)

    history_file.write('\n')
    history_file.close()

    # Add Variables to optProb
    FEASolver.addVariablesPyOpt(optProb)


    if comm.rank == 0:
        print optProb
    optProb.printSparsity()

    opt = OPT('snopt',options={
        'Major feasibility tolerance':1e-6,
        'Major optimality tolerance':1e-6,
        'Minor feasibility tolerance':1e-6,
        'Iterations limit':100000,
        'Major iterations limit':500,
        'Minor iterations limit':500,
        'Major step limit':2.0})

    # Solve

    sol = opt(optProb, sens=sens)

    # Write Files

    print_tag = []
    for key_in in JUNCTIONS + MEMBERS:
        for key_desc in load._descriptions:
            if key_in.upper() in key_desc:
                print_tag.append(load._descriptions[key_desc])

    write_files(config, FEASolver, SPs[0], corresp, load, print_tag)

#: def computeTacs()

def addDVGroups(FEASolver, load):

    # SKIN

    SKIN_FUSE_U = ['FUSE:TOP','FUSE:LFT','FUSE_R','CTAIL:LOW','CTAIL_T::1']
    SKIN_FUSE_L = ['FUSE:BOT','FLAP:UPP','FLAP:LOW','FUSE_F'] # 'FLAP_T',
    SKIN_WING_U = ['LWING:UPP','LWING_T::0']
    SKIN_WING_L = ['LWING:LOW','LWING_T::1']
    SKINS = SKIN_FUSE_U + SKIN_FUSE_L + SKIN_WING_U + SKIN_WING_L

    # JUNCTIONS

    JUNCTIONS = ['FLAP_FUSE','LWING_FUSE','CTAIL_FUSE']

    # MEMBERS

    MEMBERS = []
    for name in load._descriptions.keys():
        parts = name.split(':')
        if parts[0] in ['MFRAME','MRIBF','MRIBV','MRIBW','MSPARF','MSPARV','MSPARC','MSPARW','MSTRINGC','MSTRINGW']:
            new_name = '%s:%s' % (parts[0],parts[1])
            if new_name not in MEMBERS:
                MEMBERS.append(new_name)
        elif parts[0] in ['MLONG','MSKINC','MSKINC']:
            new_name = '%s:%s:%s' % (parts[0],parts[1],parts[2])
            if new_name not in MEMBERS:
                MEMBERS.append(new_name)     

    # # MFRAME:xx
    # FRAMES = ['MFRAME:00','MFRAME:01','MFRAME:02','MFRAME:03','MFRAME:04','MFRAME:05','MFRAME:06','MFRAME:07','MFRAME:08','MFRAME:09',
    #     'MFRAME:10','MFRAME:11','MFRAME:12']
    # # MLONG:xx:x
    # LONGERONS = ['MLONG:02:2','MLONG:00:3','MLONG:01:3','MLONG:02:3','MLONG:00:4','MLONG:01:4']
    # # MRIBF:xx MRIBV:xx MRIBW:xx
    # RIBS = ['MRIBF:00','MRIBF:01','MRIBF:02','MRIBF:03','MRIBF:04','MRIBF:05','MRIBF:06','MRIBF:07',
    #     'MRIBV:00','MRIBV:01','MRIBV:02','MRIBV:03','MRIBV:04','MRIBV:05','MRIBV:06','MRIBV:07','MRIBV:08','MRIBV:09',
    #     'MRIBW:00','MRIBW:01','MRIBW:02','MRIBW:03','MRIBW:04','MRIBW:05']
    # # MSPARF:xx MSPARV:xx MSPARC:xx MSPARW:xx
    # SPARS = ['MSPARF:00','MSPARF:01',
    #     'MSPARV:00','MSPARV:01','MSPARV:02',
    #     'MSPARC:03','MSPARC:09',
    #     'MSPARW:00','MSPARW:03','MSPARW:09']
    # # MSTRINGC:xx MSTRINGW:xx
    # STRINGERS = ['MSTRINGC:04','MSTRINGC:05','MSTRINGC:06','MSTRINGC:07','MSTRINGC:08',
    #     'MSTRINGW:01','MSTRINGW:02','MSTRINGW:04','MSTRINGW:05','MSTRINGW:06','MSTRINGW:07','MSTRINGW:08']
    # # MSKINC:a:xx MSKINC:b:xx
    # SKIN_BOX = ['MSKINC:a:03','MSKINC:a:04','MSKINC:a:05','MSKINC:a:06','MSKINC:a:07','MSKINC:a:08','MSKINC:b:03','MSKINC:b:04','MSKINC:b:05','MSKINC:b:06','MSKINC:b:07','MSKINC:b:08',]
    # MEMBERS = FRAMES + LONGERONS + RIBS + SPARS + STRINGERS + SKIN_BOX

    assert len(FEASolver.selectCompIDs(include=SKINS+JUNCTIONS+MEMBERS)[0]) == FEASolver.nComp

    corresp = [-1 for index in range(FEASolver.nComp)]
    ndv = 0;

    group_skin = False
    group_junction = True

    if group_skin:
        for i in range(len(SKINS)):
            dv_name = "SKIN_" + str(i) # SKINS[i]
            FEASolver.addDVGroup(dv_name, include = SKINS[i])
            ndv = ndv+1;
            SKIN_IDS_I = FEASolver.selectCompIDs(include=SKINS[i])[0]
            for k in range(len(SKIN_IDS_I)):
                corresp[SKIN_IDS_I[k]] = ndv;
    else:
        SKIN_IDS = FEASolver.selectCompIDs(include=SKINS)[0]
        for i in range(len(SKIN_IDS)):
            dv_name = "SKIN_" + str(i)
            FEASolver.addDVGroup(dv_name, include = SKIN_IDS[i])
            ndv = ndv+1;
            corresp[SKIN_IDS[i]] = ndv;

    if group_junction:
        for i in range(len(JUNCTIONS)):
            dv_name = "JUNCTION_" + str(i) # JUNCTIONS[i]
            FEASolver.addDVGroup(dv_name, include = JUNCTIONS[i])
            ndv = ndv+1;
            JUNCTIONS_IDS_I = FEASolver.selectCompIDs(include=JUNCTIONS[i])[0]
            for k in range(len(JUNCTIONS_IDS_I)):
                corresp[JUNCTIONS_IDS_I[k]] = ndv;
    else:
        JUNCTION_IDS = FEASolver.selectCompIDs(include=JUNCTIONS)[0]
        for i in range(len(JUNCTION_IDS)):
            dv_name = "JUNCTION_" + str(i)
            FEASolver.addDVGroup(dv_name, include = JUNCTION_IDS[i])
            ndv = ndv+1;
            corresp[JUNCTION_IDS[i]] = ndv;

    for i in range(len(MEMBERS)):
        dv_name = "MEMBER_" + str(i) # MEMBERS[i]
        FEASolver.addDVGroup(dv_name, include = MEMBERS[i])
        ndv = ndv+1;
        MEMBERS_IDS_I = FEASolver.selectCompIDs(include=MEMBERS[i])[0]
        for k in range(len(MEMBERS_IDS_I)):
            corresp[MEMBERS_IDS_I[k]] = ndv;

    return ndv, corresp, SKINS, JUNCTIONS, MEMBERS

    # ncoms = FEASolver.nComp
    # for i in range(0,ncoms):
    #     dv_name = 'stru_'+str(i)
    #     FEASolver.addDVGroup(dv_name, include = i)

#: def addDVGroups()

def computeSimpleTacs():

    bdf_file = '../beam.bdf'
    load_filename = 'load_distribution.txt'

    # Material properties

    material_rho = 2780.0
    material_E = 72.0E9
    material_nu = 0.33
    material_ys = 324.0E6
    kcorr = 5.0/6.0

    load_factor = 3.0

    F = 100000.0
    b = 10.0
    h = 5.0
    L = 110.0
    w = 9.81*load_factor*b*h*material_rho
    beam_I = b*h**3.0/12.0
    delta_F = F*L**3.0/(3.0*material_E*beam_I) # Load
    delta_G = w*L**4.0/(8.0*material_E*beam_I) # Gravity
    print "Disp F: ", delta_F
    print "Disp G: ", delta_G

    # READ

    bdf = open(bdf_file)
    coord_bdf = []
    elem_bdf = []
    elem_tag_bdf = []
    descriptions = {}
    for line in bdf:
        data = line.split()
        if (line[0]=="$" and len(data) == 3):
            descriptions[data[2].strip().split('/')[0].upper()] = int(data[1])
        if (line[0]=="G" and len(data) == 4):
            vec = [float(data[2]), float(data[3].split('*')[0])]
        elif (line[0]=="*" and len(data) == 2):
            vec.append(float(data[1]))
            coord_bdf.append(vec)
        elif (line[0]=="C" and len(data) == 8):
            elem_bdf.append([int(data[3]), int(data[4]), int(data[5]), int(data[6])])
            elem_tag_bdf.append(int(data[2]))
    bdf.close()

    nPoint_bdf = len(coord_bdf)
    nElem_bdf = len(elem_bdf)

    print nPoint_bdf
    print nElem_bdf

    # LOAD

    nDim = 3
    load_bdf = [[0.0 for iDim in range(nDim)] for iPoint_bdf in range(nPoint_bdf)]

    #load_bdf[50] = [0.0, 0.0, F]

    # WRITE

    load_file = open(load_filename,"w")
    load_file.write(str(nPoint_bdf) + " " + str(nElem_bdf) + "\n")
    for iPoint_bdf in range(nPoint_bdf):
        load_file.write(str(coord_bdf[iPoint_bdf][0]) + " " + str(coord_bdf[iPoint_bdf][1]) + " " + str(coord_bdf[iPoint_bdf][2]) + " " + str(load_bdf[iPoint_bdf][0]) + " " + str(load_bdf[iPoint_bdf][1]) + " " + str(load_bdf[iPoint_bdf][2]) + "\n")
    for iElem_bdf in range(nElem_bdf):
        load_file.write(str(elem_bdf[iElem_bdf][0]-1) + " " + str(elem_bdf[iElem_bdf][1]-1)  + " " + str(elem_bdf[iElem_bdf][2]-1) + " " + str(elem_bdf[iElem_bdf][3]-1) + "\n")
    load_file.close()

    gcomm = comm = MPI.COMM_WORLD



    t = h
    tMin = 0.0016
    tMax = 2.0*h

    SPs = [StructProblem('lc0', loadFactor=load_factor, loadFile=load_filename)]
    numLoadCases = len(SPs)

    # Create Solver

    structOptions = {'transferSize':0.5, 'transferGaussOrder':3}
    FEASolver = pytacs.pyTACS(bdf_file, comm=comm, options=structOptions)

    # Add Design Variables

    ncoms = FEASolver.nComp
    for i in range(0,ncoms):
        dv_name = 'stru_'+str(i)
        FEASolver.addDVGroup(dv_name, include = i)

    def conCallBack(dvNum, compDescripts, userDescript, specialDVs, **kargs):
        con = constitutive.isoFSDTStiffness(material_rho, material_E, material_nu, kcorr, material_ys, t, dvNum, tMin, tMax)
        scale = [100.0]
        return con, scale

    FEASolver.createTACSAssembler(conCallBack)

    # NO load_factor HERE !!!!!!!!!!!!! Already taken into account in StructProblem definition
    FEASolver.setOption('gravityVector',(9.81*numpy.array([0.0,0.0,-1.0])).tolist()) 
    for i in range(numLoadCases):
        FEASolver.addInertialLoad(SPs[i])

    for i in range(numLoadCases):
        FEASolver(SPs[i])
        disp = FEASolver.writeMeshDisplacements(SPs[i], "beam_disp.sol")
        force = FEASolver.writeMeshForces(SPs[i], "beam_force.sol")
        thickness = [0.0 for iPoint_bdf in range(nPoint_bdf)]

        write_tecplot("beam.dat",coord_bdf,elem_bdf,force,disp,thickness)

#: def computeSimpleTacs()


def computeNastran(config, load):

    ini_thickness = 0.01
    min_thickness = 0.0016
    max_thickness = 3.0

    thickness_tag = numpy.loadtxt('thickness_final.dat')

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
        string = line.strip()
        if (string[0] == '$' or string[0] == 'G' or string[0] == '*'):
            bdf_nastran.write(line)
    bdf.close()

    for iElem_bdf in range(load._nElem_bdf):

        coord_0 = load._coord_bdf[load._elem_bdf[iElem_bdf][0]-1]
        coord_1 = load._coord_bdf[load._elem_bdf[iElem_bdf][1]-1]
        coord_2 = load._coord_bdf[load._elem_bdf[iElem_bdf][2]-1]
        coord_3 = load._coord_bdf[load._elem_bdf[iElem_bdf][3]-1]

        if (index_max_angle(coord_0,coord_1,coord_2,coord_3) in [0,2]):

            write_line(bdf_nastran,'CTRIA3',r=8)
            write_line(bdf_nastran,str(2*iElem_bdf+1),r=8)
            write_line(bdf_nastran,str(load._elem_tag_bdf[iElem_bdf]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][0]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][1]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][2]),r=8)
            bdf_nastran.write('\n')

            write_line(bdf_nastran,'CTRIA3',r=8)
            write_line(bdf_nastran,str(2*iElem_bdf+2),r=8)
            write_line(bdf_nastran,str(load._elem_tag_bdf[iElem_bdf]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][2]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][3]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][0]),r=8)
            bdf_nastran.write('\n')

        else:

            write_line(bdf_nastran,'CTRIA3',r=8)
            write_line(bdf_nastran,str(2*iElem_bdf+1),r=8)
            write_line(bdf_nastran,str(load._elem_tag_bdf[iElem_bdf]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][1]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][2]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][3]),r=8)
            bdf_nastran.write('\n')

            write_line(bdf_nastran,'CTRIA3',r=8)
            write_line(bdf_nastran,str(2*iElem_bdf+2),r=8)
            write_line(bdf_nastran,str(load._elem_tag_bdf[iElem_bdf]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][3]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][0]),r=8)
            write_line(bdf_nastran,str(load._elem_bdf[iElem_bdf][1]),r=8)
            bdf_nastran.write('\n')

    # Copy Mesh
    bdf = open(config.STRUCT + '.bdf')
    for line in bdf:
        string = line.strip()
        if (string[0] == 'S'):
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
        write_line(bdf_nastran,('%.6f' % thickness_tag[load._descriptions[desc]-1]),r=8)
        #write_line(bdf_nastran,('%.5f' % ini_thickness),r=8)
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


    optimise = False

    if optimise:

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



def write_files(config, FEASolver, SP, corresp, load, print_tag):

    FEASolver.writeBDFForces(SP, "visualize_forces.bdf")
    FEASolver.writeSolution()

    x_final = numpy.zeros(FEASolver.getNumDesignVars())
    FEASolver.structure.getDesignVars(x_final)

    thickness_point = [0.0 for iPoint_bdf in range(load._nPoint_bdf)]
    thickness_point_count = [0 for iPoint_bdf in range(load._nPoint_bdf)]

    current_tag = load._elem_tag_bdf[0]
    for iElem_bdf in range(load._nElem_bdf):
        new_tag = load._elem_tag_bdf[iElem_bdf]
        if new_tag == current_tag:
            for iNode in range(load._nNode):
                thickness_point[load._elem_bdf[iElem_bdf][iNode]-1] += x_final[corresp[load._elem_tag_bdf[iElem_bdf]-1]-1]
                thickness_point_count[load._elem_bdf[iElem_bdf][iNode]-1] += 1
        else:
            for iNode in range(load._nNode):
                thickness_point[load._elem_bdf[iElem_bdf][iNode]-1] = x_final[corresp[load._elem_tag_bdf[iElem_bdf]-1]-1]
                thickness_point_count[load._elem_bdf[iElem_bdf][iNode]-1] = 1
        current_tag = new_tag

    for iPoint_bdf in range(load._nPoint_bdf):
        thickness_point[iPoint_bdf] = thickness_point[iPoint_bdf]/thickness_point_count[iPoint_bdf]

    thickness_elem = [0.0 for iElem_bdf in range(load._nElem_bdf)]
    for iElem_bdf in range(load._nElem_bdf):
        thickness_elem[iElem_bdf] = x_final[corresp[load._elem_tag_bdf[iElem_bdf]-1]-1]

    thickness_tag = [0.0 for iTag_bdf in range(max(load._elem_tag_bdf))]
    for iTag_bdf in range(max(load._elem_tag_bdf)):
        thickness_tag[iTag_bdf] = x_final[corresp[iTag_bdf]-1]

    write_sol_1('struct_thickness.sol',thickness_point)
    disp = FEASolver.writeMeshDisplacements(SP, "struct_disp.sol")
    force = FEASolver.writeMeshForces(SP, "struct_force.sol")

    write_tecplot('surface_struct_interior.dat',load._coord_bdf,load._elem_bdf,load._elem_tag_bdf,print_tag,force,disp,thickness_point)
    write_tecplot('surface_struct.dat',load._coord_bdf,load._elem_bdf,load._elem_tag_bdf,load._elem_tag_bdf,force,disp,thickness_point)

    thickness_file = open('thickness_final.dat','w')
    for iTag_bdf in range(max(load._elem_tag_bdf)):
        thickness_file.write('%f\n' % thickness_tag[iTag_bdf])
    thickness_file.close()

    dvs_file = open('x_final.dat','w')
    for iDesVar in range(len(x_final)):
        dvs_file.write('%f\n' % x_final[iDesVar])
    dvs_file.close()

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

def write_tecplot(sol_file,coord,elem,elem_tag,print_tag,force,disp,thickness):
 
    nPoint = len(coord)
    nElem = len(elem)

    nElem_print = 0
    for iElem in range(nElem):
        if elem_tag[iElem] in print_tag:
            nElem_print += 1

    file = open(sol_file,'w')
    file.write('TITLE = "Visualization of the surface solution"\n')
    file.write('VARIABLES = "x""y""z""f_x""f_y""f_z""d_x""d_y""d_z""thickness"\n')
    file.write('ZONE NODES= %d, ELEMENTS= %d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n' % (nPoint, nElem_print))
    for iPoint in range(nPoint):
        file.write(str(coord[iPoint][0]) + " " + str(coord[iPoint][1]) + " " + str(coord[iPoint][2]) + " ")
        file.write(str(force[iPoint][0]) + " " + str(force[iPoint][1]) + " " + str(force[iPoint][2]) + " ")
        file.write(str(disp[iPoint][0]) + " " + str(disp[iPoint][1]) + " " + str(disp[iPoint][2]) + " ")
        # file.write(str(rotation[iPoint][0]) + " " + str(rotation[iPoint][1]) + " " + str(rotation[iPoint][2]) + " ")
        file.write(str(thickness[iPoint]) + "\n")
    for iElem in range(nElem):
        if elem_tag[iElem] in print_tag:
            file.write(str(elem[iElem][0]) + " " + str(elem[iElem][1])  + " " + str(elem[iElem][2]) + " " + str(elem[iElem][3]) +  "\n")
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

def index_max_angle(coord_0,coord_1,coord_2,coord_3):

    D01 = ( (coord_0[0]-coord_1[0])**2.0 + (coord_0[1]-coord_1[1])**2.0 + (coord_0[2]-coord_1[2])**2.0 )**0.5 # 01
    D02 = ( (coord_0[0]-coord_2[0])**2.0 + (coord_0[1]-coord_2[1])**2.0 + (coord_0[2]-coord_2[2])**2.0 )**0.5 # 02
    D03 = ( (coord_0[0]-coord_3[0])**2.0 + (coord_0[1]-coord_3[1])**2.0 + (coord_0[2]-coord_3[2])**2.0 )**0.5 # 03
    D12 = ( (coord_1[0]-coord_2[0])**2.0 + (coord_1[1]-coord_2[1])**2.0 + (coord_1[2]-coord_2[2])**2.0 )**0.5 # 12
    D13 = ( (coord_1[0]-coord_3[0])**2.0 + (coord_1[1]-coord_3[1])**2.0 + (coord_1[2]-coord_3[2])**2.0 )**0.5 # 13
    D23 = ( (coord_2[0]-coord_3[0])**2.0 + (coord_2[1]-coord_3[1])**2.0 + (coord_2[2]-coord_3[2])**2.0 )**0.5 # 23
    D10 = D01; D21 = D12; D31 = D13; D30 = D03; D32 = D23; D20 = D02

    # angle 01 - 03
    D1 = D01; D2 = D03; D3 = D13
    angle_0 = numpy.arccos((D1**2.0 + D2**2.0 - D3**2.0) / (2.0 * D1 * D2))

    # angle 12 - 10
    D1 = D12; D2 = D10; D3 = D20
    angle_1 = numpy.arccos((D1**2.0 + D2**2.0 - D3**2.0) / (2.0 * D1 * D2))

    # angle 23 - 21
    D1 = D23; D2 = D21; D3 = D31
    angle_2 = numpy.arccos((D1**2.0 + D2**2.0 - D3**2.0) / (2.0 * D1 * D2))

    # angle 30 - 32
    D1 = D30; D2 = D32; D3 = D02
    angle_3 = numpy.arccos((D1**2.0 + D2**2.0 - D3**2.0) / (2.0 * D1 * D2))

    angles = [angle_0, angle_1, angle_2, angle_3]

    return angles.index(max(angles))

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


