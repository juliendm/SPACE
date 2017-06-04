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
# from load      import load      as spaceload

sys.path.append(os.environ['SPACE_RUN'])
import SPACE
from SPACE import io   as spaceio
from SPACE import util as spaceutil
from SPACE.run.interface import INT       as SPACE_INT
from SPACE.run.load      import load      as spaceload


# ----------------------------------------------------------------------
#  Structure Simulation
# ----------------------------------------------------------------------

def structure(config):

    gcomm = comm = MPI.COMM_WORLD

    # local copy
    konfig = copy.deepcopy(config)

    # Load

    spaceutil.surf2sol(konfig)
    SPACE_INT(konfig)

    nx = float(konfig.ACCELERATION_X)
    ny = float(konfig.ACCELERATION_Y)
    nz = float(konfig.ACCELERATION_Z)

    loadFactor = numpy.sqrt(nx*nx+ny*ny+nz*nz)
    gravityVector = -9.81 * numpy.array([-nx,-nz,-ny])/numpy.sqrt(nx*nx+ny*ny+nz*nz) # Change of Frame: to Structure Frame

    info = spaceload(konfig, loadFactor, gravityVector)

    # Material properties

    material_rho = float(konfig.MATERIAL_DENSITY)
    material_E = float(konfig.MATERIAL_YOUNG_MODULUS)
    material_ys = 2e7 #float(konfig.MATERIAL_YIELD_STRENGTH)
    material_nu = float(konfig.MATERIAL_POISSON_RATIO)
    kcorr = 5.0/6.0

    t = 1.0
    tMin = 0.0016 # 0.0016
    tMax = 3.0 # 0.020

    KSWeight = 80.0
    evalFuncs = ['mass','ks0','mf0']
    SPs = [StructProblem('lc0', loadFactor=loadFactor, loadFile=konfig.LOAD_FILENAME, evalFuncs=evalFuncs)]
    #SPs = [StructProblem('lc0', loadFactor=loadFactor, loadFile=konfig.LOAD_FILENAME, evalFuncs=['mass','ks0','ks1','ks2'])]
    numLoadCases = len(SPs)

    # Create Solver

    structOptions = {'transferSize':0.5, 'transferGaussOrder':3}
    FEASolver = pytacs.pyTACS(konfig.STRUCT + '.bdf', comm=comm, options=structOptions)

    # Add Design Variables

    corresp = addDVGroups(FEASolver)

    def conCallBack(dvNum, compDescripts, userDescript, specialDVs, **kargs):
        con = constitutive.isoFSDTStiffness(material_rho, material_E, material_nu, kcorr, material_ys, t, dvNum, tMin, tMax)
        # if userDescript in ['JUNCTIONS','FRAMES','LONGERONS','SPARS']:
        #     con = constitutive.isoFSDTStiffness(material_rho, material_E, material_nu, kcorr, material_ys, 0.05, dvNum, 0.0016, 0.05)
        scale = [100.0]
        return con, scale

    FEASolver.createTACSAssembler(conCallBack)

    # Add Functions

    # Mass Functions
    FEASolver.addFunction('mass', functions.StructuralMass)

    # KS Functions
    ks0 = FEASolver.addFunction('ks0', functions.AverageKSFailure, KSWeight=KSWeight, loadFactor=1.0)
    #ks0 = FEASolver.addFunction('ks0', functions.AverageKSFailure, KSWeight=KSWeight, include=RIBS+SPARS+FRAMES+LONGERONS+WING_BOX, loadFactor=loadFactor)
    #ks1 = FEASolver.addFunction('ks1', functions.AverageKSFailure, KSWeight=KSWeight, include=SKIN_U+STRINGERS_U, loadFactor=loadFactor)
    #ks2 = FEASolver.addFunction('ks2', functions.AverageKSFailure, KSWeight=KSWeight, include=SKIN_L+STRINGERS_L, loadFactor=loadFactor)

    #ksef0 = FEASolver.addFunction('ksef0', functions.KSElementFailure, KSWeight=KSWeight)
    #ksf0 = FEASolver.addFunction('ksf0', functions.KSFailure, KSWeight=KSWeight)
    mf0 = FEASolver.addFunction('mf0', functions.AverageMaxFailure)

    #ad0 = FEASolver.addFunction('ad0', functions.AggregateDisplacement)

    # Load Factor

    FEASolver.setOption('gravityVector',gravityVector.tolist())
    for i in range(numLoadCases):
        FEASolver.addInertialLoad(SPs[i])

    history_filename = 'history_structure.dat'
    history_iteration = {'val':0}

    def obj(x):
        '''Evaluate the objective and constraints'''
        funcs = {}
        FEASolver.setDesignVars(x)
        for i in range(numLoadCases):
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
        'Major iterations limit':3,
        'Minor iterations limit':500,
        'Major step limit':2.0})
    sol = opt(optProb, sens=sens) #NULL result without error in PyObject_Call

    # Write Files

    write_files(config, FEASolver, SPs[0])

    # get history and objectives

    history      = spaceio.read_history( history_filename )
    outputs      = spaceutil.ordered_bunch()
    outputs.MASS = history[obj_name][-1]

    # info out

    info = spaceio.State(info)
    info.FUNCTIONS.update(outputs)

    return info

#: def structure()

def addDVGroups(FEASolver):

    # SKIN

    SKIN_FUSE_U = ['FUSE:TOP','FUSE:LFT','FUSE_R'] # 'CTAIL:LOW','CTAIL_T',
    SKIN_FUSE_L = ['FUSE:BOT','FLAP:UPP','FLAP:LOW','FLAP_T','FUSE_F']
    SKIN_WING_U = ['LWING:UPP','LWING_T::0']
    SKIN_WING_L = ['LWING:LOW','LWING_T::1']
    SKINS = SKIN_FUSE_U + SKIN_FUSE_L + SKIN_WING_U + SKIN_WING_L + ['MSKINC:a','MSKINC:b']

    # JUNCTIONS

    JUNCTIONS = ['FLAP_FUSE','LWING_FUSE'] # 'CTAIL_FUSE'

    # MEMBERS

    FRAMES = ['MFRAME:00','MFRAME:01','MFRAME:02','MFRAME:03','MFRAME:04','MFRAME:05','MFRAME:06','MFRAME:07','MFRAME:08','MFRAME:09',
        'MFRAME:10','MFRAME:11','MFRAME:12']
    LONGERONS = ['MLONG:02:2','MLONG:00:3','MLONG:01:3','MLONG:02:3','MLONG:00:4','MLONG:01:4']
    RIBS = ['MRIBF:00','MRIBF:01','MRIBF:02','MRIBF:03','MRIBF:04','MRIBF:05','MRIBF:06','MRIBF:07',
        # 'MRIBV:00','MRIBV:01','MRIBV:02','MRIBV:03','MRIBV:04','MRIBV:05','MRIBV:06','MRIBV:07','MRIBV:08','MRIBV:09',
        'MRIBW:00','MRIBW:01','MRIBW:02','MRIBW:03','MRIBW:04','MRIBW:05']
    SPARS = ['MSPARF:00','MSPARF:01', # 'MSPARF:02','MSPARF:03',
        # 'MSPARV:00','MSPARV:01',
        'MSPARC:00','MSPARC:04','MSPARC:05',
        'MSPARW:00','MSPARW:02','MSPARW:08'] # 'MSPARW:09'
    STRINGERS = ['MSTRINGC:01','MSTRINGC:02','MSTRINGC:03',
        'MSTRINGW:01','MSTRINGW:03','MSTRINGW:04','MSTRINGW:05','MSTRINGW:06','MSTRINGW:07']
    MEMBERS = FRAMES + LONGERONS + RIBS + SPARS + STRINGERS

    assert len(FEASolver.selectCompIDs(include=SKINS+JUNCTIONS+MEMBERS)[0]) == FEASolver.nComp

    corresp = [-1 for index in range(FEASolver.nComp)]
    dv = -1;

    SKIN_IDS = FEASolver.selectCompIDs(include=SKINS)[0]
    for i in range(len(SKIN_IDS)):
        dv_name = "SKIN_" + str(i)
        FEASolver.addDVGroup(dv_name, include = SKIN_IDS[i])
        dv = dv+1;
        corresp[SKIN_IDS[i]] = dv;

    JUNCTION_IDS = FEASolver.selectCompIDs(include=JUNCTIONS)[0]
    for i in range(len(JUNCTION_IDS)):
        dv_name = "JUNCTION_" + str(i)
        FEASolver.addDVGroup(dv_name, include = JUNCTION_IDS[i])
        dv = dv+1;
        corresp[JUNCTION_IDS[i]] = dv;

    for i in range(len(MEMBERS)):
        dv_name = MEMBERS[i]
        FEASolver.addDVGroup(dv_name, include = MEMBERS[i])
        dv = dv+1;
        MEMBERS_IDS_I = FEASolver.selectCompIDs(include=MEMBERS[i])[0]
        for k in range(len(MEMBERS_IDS_I)):
            corresp[MEMBERS_IDS_I[k]] = dv;

    return corresp

    # ncoms = FEASolver.nComp
    # for i in range(0,ncoms):
    #     dv_name = 'stru_'+str(i)
    #     FEASolver.addDVGroup(dv_name, include = i)


def write_files(config, FEASolver, SP):

    FEASolver.writeBDFForces(SP, "visualize_forces.bdf")
    FEASolver.writeMeshDisplacements(SP, "struct_tacs.sol")
    FEASolver.writeSolution()

    x_final = numpy.zeros(FEASolver.nComp)
    FEASolver.structure.getDesignVars(x_final)

    n_point_bdf = 0
    elem_bdf = []
    elem_tag_bdf = []
    bdf = open(config.STRUCT + '.bdf')
    for line in bdf:
        data = line.split()
        if (line[0]=="G" and len(data) == 6):
            n_point_bdf += 1
        elif (line[0]=="C" and len(data) == 7):
            elem_bdf.append([int(data[3]), int(data[4]), int(data[5]), int(data[6])])
            elem_tag_bdf.append(int(data[2]))
    bdf.close()
    n_elem_bdf = len(elem_bdf)
    thickness_point = [0.0 for i_point_bdf in range(n_point_bdf)]
    thickness_point_count = [0 for i_point_bdf in range(n_point_bdf)]
    for i_elem_bdf in range(n_elem_bdf):
        for i_node in range(4):
            thickness_point[elem_bdf[i_elem_bdf][i_node]-1] += x_final[elem_tag_bdf[i_elem_bdf]-1]
            thickness_point_count[elem_bdf[i_elem_bdf][i_node]-1] += 1
    for i_point_bdf in range(n_point_bdf):
        thickness_point[i_point_bdf] = thickness_point[i_point_bdf]/thickness_point_count[i_point_bdf]
    write_sol_1('thicknesses.sol',thickness_point)

    dvs_file = 'x_final.dat'
    dvs = open(dvs_file,'w')
    for i in range(len(x_final)):
        dvs.write('%f\n' % x_final[i])
    dvs.close()
    postprocess(config, x_final)

def write_sol_1(sol_file,solution):
 
    sol = open(sol_file,'w')
    sol.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nSolAtVertices\n' + str(len(solution)) + '\n1 1\n')
    for i_point in range(len(solution)):
        sol.write(str(solution[i_point]) + "\n")
    sol.write('\nEnd\n')
    sol.close()

def postprocess(config, x_final):

    # local copy
    konfig = copy.deepcopy(config)

    postpro_file = 'postpro.dat'
    postpro = open(postpro_file,'w')

    nx = float(konfig.ACCELERATION_X)
    ny = float(konfig.ACCELERATION_Y)
    nz = float(konfig.ACCELERATION_Z)

    loadFactor = numpy.sqrt(nx*nx+ny*ny+nz*nz)
    gravityVector = -9.81 * numpy.array([-nx/loadFactor,-nz/loadFactor,-ny/loadFactor]) # Change of Frame: to Structure Frame

    # Read bdf

    coord_bdf = []
    elem_bdf = []
    elem_tag_bdf = []
    descriptions = {}
    nDv_bdf = 0
    bdf = open(konfig.STRUCT + '.bdf')
    for line in bdf:
        data = line.split()
        if (line[0]=="$" and len(data) == 3):
            nDv_bdf += 1
            descriptions[data[2].strip().split('/')[0].upper()] = int(data[1])
        if (line[0]=="G" and len(data) == 6):
            vec = [float(data[3]), float(data[4].strip('*'))]
        elif (line[0]=="*" and len(data) == 5):
            vec.append(float(data[2]))
            coord_bdf.append(vec)
        elif (line[0]=="C" and len(data) == 7):
            elem_bdf.append([int(data[3]), int(data[4]), int(data[5]), int(data[6])])
            elem_tag_bdf.append(int(data[2]))
    bdf.close()

    tag_members = []
    for desc in descriptions.keys():
        if desc.split(":")[0] in ['RIB','SPAR','STRING','SKIN','FRAME','LONG']:
            tag_members.append(descriptions[desc])

    nPoint_bdf = len(coord_bdf)
    nElem_bdf = len(elem_bdf)

    # Read load

    nDim = 3

    load = open(konfig.LOAD_FILENAME)

    data = load.readline().split()
    nPoint_load = int(data[0])

    forces = [[0.0 for iDim in range(nDim)] for iPoint_load in range(nPoint_load)]
    for iPoint_load in range(nPoint_load):
        data = load.readline().split()
        forces[iPoint_load][0] = float(data[3])
        forces[iPoint_load][1] = float(data[4])
        forces[iPoint_load][2] = float(data[5])

    load.close()

    # Area

    nNode = 4

    normal_elem_bdf = [[0.0 for iDim in range(nDim)] for iElem_bdf in range(nElem_bdf)]
    area_elem_bdf = [0.0 for iElem_bdf in range(nElem_bdf)]
    vec_a = [0.0 for iDim in range(nDim)] 
    vec_b = [0.0 for iDim in range(nDim)]

    center_elem_bdf = [[0.0 for iDim in range(nDim)] for iElem_bdf in range(nElem_bdf)]

    for iElem_bdf in range(nElem_bdf):
        iPoint_0 = elem_bdf[iElem_bdf][0]-1
        iPoint_1 = elem_bdf[iElem_bdf][1]-1
        iPoint_2 = elem_bdf[iElem_bdf][2]-1
        iPoint_3 = elem_bdf[iElem_bdf][3]-1
        for iDim in range(nDim):
            vec_a[iDim] = coord_bdf[iPoint_0][iDim]-coord_bdf[iPoint_1][iDim]
            vec_b[iDim] = coord_bdf[iPoint_2][iDim]-coord_bdf[iPoint_1][iDim]
            center_elem_bdf[iElem_bdf][iDim] = 0.25*(coord_bdf[iPoint_0][iDim]+coord_bdf[iPoint_1][iDim]+coord_bdf[iPoint_2][iDim]+coord_bdf[iPoint_3][iDim])
        normal_elem_bdf[iElem_bdf][0] += 0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])
        normal_elem_bdf[iElem_bdf][1] += -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0])
        normal_elem_bdf[iElem_bdf][2] += 0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])
        for iDim in range(nDim):
            vec_a[iDim] = coord_bdf[iPoint_2][iDim]-coord_bdf[iPoint_3][iDim]
            vec_b[iDim] = coord_bdf[iPoint_0][iDim]-coord_bdf[iPoint_3][iDim]
        normal_elem_bdf[iElem_bdf][0] += 0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])
        normal_elem_bdf[iElem_bdf][1] += -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0])
        normal_elem_bdf[iElem_bdf][2] += 0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])
        area_elem_bdf[iElem_bdf] += numpy.sqrt(normal_elem_bdf[iElem_bdf][0]*normal_elem_bdf[iElem_bdf][0] + normal_elem_bdf[iElem_bdf][1]*normal_elem_bdf[iElem_bdf][1] + normal_elem_bdf[iElem_bdf][2]*normal_elem_bdf[iElem_bdf][2])

    # Mass and Center Of Mass

    material_rho = float(konfig.MATERIAL_DENSITY)

    structure_mass = 0.0
    com = [0.0 for iDim in range(nDim)]

    for iElem_bdf in range(nElem_bdf):
        elem_thickess = x_final[elem_tag_bdf[iElem_bdf]-1]
        elem_mass = area_elem_bdf[iElem_bdf]*elem_thickess*material_rho
        structure_mass += elem_mass
        for iDim in range(nDim):
            com[iDim] += elem_mass*center_elem_bdf[iElem_bdf][iDim]

    additional_mass = 0.0
    # for iAdd_mass in range(len(additional_masses)):
    #   i_mass = additional_masses[iAdd_mass][3]
    #   additional_mass += i_mass
    #   for iDim in range(nDim):
    #     com[iDim] += i_mass*additional_masses[iAdd_mass][iDim]

    for iDim in range(nDim):
        com[iDim] /= (structure_mass+additional_mass)


    postpro.write('Structural Mass: %f\n' % structure_mass)
    postpro.write('Center of Mass: %f,%f,%f\n' % (com[0],com[1],com[2]))

    # Forces

    external_forces = [0.0 for iDim in range(nDim)] 
    for iPoint_load in range(nPoint_load):
        for iDim in range(nDim):
            external_forces[iDim] += forces[iPoint_load][iDim]

    # Moments (Inertial Forces included in external_forces (additional_mass) will cancel out with the structural_mass Inertial Forces)

    pitch_moment = 0.0
    dist = [0.0 for iDim in range(nDim)]

    for iPoint_load in range(nPoint_load):
        for iDim in range(nDim):
            dist[iDim] = coord_bdf[iPoint_load][iDim]-com[iDim]
        pitch_moment += (forces[iPoint_load][1]*dist[0]-forces[iPoint_load][0]*dist[1])

    pitch_moment_elem = 0.0
    for iElem_bdf in range(nElem_bdf):                             # Contribution of this - is nul if cog computed without additional masses (so external forces cancel out by themselves)
        elem_thickess = x_final[elem_tag_bdf[iElem_bdf]-1]         #                      - will cancel out Inertial Forces included in external_forces if computed with additional masses
        elem_mass = area_elem_bdf[iElem_bdf]*elem_thickess*material_rho
        local_force = elem_mass*gravityVector*loadFactor
        for iDim in range(nDim):
            dist[iDim] = center_elem_bdf[iElem_bdf][iDim]-com[iDim]
        pitch_moment_elem += (local_force[1]*dist[0]-local_force[0]*dist[1])

    # Conclusion: its fine to compute the pitch with the cog of the structre only !!!!
    # The inertial forces included in the External loads (which should have normally have no effect on the pitch) will correct for the offset wrt to the actual cog

    # About accounting for point masses as External Load
    # 3 options to get the pitch:
    #   - Compute actual cog of spaceship and compute pitch due to External forces only
    #   - Compute cog of structure only and compute pitch due to External forces + Inertial forces of left asides masses (tanks, etc ...) : pitch due to Inertial forces will correct the wrong placed cog
    #   - Compute actual cog of spaceship and compute pitch due to External forces + Inertial forces of left asides masses (tanks, etc ...) : pitch due to Inertial forces will cancel out

    # Options 2 and 3 give the same pitch for 2 different cog but same forces distribution ...

    postpro.write('Pitch Moment due to Loads: %f\n' % pitch_moment)
    postpro.write('Pitch Moment due to Gravity (Expected Zero): %f\n' % pitch_moment_elem)

    inertial_forces = structure_mass*gravityVector*loadFactor

    forces_in_non_inertial_frame = external_forces+inertial_forces
    postpro.write('External + Inertial Forces (Expected Zero): %f,%f,%f\n' % (forces_in_non_inertial_frame[0],forces_in_non_inertial_frame[1],forces_in_non_inertial_frame[2]))

    postpro.close()

#: def postprocess()

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


