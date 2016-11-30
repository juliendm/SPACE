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
    spaceutil.surf2pressuresol(konfig)
    SPACE_INT(konfig)
    info = spaceload(konfig)

    # Material properties

    material_rho = float(konfig.MATERIAL_DENSITY)
    material_E = float(konfig.MATERIAL_YOUNG_MODULUS)
    material_ys = 2e7 #float(konfig.MATERIAL_YIELD_STRENGTH)
    material_nu = float(konfig.MATERIAL_POISSON_RATIO)
    kcorr = 5.0/6.0

    t = 1.0
    tMin = 0.0016 # 0.0016
    tMax = 3.0 # 0.020

    nx = float(konfig.ACCELERATION_X)
    ny = float(konfig.ACCELERATION_Y)
    nz = float(konfig.ACCELERATION_Z)

    loadFactor = numpy.sqrt(nx*nx+ny*ny+nz*nz)
    gravityVector = -9.81 * numpy.array([-nx,-nz,-ny])/numpy.sqrt(nx*nx+ny*ny+nz*nz) # Change of Frame: to Structure Frame



    KSWeight = 80.0
    SPs = [StructProblem('lc0', loadFactor=loadFactor, loadFile=konfig.LOAD_FILENAME, evalFuncs=['mass','ks0','mf0'])]
    #SPs = [StructProblem('lc0', loadFactor=loadFactor, loadFile=konfig.LOAD_FILENAME, evalFuncs=['mass','ks0','ks1','ks2'])]
    numLoadCases = len(SPs)

    structOptions = {'transferSize':0.5, 'transferGaussOrder':3}

    FEASolver = pytacs.pyTACS(konfig.STRUCT + '.bdf', comm=comm, options=structOptions)

    # SKIN_FUSE_U = ['FUSE:TOP','FUSE:LFT','CTAIL:LOW','CTAIL_T','FUSE_R']
    # SKIN_FUSE_L = ['FUSE:BOT','FLAP:UPP','FLAP:LOW','FLAP_T','FUSE_F']
    # SKIN_WING_U = ['LWING:UPP','LWING_T::0']
    # SKIN_WING_L = ['LWING:LOW','LWING_T::1']
    # SKINS = SKIN_FUSE_U + SKIN_FUSE_L + SKIN_WING_U + SKIN_WING_L + ['MSKINC:a','MSKINC:b']

    # JUNCTIONS = ['FLAP_FUSE','LWING_FUSE','CTAIL_FUSE']

    # FRAMES = ['MFRAME:00','MFRAME:01','MFRAME:02','MFRAME:03','MFRAME:04','MFRAME:05','MFRAME:06','MFRAME:07','MFRAME:08','MFRAME:09',
    #     'MFRAME:10','MFRAME:11','MFRAME:12']
    # LONGERONS = ['MLONG:02:2','MLONG:00:3','MLONG:01:3','MLONG:02:3','MLONG:00:4','MLONG:01:4']
    # RIBS = ['MRIBF:00','MRIBF:01','MRIBF:02','MRIBF:03','MRIBF:04','MRIBF:05','MRIBF:06','MRIBF:07',
    #     'MRIBV:00','MRIBV:01','MRIBV:02','MRIBV:03','MRIBV:04','MRIBV:05','MRIBV:06','MRIBV:07','MRIBV:08','MRIBV:09',
    #     'MRIBW:00','MRIBW:01','MRIBW:02','MRIBW:03','MRIBW:04','MRIBW:05']
    # SPARS = ['MSPARF:02','MSPARF:03',
    #     'MSPARV:00','MSPARV:01',
    #     'MSPARC:00','MSPARC:04','MSPARC:05',
    #     'MSPARW:00','MSPARW:02','MSPARW:08','MSPARW:09']
    # STRINGERS = ['MSTRINGC:01','MSTRINGC:02','MSTRINGC:03',
    #     'MSTRINGW:01','MSTRINGW:03','MSTRINGW:04','MSTRINGW:05','MSTRINGW:06','MSTRINGW:07']
    # MEMBERS = FRAMES + LONGERONS + RIBS + SPARS + STRINGERS

    # corresp = [-1 for index in range(len(FEASolver.selectCompIDs(include=SKINS+JUNCTIONS+MEMBERS)[0]))]
    # dv = -1;
    # SKIN_IDS = FEASolver.selectCompIDs(include=SKINS)[0]
    # for i in range(len(SKIN_IDS)):
    #     dv_name = "SKIN_" + str(i)
    #     FEASolver.addDVGroup(dv_name, include = SKIN_IDS[i])
    #     dv = dv+1;
    #     corresp[SKIN_IDS[i]] = dv;
    # JUNCTION_IDS = FEASolver.selectCompIDs(include=JUNCTIONS)[0]
    # for i in range(len(JUNCTION_IDS)):
    #     dv_name = "JUNCTION_" + str(i)
    #     FEASolver.addDVGroup(dv_name, include = JUNCTION_IDS[i])
    #     dv = dv+1;
    #     corresp[JUNCTION_IDS[i]] = dv;
    # for i in range(len(MEMBERS)):
    #     dv_name = MEMBERS[i]
    #     FEASolver.addDVGroup(dv_name, include = MEMBERS[i])
    #     dv = dv+1;
    #     MEMBERS_IDS_I = FEASolver.selectCompIDs(include=MEMBERS[i])[0]
    #     for k in range(len(MEMBERS_IDS_I)):
    #         corresp[MEMBERS_IDS_I[k]] = dv;

    #print(corresp)
    #print(len(FEASolver.selectCompIDs(include=SKINS+JUNCTIONS+MEMBERS)[0]))

    ncoms = FEASolver.nComp
    for i in range(0,ncoms):
        dv_name = 'stru_'+str(i)
        FEASolver.addDVGroup(dv_name, include = i)


    def conCallBack(dvNum, compDescripts, userDescript, specialDVs, **kargs):
        con = constitutive.isoFSDTStiffness(material_rho, material_E, material_nu, kcorr, material_ys, t, dvNum, tMin, tMax)
        # if userDescript in ['JUNCTIONS','FRAMES','LONGERONS','SPARS']:
        #     con = constitutive.isoFSDTStiffness(material_rho, material_E, material_nu, kcorr, material_ys, 0.05, dvNum, 0.0016, 0.05)
        scale = [100.0]
        return con, scale

    FEASolver.createTACSAssembler(conCallBack)

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



    def obj(x):
       '''Evaluate the objective and constraints'''
       funcs = {}
       FEASolver.setDesignVars(x)
       for i in range(numLoadCases):
           FEASolver(SPs[i])
           FEASolver.evalFunctions(SPs[i], funcs)
       if comm.rank == 0:
           print funcs
       return funcs, False

    def sens(x, funcs):
       '''Evaluate the objective and constraint sensitivities'''
       funcsSens = {}
       for i in range(numLoadCases):
           FEASolver.evalFunctionsSens(SPs[i], funcsSens)
       return funcsSens, False

    # Set up the optimization problem
    optProb = Optimization('Mass min', obj)
    optProb.addObj('lc0_mass')
    FEASolver.addVariablesPyOpt(optProb)

    for i in range(numLoadCases):
       for j in xrange(1):
           optProb.addCon('%s_ks%d'% (SPs[i].name, j), lower=1.0, upper=1.0)
           #optProb.addCon('%s_mf%d'% (SPs[i].name, j), lower=1.0, upper=1.0)
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
    sol = opt(optProb, sens=sens) #NULL result without error in PyObject_Call





    # # Optimize

    # def obj(x):
    #     '''Evaluate the objective and constraints'''
    #     funcs = {}
    #     FEASolver.setDesignVars(x)
    #     for i in range(numLoadCases):
    #         FEASolver(SPs[i])
    #         FEASolver.evalFunctions(SPs[i], funcs)
    #     if comm.rank == 0:
    #         print funcs
    #     f = funcs['lc0_mass']
    #     g = []
    #     g.append(funcs['lc0_ks0'])
    #     #g.append(funcs['lc0_ks1'])
    #     #g.append(funcs['lc0_ks2'])
    #     fail = 0
    #     return f, g, fail

    # def sens(x, f, g):
    #     '''Evaluate the objective and constraint sensitivities'''
    #     funcsSens = {}
        
    #     x_cur = numpy.zeros(len(x))
    #     FEASolver.structure.getDesignVars(x_cur)
    #     if not numpy.array_equal(x_cur,x):
    #         print 'need recompute'
    #         FEASolver.setDesignVars(x)
    #         for i in range(numLoadCases):
    #             FEASolver(SPs[i])

    #     for i in range(numLoadCases):
    #         FEASolver.evalFunctionsSens(SPs[i], funcsSens)

    #     df1 = funcsSens['lc0_mass'][FEASolver.varSet]

    #     dg1 = []
    #     dg1.append(funcsSens['lc0_ks0'][FEASolver.varSet])
    #     #dg1.append(funcsSens['lc0_ks1'][FEASolver.varSet])
    #     #dg1.append(funcsSens['lc0_ks2'][FEASolver.varSet])

    #     df = [0.0]*len(df1)
    #     for i in range(0,len(df1)):
    #         df[i] = df1[i] # NEEDED ?????????????????

    #     dg = numpy.zeros([len(dg1),len(dg1[0])])
    #     for i in range(0,len(dg1[0])):
    #         dg[0][i] = dg1[0][i]
    #         #dg[1][i] = dg1[1][i]
    #         #dg[2][i] = dg1[2][i]

    #     fail = 0

    #     return df, dg, fail


    # optimize = True
    # if optimize:

    #     # Set up the optimization problem
    #     optProb = Optimization('Mass min', obj)
    #     optProb.addObj('lc0_mass')
    #     FEASolver.addVariablesPyOpt(optProb)

    #     for i in range(numLoadCases):
    #         for j in xrange(1):   ###################################################
    #             optProb.addCon('%s_ks%d'% (SPs[i].name, j), upper=1.0)

    #     if comm.rank == 0:
    #         print optProb
    #     #optProb.printSparsity()


    #     # SNOPT
    #     snopt = SNOPT()
    #     #snopt.setOption('Major iterations limit',3)
    #     [obj_fun, x_dvs, inform] = snopt(optProb, sens_type=sens)

    #     # # IPOPT
    #     # ipopt = IPOPT()
    #     # [obj_fun, x_dvs, inform] = snopt(optProb, sens_type=sens)

    #     # # SLSQP
    #     # slsqp = SLSQP()
    #     # [obj_fun, x_dvs, inform] = slsqp(optProb, sens_type=sens)





    #     #print optProb.solution(0)
    #     optProb.write2file(outfile='pySNOPT.txt', disp_sols=False, solutions=[0])

    #     # Setting the final solution

    #     FEASolver.setDesignVars(x_dvs) # NEEDED ?????????????????

    #     # dvs_file = 'dvs.dat'
    #     # dvs = open(dvs_file,'w')
    #     # for i in range(len(corresp)):
    #     #     dvs.write('%f\n' % x_dvs[corresp[i]])
    #     # dvs.close()





    # Getting the final solution

    for i in range(numLoadCases):
        FEASolver(SPs[i]) # NEEDED ?????????????????



    # funcs = {}
    # FEASolver.evalFunctions(SPs[0], funcs)

    FEASolver.writeBDFForces(SPs[0], "visualize_forces.bdf")
    FEASolver.writeSolution()
    # #FEASolver.writeBDF("crm_wing_design_tacs_final.bdf")


    #postprocess(config)

    # info out
    info = spaceio.State(info)
    #info.FILES.MASS = konfig.MASS
    return info

#: def structure()

def postprocess(config):

    # local copy
    konfig = copy.deepcopy(config)




    nx = float(konfig.ACCELERATION_X)
    ny = float(konfig.ACCELERATION_Y)
    nz = float(konfig.ACCELERATION_Z)

    loadFactor = np.sqrt(nx*nx+ny*ny+nz*nz)
    gravityVector = -9.81 * np.array([-nx/loadFactor,-nz/loadFactor,-ny/loadFactor]) # Change of Frame: to Structure Frame

    # Read bdf

    coord_bdf = []
    elem_bdf = []
    elem_tag_bdf = []
    descriptions = {}
    nDv_bdf = 0
    bdf = open(konfig.STRUCT_SURFACE + '.bdf')
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

    # Read dvs

    dvs_file = 'dvs.dat'
    dvs = open(dvs_file)
    x_dvs = np.loadtxt(dvs);
    dvs.close()
    print len(x_dvs)

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
        area_elem_bdf[iElem_bdf] += np.sqrt(normal_elem_bdf[iElem_bdf][0]*normal_elem_bdf[iElem_bdf][0] + normal_elem_bdf[iElem_bdf][1]*normal_elem_bdf[iElem_bdf][1] + normal_elem_bdf[iElem_bdf][2]*normal_elem_bdf[iElem_bdf][2])

    # Mass and Center Of Mass

    material_rho = float(konfig.MATERIAL_DENSITY)

    structure_mass = 0.0
    com = [0.0 for iDim in range(nDim)]

    for iElem_bdf in range(nElem_bdf):
        elem_thickess = x_dvs[elem_tag_bdf[iElem_bdf]-1]
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

    print structure_mass
    print com

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
    for iElem_bdf in range(nElem_bdf):                           # Contribution of this is nul if com computed without additional masses (so external forces cancel out by themselves)
        elem_thickess = x_dvs[elem_tag_bdf[iElem_bdf]-1]         #                      will cancel out Inertial Forces included in external_forces if computed with additional masses
        elem_mass = area_elem_bdf[iElem_bdf]*elem_thickess*material_rho
        local_force = elem_mass*gravityVector*loadFactor
        for iDim in range(nDim):
            dist[iDim] = center_elem_bdf[iElem_bdf][iDim]-com[iDim]
        pitch_moment_elem += (local_force[1]*dist[0]-local_force[0]*dist[1])

    # Conclusion: its fine to compute the pitch with the cog of the structre only !!!!

    print pitch_moment, pitch_moment_elem

    inertial_forces = structure_mass*gravityVector*loadFactor

    print external_forces+inertial_forces

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


