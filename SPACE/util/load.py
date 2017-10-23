#!/usr/bin/env python2.7

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, time, sys, shutil, copy
import numpy as np

# ----------------------------------------------------------------------
#  Load Simulation
# ----------------------------------------------------------------------

class Load(object):

    def __init__(self, config, loadFactor, gravity_vector, pdyn_inf, half_thrust, thrust_angle = 0.0, fuel_percentage = 1.0, 
        safetyFactor_thrust = 1.0, safetyFactor_inertial = 1.0, safetyFactor_non_inertial = 1.0):

        self._nDim = 3
        self._nNode = 4

        self._nFrame = 13
        self._nLongeron = 4

        self._thrust_frames   = [10]
        self._avionics_frames = [0,1]
        self._lox_frames      = [7,10]
        self._kero_frames     = [1,3]
        self._engine_frames   = [10,11,16]
        self._payload_frames  = [3,4,5,6,7]

        self._config = copy.deepcopy(config)

        self._loadFactor = loadFactor
        self._gravity_vector = gravity_vector

        self._half_thrust_newtons = half_thrust
        self._thrust_angle = thrust_angle # 0 degrees for axial thrust

        self._pdyn_inf = pdyn_inf
        self._fuel_percentage = fuel_percentage

        self._safetyFactor_thrust = safetyFactor_thrust
        self._safetyFactor_inertial = safetyFactor_inertial
        self._safetyFactor_non_inertial = safetyFactor_non_inertial


        self._material_rho = float(self._config.MATERIAL_DENSITY)

        self._half_mass_payload = float(self._config.HALF_PAYLOAD_MASS)
        self._half_mass_fuel_kero = float(self._config.HALF_KERO_MASS)
        self._half_mass_fuel_lox = float(self._config.HALF_LOX_MASS)

        # Read Meshes

        self._coord_bdf, self._elem_bdf, self._elem_tag_bdf, self._descriptions = read_bdf(self._config , self._nDim, self._nNode) # OML + Structure
        self._coord, self._elem, self._bdf_corresp                              = read_mesh(self._config, self._nDim, self._nNode) # OML only

        self._nPoint_bdf = len(self._coord_bdf)
        self._nPoint = len(self._coord)

        self._nElem_bdf = len(self._elem_bdf)
        self._nElem = len(self._elem)

        # Read Solution

        pressureCoeff, frictionCoeff = read_sol(self._config, self._nDim, self._nNode)  # OML only
        self._pressureCoeff_bdf = [0.0 for iPoint_bdf in range(self._nPoint_bdf)]
        self._frictionCoeff_bdf = [[0.0 for iDim in range(self._nDim)] for iPoint_bdf in range(self._nPoint_bdf)]
        for iPoint in range(self._nPoint):
            self._pressureCoeff_bdf[self._bdf_corresp[iPoint]] = pressureCoeff[iPoint]
            for iDim in range(self._nDim):
                self._frictionCoeff_bdf[self._bdf_corresp[iPoint]][iDim] = frictionCoeff[iPoint][iDim]

        # Compute Area and Normal

        area, normal = areas_normals_voronoi(self._nDim, self._nNode, self._coord, self._elem)
        self._area_voronoi = [0.0 for iPoint_bdf in range(self._nPoint_bdf)]
        self._normal_voronoi = [[0.0 for iDim in range(self._nDim)] for iPoint_bdf in range(self._nPoint_bdf)]
        for iPoint in range(self._nPoint):
            self._area_voronoi[self._bdf_corresp[iPoint]] = area[iPoint]
            for iDim in range(self._nDim):
                self._normal_voronoi[self._bdf_corresp[iPoint]][iDim] = normal[iPoint][iDim]

        self._area_elem, self._normal_elem, self._center_elem = areas_normals_elem(self._nDim, self._nNode, self._coord_bdf, self._elem_bdf)

        # Vehicle Dimensions

        coord_array = np.array(self._coord)
        self._body_length = 18.0 # abs(max(coord_array[:,0])-min(coord_array[:,0])) # FIXED TO NOT INCLUDE VARIATION DUE TO BODYFLAP !!!!!!!!!!
        self._wing_span = 2.0*max(abs(coord_array[:,1]))

        # Apply

        self.__compute_apply_noninertial()
        self.__compute_apply_inertial()

        # Non Inertial Load

        self.__compte_load_noninertial()

        # Structural Mass

        self._structural_mass_bdf = [0.0 for iPoint_bdf in range(len(self._coord_bdf))]

        # sum_x = 0.0
        # sum_y = 0.0
        # sum_z = 0.0
        # sum_area = 0.0
        # for iPoint_bdf in range(nPoint_bdf):
        #     sum_x += normal_bdf[iPoint_bdf][0]
        #     sum_y += normal_bdf[iPoint_bdf][1]
        #     sum_z += normal_bdf[iPoint_bdf][2]
        #     sum_area += area_bdf[iPoint_bdf]
        # print sum_x
        # print sum_y
        # print sum_z
        # print sum_area

        # # Control Forces

        # aero_forces = [0.0 for iDim in range(nDim)]
        # for iPoint_bdf in range(nPoint_bdf):

        #     # Pressure
            
        #     pressure = float(konfig.P_DYN_INF)*pressureCoeff_bdf[iPoint_bdf] # + float(konfig.P_INF) # NOT ADDING p_inf CAUSE SPACEPLANE IS NOT PRESSURIZED
        #     for iDim in range(nDim):
        #         aero_forces[iDim] += normal_bdf[iPoint_bdf][iDim]*pressure

        #     # Friction
            
        #     for iDim in range(nDim):
        #         shear_stress = float(konfig.P_DYN_INF)*frictionCoeff_bdf[iPoint_bdf][iDim]
        #         aero_forces[iDim] += area_bdf[iPoint_bdf]*shear_stress

        # for iDim in range(nDim):
        #     aero_forces[iDim] /= float(konfig.P_DYN_INF) * float(konfig.REF_AREA)

        # print "AERO FORCES: ", aero_forces

        # aoa = float(konfig.AoA)*np.pi/180.0

        # lift_coeff = ( -aero_forces[0]*np.sin(aoa) + aero_forces[1]*np.cos(aoa) )
        # drag_coeff = ( aero_forces[0]*np.cos(aoa) + aero_forces[1]*np.sin(aoa) )

        # print lift_coeff, drag_coeff


    def update(self, half_mass_kg_current_step):

        self._load_bdf = [[0.0 for iDim in range(self._nDim)] for iPoint_bdf in range(self._nPoint_bdf)]

        # NON INERTIAL

        for iPoint_bdf in range(self._nPoint_bdf):
            for iDim in range(self._nDim):
                self._load_bdf[iPoint_bdf][iDim] = self._load_noninertial_bdf[iPoint_bdf][iDim]

        # INERTIAL

        self.__update_additional_mass(half_mass_kg_current_step)

        for iPoint_bdf in range(self._nPoint_bdf):
            for iDim in range(self._nDim):

                # # DO NOT ADD THIS WHEN USING NASTRAN, IT WILL TAKE CARE OF THIS
                # # Structural Mass
                # self._load_bdf[iPoint_bdf][iDim] += self._structural_mass_bdf[iPoint_bdf]*self._gravity_vector[iDim]*self._loadFactor*self._safetyFactor_inertial

                # Additional Mass
                self._load_bdf[iPoint_bdf][iDim] += self._additional_mass_bdf[iPoint_bdf]*self._gravity_vector[iDim]*self._loadFactor*self._safetyFactor_inertial

        # Write load
        write_load(self._config.LOAD_FILENAME,self._load_bdf,self._coord_bdf,self._elem_bdf,self._nNode)

        # Write Check
        write_check(self._load_bdf,self._coord_bdf,self._elem_bdf,self._normal_voronoi,self._nNode)

    #: def load()

    def postprocess(self, x_dvs, corresp):

        # # DO NOT ADD THIS WHEN USING NASTRAN, IT WILL TAKE CARE OF THIS
        # # Structural Mass
        # self.__update_structural_mass(x_dvs, corresp)

        # Mass and Center Of Mass

        self._center_of_mass = [0.0 for iDim in range(self._nDim)]

        self._half_structure_mass = 0.0
        for iElem_bdf in range(self._nElem_bdf):
            elem_thickess = x_dvs[corresp[self._elem_tag_bdf[iElem_bdf]-1]-1]
            elem_mass = self._area_elem[iElem_bdf]*elem_thickess*self._material_rho
            self._half_structure_mass += elem_mass
            for iDim in range(self._nDim):
                self._center_of_mass[iDim] += elem_mass*self._center_elem[iElem_bdf][iDim]

        self._half_additional_mass = 0.0
        for iPoint_bdf in range(self._nPoint_bdf):
            point_mass = self._additional_mass_bdf[iPoint_bdf]
            self._half_additional_mass += point_mass
            for iDim in range(self._nDim):
                self._center_of_mass[iDim] += point_mass*self._coord_bdf[iPoint_bdf][iDim]

        for iDim in range(self._nDim):
            self._center_of_mass[iDim] /= (self._half_structure_mass+self._half_additional_mass)

        # Forces

        inertial_forces = (self._half_structure_mass+self._half_additional_mass)*self._gravity_vector*self._loadFactor*self._safetyFactor_inertial

        non_inertial_forces = [0.0 for iDim in range(self._nDim)] 
        for iPoint_bdf in range(self._nPoint_bdf):
            for iDim in range(self._nDim):
                non_inertial_forces[iDim] += self._load_noninertial_bdf[iPoint_bdf][iDim]

        # Moments (Inertial Forces included in non_inertial_forces (additional_mass) would cancel out with the structural_mass Inertial Forces)

        self._pitch_moment = 0.0
        dist = [0.0 for iDim in range(self._nDim)]
        for iPoint_bdf in range(self._nPoint_bdf):
            for iDim in range(self._nDim):
                dist[iDim] = self._coord_bdf[iPoint_bdf][iDim]-self._center_of_mass[iDim]
            self._pitch_moment += (self._load_noninertial_bdf[iPoint_bdf][2]*dist[0]-self._load_noninertial_bdf[iPoint_bdf][0]*dist[2])

        # Center Of Pressure

        self._center_of_pressure = [0.0 for iDim in range(self._nDim)]
        for iPoint_bdf in range(self._nPoint_bdf):
            self._center_of_pressure[0] += self._load_noninertial_bdf[iPoint_bdf][2]*self._coord_bdf[iPoint_bdf][0]
            self._center_of_pressure[2] += self._load_noninertial_bdf[iPoint_bdf][0]*self._coord_bdf[iPoint_bdf][2]
        self._center_of_pressure[0] /= non_inertial_forces[2]
        self._center_of_pressure[2] /= non_inertial_forces[0]

        # pitch_moment_elem = 0.0
        # for iElem_bdf in range(self._nElem_bdf):                             # Contribution of this - is nul if cog computed without additional masses (so external forces cancel out by themselves)
        #     elem_thickess = x_final[self._elem_tag_bdf[iElem_bdf]-1]         #                      - will cancel out Inertial Forces included in non_inertial_forces if computed with additional masses
        #     elem_mass = self._area_elem[iElem_bdf]*elem_thickess*material_rho
        #     local_force = elem_mass*self._gravityVector*self._loadFactor*self._safetyFactor_inertial
        #     for iDim in range(self._nDim):
        #         dist[iDim] = self._center_elem[iElem_bdf][iDim]-com[iDim]
        #     pitch_moment_elem += (local_force[2]*dist[0]-local_force[0]*dist[2])

        # Conclusion: its fine to compute the pitch with the cog of the structre only !!!!
        # The inertial forces included in the External loads (which should have normally have no effect on the pitch) will correct for the offset wrt to the actual cog

        # About accounting for point masses as External Load
        # 3 options to get the pitch:
        #   - Compute actual cog of spaceship and compute pitch due to External forces only
        #   - Compute cog of structure only and compute pitch due to External forces + Inertial forces of left asides masses (tanks, etc ...) : pitch due to Inertial forces will correct the wrong placed cog
        #   - Compute actual cog of spaceship and compute pitch due to External forces + Inertial forces of left asides masses (tanks, etc ...) : pitch due to Inertial forces will cancel out

        # Options 2 and 3 give the same pitch for 2 different cog but same forces distribution ...

        postpro_file = 'postpro.dat'
        postpro = open(postpro_file,'w')

        postpro.write('Half Structural Mass: %f\n' % self._half_structure_mass)
        postpro.write('Half Additional Mass: %f\n' % self._half_additional_mass)
        postpro.write('Half Mass           : %f\n' % (self._half_structure_mass+self._half_additional_mass))
        postpro.write('\n')

        postpro.write('Center of Mass    : %f,%f\n' % (self._center_of_mass[0],self._center_of_mass[2]))
        postpro.write('Center of Pressure: %f,%f\n' % (self._center_of_pressure[0],self._center_of_pressure[2]))
        postpro.write('Distance CoG-CoP x: %f\n'    % (np.sqrt((self._center_of_mass[0]-self._center_of_pressure[0])**2.0)))
        postpro.write('Distance CoG-CoP z: %f\n'    % (np.sqrt((self._center_of_mass[2]-self._center_of_pressure[2])**2.0)))
        postpro.write('\n')

        forces_in_non_inertial_frame = non_inertial_forces+inertial_forces
        postpro.write('Sum Forces  (Expected Zero): %f,%f\n' % (forces_in_non_inertial_frame[0],forces_in_non_inertial_frame[2]))
        postpro.write('Sum Moments (Expected Zero): %f\n'    % self._pitch_moment)
        postpro.write('\n')

        postpro.close()

    #: def postprocess()

    def __update_structural_mass(self, x_dvs, corresp):

        point_thickess = [0.0 for iPoint_bdf in range(len(self._coord_bdf))]
        point_thickess_count = [0 for iPoint_bdf in range(len(self._coord_bdf))]

        for iElem_bdf in range(self._nElem_bdf):
            for iNode in range(self._nNode):
                iPoint_bdf = self._elem_bdf[iElem_bdf][iNode]-1
                point_thickess[iPoint_bdf] += x_dvs[corresp[self._elem_tag_bdf[iElem_bdf]-1]-1]
                point_thickess_count[iPoint_bdf] += 1

        for iPoint_bdf in range(len(self._coord_bdf)):
            point_thickess[iPoint_bdf] /= point_thickess_count[iPoint_bdf]

        for iPoint_bdf in range(len(self._coord_bdf)):
            point_mass = self._area_voronoi[iPoint_bdf]*point_thickess[iPoint_bdf]*self._material_rho
            self._structural_mass_bdf[iPoint_bdf] = point_mass

    def __update_additional_mass(self, half_mass_kg_current_step):

        self.__update_half_additional_mass_kg_next_step(half_mass_kg_current_step)

        self._additional_mass_bdf = [0.0 for iPoint_bdf in range(len(self._coord_bdf))]

        surface_equip = get_apply_surface(self._apply_equip, self._area_voronoi)

        for iPoint_bdf in range(len(self._coord_bdf)):

            # Landing Gear

            if iPoint_bdf in self._apply_gear:

                self._additional_mass_bdf[iPoint_bdf] += self._half_mass_gear/len(self._apply_gear)

            # Hydraulic

            if iPoint_bdf in self._apply_hydraulic:

                self._additional_mass_bdf[iPoint_bdf] += self._half_mass_hydraulic/len(self._apply_hydraulic)

            # Avionics

            if iPoint_bdf in self._apply_avionics:

                self._additional_mass_bdf[iPoint_bdf] += self._half_mass_avionics/len(self._apply_avionics)

            # Electrical System

            if iPoint_bdf in self._apply_elec:

                self._additional_mass_bdf[iPoint_bdf] += self._half_mass_elec/len(self._apply_elec)

            # Equipment

            if iPoint_bdf in self._apply_equip:

                self._additional_mass_bdf[iPoint_bdf] += self._area_voronoi[iPoint_bdf]*self._half_mass_equip/surface_equip

            # Tank LOX

            if iPoint_bdf in self._apply_lox:

                self._additional_mass_bdf[iPoint_bdf] += self._half_mass_tank_lox/len(self._apply_lox)

            # Fuel LOX - Note: Assume fluid-air interface always perpendicular to longitudinal axe

            if iPoint_bdf in self._apply_lox_front:
                self._additional_mass_bdf[iPoint_bdf] += (self._fuel_percentage*0.5)*self._fuel_percentage*self._half_mass_fuel_lox/len(self._apply_lox_front)
            if iPoint_bdf in self._apply_lox_rear:
                self._additional_mass_bdf[iPoint_bdf] += (1.0-self._fuel_percentage*0.5)*self._fuel_percentage*self._half_mass_fuel_lox/len(self._apply_lox_rear)

            # Tank KERO

            if iPoint_bdf in self._apply_kero:

                self._additional_mass_bdf[iPoint_bdf] += self._half_mass_tank_kero/len(self._apply_kero)

            # Fuel KERO - Note: Assume fluid-air interface always perpendicular to longitudinal axe

            if iPoint_bdf in self._apply_kero_front:
                self._additional_mass_bdf[iPoint_bdf] += (self._fuel_percentage*0.5)*self._fuel_percentage*self._half_mass_fuel_kero/len(self._apply_kero_front)
            if iPoint_bdf in self._apply_kero_rear:
                self._additional_mass_bdf[iPoint_bdf] += (1.0-self._fuel_percentage*0.5)*self._fuel_percentage*self._half_mass_fuel_kero/len(self._apply_kero_rear)

            # Engine

            if iPoint_bdf in self._apply_engine:

                self._additional_mass_bdf[iPoint_bdf] += self._half_mass_engine/len(self._apply_engine)

            # Thermal Protection System

            if iPoint_bdf in self._apply_tps:

                # https://science.ksc.nasa.gov/shuttle/technology/sts-newsref/sts-tps.html
                # 9 pounds per cubic foot = 144.166 kg/m^3
                # thickness from 1 inch to 5 inches -> 3 inches = 0.0762 m
                # generally, the HRSI tiles are thicker at the forward areas of the orbiter and thinner toward the aft end

                density_tps = 144.166 # kg/m^3
                thickness_tps = 0.0762 # m
                self._additional_mass_bdf[iPoint_bdf] += self._area_voronoi[iPoint_bdf]*thickness_tps*density_tps

            # Payload

            if iPoint_bdf in self._apply_payload:

                self._additional_mass_bdf[iPoint_bdf] += self._half_mass_payload/len(self._apply_payload)



    def __update_half_additional_mass_kg_next_step(self, half_mass_kg_current_step):

        pounds_to_kg = 0.453592
        newtons_to_pounds = 0.224809
        kg_to_pounds = 2.20462
        meters_to_feet = 3.28084

        weight_pounds_current_step = 2.0*half_mass_kg_current_step*kg_to_pounds
        thrust_pounds = 2.0*self._half_thrust_newtons*newtons_to_pounds

        body_length_feet = self._body_length*meters_to_feet
        wing_span_feet = self._wing_span*meters_to_feet

        wing_body_elevon_surface_square_feet = 2.0*get_apply_surface(self._apply_wing_body_elevon, self._area_voronoi)*meters_to_feet*meters_to_feet
        body_flap_surface_square_feet        = 2.0*get_apply_surface(self._apply_body_flap, self._area_voronoi)*meters_to_feet*meters_to_feet
        vertical_tail_surface_square_feet    = 2.0*get_apply_surface(self._apply_vertical_tail, self._area_voronoi)*meters_to_feet*meters_to_feet


        lox_density = 1141 # kg/m3
        kero_density = 810 # kg/m3

        N_engines = 1
        rocket_expansion_ratio = 77.5 #TODO: CHECK

        # Space Shuttle External Tank
        # Length: 153.8 ft (46.9 m)
        # Diameter: 27.6 ft (8.4 m)
        # Empty Weight: 58,500 lb (26,500 kg)
        # Gross Liftoff Weight: 1,680,000 lb (760,000 kg)
        sts_tank_radius = 4.2 # m
        sts_tank_height = 46.9 # m
        sts_tank_area = 2.0*np.pi*sts_tank_radius*sts_tank_height + 2.0*np.pi*sts_tank_radius**2.0 # m2 # Assume: Right Cylinder
        sts_tank_empty_mass = 26500.0 # kg
        tank_mass_per_area = sts_tank_empty_mass/sts_tank_area # kg/m2

        # Landing Gear Weight

        weight_gear = 0.00916*weight_pounds_current_step**1.124
        self._half_mass_gear = weight_gear*0.5*pounds_to_kg

        # Hydraulic Weight

        weight_hydraulic = 2.64 * ( ( (wing_body_elevon_surface_square_feet + body_flap_surface_square_feet + vertical_tail_surface_square_feet)*self._pdyn_inf/1000.0)**0.334 * (body_length_feet + wing_span_feet)**0.5 )
        self._half_mass_hydraulic = weight_hydraulic*0.5*pounds_to_kg

        # Avionics Weight

        weight_avionics = 66.37*weight_pounds_current_step**0.361
        self._half_mass_avionics = weight_avionics*0.5*pounds_to_kg

        # Electrical System Weight

        weight_elec = 1.167*weight_pounds_current_step**0.5*body_length_feet**0.25
        self._half_mass_elec = weight_elec*0.5*pounds_to_kg

        # Equipment Weight

        weight_equip = 1000.0 + 0.01*weight_pounds_current_step ############## Maybe reduce fixed value
        self._half_mass_equip = weight_equip*0.5*pounds_to_kg

        # Tank LOX Weight

        volume_lox = 2.0*self._half_mass_fuel_lox/lox_density # m3
        lox_tank_radius = 1.65 # m
        lox_tank_height = volume_lox/np.pi/lox_tank_radius/lox_tank_radius
        lox_tank_area = 2.0*np.pi*lox_tank_radius*lox_tank_height + 2.0*np.pi*lox_tank_radius**2.0 # m2 # Assume: Right Cylinder

        self._half_mass_tank_lox = 0.5*tank_mass_per_area*lox_tank_area # kg

        # Tank KERO Weight

        volume_kero = 2.0*self._half_mass_fuel_kero/kero_density # m3
        kero_tank_radius = 1.65 # m
        kero_tank_height = volume_kero/np.pi/kero_tank_radius/kero_tank_radius
        kero_tank_area = 2.0*np.pi*kero_tank_radius*kero_tank_height + 2.0*np.pi*kero_tank_radius**2.0 # m2 # Assume: Right Cylinder

        self._half_mass_tank_kero = 0.5*tank_mass_per_area*kero_tank_area # kg

        # Engine Weight

        weight_engine = 0.00766*thrust_pounds + 0.00033*thrust_pounds*rocket_expansion_ratio**0.5 + 130.0*N_engines
        self._half_mass_engine = weight_engine*0.5*pounds_to_kg


    def __compte_load_noninertial(self):

        self._load_noninertial_bdf = [[0.0 for iDim in range(self._nDim)] for iPoint_bdf in range(self._nPoint_bdf)]

        # NON INERTIAL

        for iPoint_bdf in range(self._nPoint_bdf):

            # AERO

            if not iPoint_bdf in self._apply_fuse_r:

                # Pressure
                
                pressure = self._pdyn_inf*self._pressureCoeff_bdf[iPoint_bdf] # + float(konfig.P_INF) # NOT ADDING p_inf CAUSE SPACEPLANE IS NOT PRESSURIZED
                for iDim in range(self._nDim):
                    self._load_noninertial_bdf[iPoint_bdf][iDim] -= self._normal_voronoi[iPoint_bdf][iDim]*pressure * self._safetyFactor_non_inertial # MINUS SIGN CAUSE PRESSURE PUSHES INWARD

                # Friction
                
                for iDim in range(self._nDim):
                    shear_stress = self._pdyn_inf*self._frictionCoeff_bdf[iPoint_bdf][iDim]
                    self._load_noninertial_bdf[iPoint_bdf][iDim] += self._area_voronoi[iPoint_bdf]*shear_stress * self._safetyFactor_non_inertial

            # THRUST

            # Truss Problem

            # joints are assumed to act like hinges, permitting free rotation of the bars around the joint
            # furthermore assumed that the truss structure is only loaded by concentrated forces acting at the joints
            # as a consequence of the assumption of hinges the bar elements can only support an axial force

            # Same signe conventions as for Elevons and Body Flap
            force_x = -self._half_thrust_newtons*self._safetyFactor_thrust*np.cos(self._thrust_angle*np.pi/180.0)
            force_z = self._half_thrust_newtons*self._safetyFactor_thrust*np.sin(self._thrust_angle*np.pi/180.0)

            # Assume hinges are at 45 degrees
            F0 = 0.5*(force_x-force_z)
            F1 = 0.5*(force_x+force_z)

            if iPoint_bdf in self._apply_thrust_0:
                self._load_noninertial_bdf[iPoint_bdf][0] += F0/len(self._apply_thrust_0)
                self._load_noninertial_bdf[iPoint_bdf][2] -= F0/len(self._apply_thrust_0)
            if iPoint_bdf in self._apply_thrust_1:
                self._load_noninertial_bdf[iPoint_bdf][0] += F1/len(self._apply_thrust_1)
                self._load_noninertial_bdf[iPoint_bdf][2] += F1/len(self._apply_thrust_1)


    def __compute_apply_noninertial(self):

        # THRUST

        tag_thrust_frames = []
        for i in range(self._nLongeron-1):
            for j in self._thrust_frames:
                for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                    if desc in self._descriptions.keys():
                        tag_thrust_frames.append(self._descriptions[desc])

        tag_longerons_0 = []
        for i in [0]:
            for j in range(self._nFrame-1):
                for desc in ['MLONG:%02d:3:%02d' % (i,j)]:
                    if desc in self._descriptions.keys():
                        tag_longerons_0.append(self._descriptions[desc])

        tag_longerons_1 = []
        for i in [0]:
            for j in range(self._nFrame-1):
                for desc in ['MLONG:%02d:4:%02d' % (i,j)]:
                    if desc in self._descriptions.keys():
                        tag_longerons_1.append(self._descriptions[desc])

        tag_members = []
        tag_skin = []
        for desc in self._descriptions.keys():
            if desc.split(":")[0] in ['MRIBF','MRIBV','MRIBW','MSPARF','MSPARV','MSPARC','MSPARW','MSTRINGC','MSTRINGW','MSKINC','MFRAME','MLONG']:
                tag_members.append(self._descriptions[desc])
            else:
                tag_skin.append(self._descriptions[desc])


        apply_thrust_frames = tag_to_apply(tag_thrust_frames, self._elem_bdf, self._elem_tag_bdf)
        apply_longerons_0   = tag_to_apply(tag_longerons_0  , self._elem_bdf, self._elem_tag_bdf)
        apply_longerons_1   = tag_to_apply(tag_longerons_1  , self._elem_bdf, self._elem_tag_bdf)
        apply_skin          = tag_to_apply(tag_skin         , self._elem_bdf, self._elem_tag_bdf)

        #self._apply_thrust = intersection(intersection(apply_thrust_frames, apply_longerons), apply_skin)
        self._apply_thrust_0 = intersection(apply_thrust_frames, apply_longerons_0)
        self._apply_thrust_1 = intersection(apply_thrust_frames, apply_longerons_1)

        # EXHAUST

        tag_fuse_r = []
        for desc in self._descriptions.keys():
            if 'FUSE_R' in desc:
                tag_fuse_r.append(self._descriptions[desc])

        self._apply_fuse_r = tag_to_apply(tag_fuse_r, self._elem_bdf, self._elem_tag_bdf)

    def __compute_apply_inertial(self):

        # TAGS

        tag_longerons = []
        for i in [0]:
            for j in range(self._nFrame-1):
                for desc in ['MLONG:%02d:3:%02d' % (i,j), 'MLONG:%02d:4:%02d' % (i,j)]:
                    if desc in self._descriptions.keys():
                        tag_longerons.append(self._descriptions[desc])

        tag_longerons_payload = []
        for i in [2]:
            for j in range(self._nFrame-1):
                for desc in ['MLONG:%02d:3:%02d' % (i,j)]:
                    if desc in self._descriptions.keys():
                        tag_longerons_payload.append(self._descriptions[desc])

        tag_avionics_frames = []
        tag_engine_frames = []
        tag_payload_frames = []

        tag_lox_front_frames = []
        tag_kero_front_frames = []
        tag_lox_rear_frames = []
        tag_kero_rear_frames = []

        for i in range(self._nLongeron-1):

            for j in self._avionics_frames:
                for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                    if desc in self._descriptions.keys():
                        tag_avionics_frames.append(self._descriptions[desc])
            for j in self._engine_frames:
                for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                    if desc in self._descriptions.keys():
                        tag_engine_frames.append(self._descriptions[desc])
            for j in self._payload_frames:
                for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                    if desc in self._descriptions.keys():
                        tag_payload_frames.append(self._descriptions[desc])

            j = self._lox_frames[0]
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in self._descriptions.keys():
                    tag_lox_front_frames.append(self._descriptions[desc])
            j = self._lox_frames[1]
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in self._descriptions.keys():
                    tag_lox_rear_frames.append(self._descriptions[desc])
            j = self._kero_frames[0]
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in self._descriptions.keys():
                    tag_kero_front_frames.append(self._descriptions[desc])
            j = self._kero_frames[1]
            for desc in ['MFRAME:%02d:2:%02d' % (j,i), 'MFRAME:%02d:3:%02d' % (j,i), 'MFRAME:%02d:4:%02d' % (j,i)]:
                if desc in self._descriptions.keys():
                    tag_kero_rear_frames.append(self._descriptions[desc])

        tag_wing_body_elevon = []
        SKIN_FUSE_L = ['FUSE:BOT','FUSE_F']
        SKIN_WING_L = ['LWING:LOW','LWING_T::1']
        SKINS_L = SKIN_FUSE_L + SKIN_WING_L 
        for desc in self._descriptions.keys():
            for val in SKINS_L:
                if val in desc:
                    tag_wing_body_elevon.append(self._descriptions[desc])

        tag_body_flap = []
        for desc in self._descriptions.keys():
            for val in ['FLAP:UPP']:
                if val in desc:
                    tag_body_flap.append(self._descriptions[desc])

        tag_vertical_tail = []
        for desc in self._descriptions.keys():
            for val in ['CTAIL:LOW','CTAIL_T::1']:
                if val in desc:
                    tag_vertical_tail.append(self._descriptions[desc])

        tag_hydraulic = []
        for desc in self._descriptions.keys():
            for val in ['MSPARW:02','MRIBF:00','MSPARV:01']:
                if val in desc:
                    tag_hydraulic.append(self._descriptions[desc])

        tag_gear = []
        for desc in self._descriptions.keys():
            for val in ['MSTRINGW:03:A:L:02','MSTRINGW:03:A:L:03',
                        'MSTRINGW:03:B:L:02','MSTRINGW:03:B:L:03',
                        'MSTRINGW:04:A:L:02','MSTRINGW:04:A:L:03',
                        'MSTRINGW:04:B:L:02','MSTRINGW:04:B:L:03',
                        'MSTRINGW:05:A:L:02','MSTRINGW:05:A:L:03',
                        'MSTRINGW:05:B:L:02','MSTRINGW:05:B:L:03',
                        'MFRAME:02:4:01','MFRAME:01:4:01']:
                if val in desc:
                    tag_gear.append(self._descriptions[desc])

        # APPLY

        self._apply_wing_body_elevon = tag_to_apply(tag_wing_body_elevon, self._elem_bdf, self._elem_tag_bdf)
        self._apply_body_flap        = tag_to_apply(tag_body_flap    , self._elem_bdf, self._elem_tag_bdf)
        self._apply_vertical_tail    = tag_to_apply(tag_vertical_tail, self._elem_bdf, self._elem_tag_bdf)


        self._apply_tps        = join(self._apply_wing_body_elevon, self._apply_body_flap)
        self._apply_equip      = self._apply_wing_body_elevon

        apply_longerons         = tag_to_apply(tag_longerons        , self._elem_bdf, self._elem_tag_bdf)
        apply_longerons_payload = tag_to_apply(tag_longerons_payload, self._elem_bdf, self._elem_tag_bdf)

        apply_avionics_frames   = tag_to_apply(tag_avionics_frames, self._elem_bdf, self._elem_tag_bdf)

        apply_lox_front_frames  = tag_to_apply(tag_lox_front_frames    , self._elem_bdf, self._elem_tag_bdf)
        apply_lox_rear_frames   = tag_to_apply(tag_lox_rear_frames     , self._elem_bdf, self._elem_tag_bdf)
        apply_kero_front_frames = tag_to_apply(tag_kero_front_frames   , self._elem_bdf, self._elem_tag_bdf)
        apply_kero_rear_frames  = tag_to_apply(tag_kero_rear_frames    , self._elem_bdf, self._elem_tag_bdf)

        apply_engine_frames    = tag_to_apply(tag_engine_frames  , self._elem_bdf, self._elem_tag_bdf)

        apply_payload_frames   = tag_to_apply(tag_payload_frames , self._elem_bdf, self._elem_tag_bdf)


        self._apply_gear       = tag_to_apply(tag_gear, self._elem_bdf, self._elem_tag_bdf)
        self._apply_hydraulic  = tag_to_apply(tag_hydraulic, self._elem_bdf, self._elem_tag_bdf)
        self._apply_avionics   = intersection(apply_avionics_frames, apply_longerons)
        self._apply_elec       = self._apply_avionics

        self._apply_lox_front  = intersection(apply_lox_front_frames     , apply_longerons)
        self._apply_lox_rear   = intersection(apply_lox_rear_frames     , apply_longerons)
        self._apply_lox        = join(self._apply_lox_front, self._apply_lox_rear)

        self._apply_kero_front = intersection(apply_kero_front_frames    , apply_longerons)
        self._apply_kero_rear  = intersection(apply_kero_rear_frames    , apply_longerons)
        self._apply_kero       = join(self._apply_kero_front, self._apply_kero_rear)

        self._apply_engine     = intersection(apply_engine_frames  , apply_longerons)

        self._apply_payload    = intersection(apply_payload_frames , apply_longerons_payload)


def isInt(s):

    try: 
        int(s)
        return True
    except ValueError:
        return False

#: def isInt()

def intersection(vec_a, vec_b):

    return [val for val in vec_a if val in vec_b]

#: def intersection()

def join(vec_a, vec_b):

    return np.unique(vec_a + vec_b).tolist()

#: def join()

def tag_to_apply(tag, elem_bdf, elem_tag_bdf):

    apply_list = []

    for iElem_bdf in range(len(elem_bdf)):
        if elem_tag_bdf[iElem_bdf] in tag:
            apply_list += (np.array(elem_bdf[iElem_bdf])-1).tolist()

    return np.unique(apply_list).tolist()

#: def tag_to_apply()

def get_apply_surface(apply_list, area_bdf):

    surface = 0.0

    for iPoint_bdf in apply_list:
        surface += area_bdf[iPoint_bdf]

    return surface

#: get_apply_surface()

def read_bdf(config, nDim, nNode):

    # Read bdf
    bdf = open(config.STRUCT + '.bdf')
    coord_bdf = []
    elem_bdf = []
    elem_tag_bdf = []
    descriptions = {}
    for line in bdf:
        data = line.split()
        if (line[0]=="$" and len(data) == 3):
            descriptions[data[2].strip().split('/')[0].upper()] = int(data[1])
        elif (line[0]=="G" and len(data) == 6):
            vec = [float(data[3]), float(data[4].strip('*G'))]
        elif (line[0]=="*" and len(data) == 5):
            vec.append(float(data[2]))
            coord_bdf.append(vec)
        elif (line[0]=="C" and len(data) == nNode+3):
            if (nNode > 3):
                elem_bdf.append([int(data[3]), int(data[4]), int(data[5]), int(data[6])])
            else:
                elem_bdf.append([int(data[3]), int(data[4]), int(data[5])])
            elem_tag_bdf.append(int(data[2]))
    bdf.close()

    return coord_bdf, elem_bdf, elem_tag_bdf, descriptions

#: def read_bdf()

def read_mesh(config, nDim, nNode):

    # Read mesh

    mesh = open(config.STRUCT + '_surface.mesh')
    line = mesh.readline()
    while not line.strip() == 'Vertices':
        line = mesh.readline()
    nPoint = int(mesh.readline())
    mesh.readline()
    coord = [[0.0 for iDim in range(nDim)] for iPoint in range(nPoint)]
    bdf_corresp = [0]*nPoint
    for iPoint in range(nPoint):
        data = mesh.readline().split()
        coord[iPoint][0] = float(data[0])
        coord[iPoint][1] = float(data[1])
        coord[iPoint][2] = float(data[2])
        bdf_corresp[iPoint] = int(data[3])-1
    line = mesh.readline()
    while not isInt(line):
        line = mesh.readline()
    nElem = int(line)
    mesh.readline()
    elem = [[0 for iNode in range(nNode)] for iElem in range(nElem)]
    for iElem in range(nElem):
        data = mesh.readline().split()
        for iNode in range(nNode):
            elem[iElem][iNode] = int(data[iNode])-1
    mesh.close()

    return coord, elem, bdf_corresp

#: def read_mesh()

def read_sol(config, nDim, nNode):

    # Read sol
    sol = open(config.STRUCT + '_surface.sol')
    line = sol.readline()
    while not line.strip() == 'SolAtVertices':
        line = sol.readline()
    nPoint = int(sol.readline())
    sol.readline()
    sol.readline()

    pressureCoeff = [0.0 for iPoint in range(nPoint)]
    frictionCoeff = [[0.0 for iDim in range(nDim)] for iPoint in range(nPoint)]

    for iPoint in range(nPoint):
        data = sol.readline().split()
        pressureCoeff[iPoint] = float(data[0])
        frictionCoeff[iPoint][0] = float(data[0+1])
        frictionCoeff[iPoint][1] = float(data[1+1])
        frictionCoeff[iPoint][2] = float(data[2+1])
    sol.close()

    return pressureCoeff, frictionCoeff

#: def read_sol()


def write_load(filename,load_bdf,coord_bdf,elem_bdf,nNode):

    nPoint_bdf = len(coord_bdf)
    nElem_bdf = len(elem_bdf)

    load = open(filename,'w')
    load.write(str(nPoint_bdf) + " " + str(nElem_bdf) + "\n")
    for iPoint_bdf in range(nPoint_bdf):
        load.write(str(coord_bdf[iPoint_bdf][0]) + " " + str(coord_bdf[iPoint_bdf][1]) + " " + str(coord_bdf[iPoint_bdf][2]) + " " + str(load_bdf[iPoint_bdf][0]) + " " + str(load_bdf[iPoint_bdf][1]) + " " + str(load_bdf[iPoint_bdf][2]) + "\n")
    if (nNode > 3):
        for iElem_bdf in range(nElem_bdf):
            load.write(str(elem_bdf[iElem_bdf][0]-1) + " " + str(elem_bdf[iElem_bdf][1]-1)  + " " + str(elem_bdf[iElem_bdf][2]-1) + " " + str(elem_bdf[iElem_bdf][3]-1) + "\n")
    else:
        for iElem_bdf in range(nElem_bdf):
            load.write(str(elem_bdf[iElem_bdf][0]-1) + " " + str(elem_bdf[iElem_bdf][1]-1)  + " " + str(elem_bdf[iElem_bdf][2]-1) + "\n")

    load.close()

#: def write_load()

def write_check(load_bdf,coord_bdf,elem_bdf,normal_voronoi,nNode):

    nPoint_bdf = len(coord_bdf)
    nElem_bdf = len(elem_bdf)

    load_mesh = open('struct_load.mesh', 'w')
    load_mesh.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nVertices\n' + str(nPoint_bdf) + '\n\n')
    for iPoint_bdf in range(nPoint_bdf):
        load_mesh.write(str(coord_bdf[iPoint_bdf][0]) + " " + str(coord_bdf[iPoint_bdf][1]) + " " + str(coord_bdf[iPoint_bdf][2]) + " " + str(iPoint_bdf+1) + "\n")
    
    if (nNode > 3):
        load_mesh.write('\nQuadrilaterals\n' + str(nElem_bdf) + '\n\n')
        for iElem_bdf in range(nElem_bdf):
            load_mesh.write(str(elem_bdf[iElem_bdf][0]) + " " + str(elem_bdf[iElem_bdf][1])  + " " + str(elem_bdf[iElem_bdf][2]) + " " + str(elem_bdf[iElem_bdf][3]) + " 0\n")
    else:
        load_mesh.write('\nTriangles\n' + str(nElem_bdf) + '\n\n')
        for iElem_bdf in range(nElem_bdf):
            load_mesh.write(str(elem_bdf[iElem_bdf][0]) + " " + str(elem_bdf[iElem_bdf][1])  + " " + str(elem_bdf[iElem_bdf][2]) + " 0\n")
    
    load_mesh.write('\nEnd\n')
    load_mesh.close()

    load_sol = open('struct_load.sol', 'w')
    load_sol.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nSolAtVertices\n' + str(nPoint_bdf) + '\n3 1 1 1\n')
    for iPoint_bdf in range(nPoint_bdf):
        load_sol.write(str(load_bdf[iPoint_bdf][0]) + " " + str(load_bdf[iPoint_bdf][1]) + " " + str(load_bdf[iPoint_bdf][2]) + "\n")
    load_sol.write('\nEnd\n')
    load_sol.close()

    load_dat = open('struct_load.dat', 'w')
    load_dat.write('TITLE = "Visualization of the surface solution"\n')
    load_dat.write('VARIABLES = "x""y""z""nx""ny""nz"\n')
    load_dat.write('ZONE NODES= %d, ELEMENTS= %d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n' % (nPoint_bdf, nElem_bdf))
    for iPoint_bdf in range(nPoint_bdf):
        load_dat.write(str(coord_bdf[iPoint_bdf][0]) + " " + str(coord_bdf[iPoint_bdf][1]) + " " + str(coord_bdf[iPoint_bdf][2]) + " " + str(normal_voronoi[iPoint_bdf][0]) + " " + str(normal_voronoi[iPoint_bdf][1]) + " " + str(normal_voronoi[iPoint_bdf][2]) + "\n")
    if (nNode > 3):
        for iElem_bdf in range(nElem_bdf):
            load_dat.write(str(elem_bdf[iElem_bdf][0]) + " " + str(elem_bdf[iElem_bdf][1])  + " " + str(elem_bdf[iElem_bdf][2]) + " " + str(elem_bdf[iElem_bdf][3]) + "\n")
    else:
        for iElem_bdf in range(nElem_bdf):
            load_dat.write(str(elem_bdf[iElem_bdf][0]) + " " + str(elem_bdf[iElem_bdf][1])  + " " + str(elem_bdf[iElem_bdf][2]) + " " + str(elem_bdf[iElem_bdf][2]) + "\n")
    load_dat.close()

#: def write_check()

def areas_normals_elem(nDim, nNode, coord_bdf, elem_bdf):

    nPoint_bdf = len(coord_bdf)
    nElem_bdf = len(elem_bdf)

    normal = [[0.0 for iDim in range(nDim)] for iElem_bdf in range(nElem_bdf)]
    area = [0.0 for iElem_bdf in range(nElem_bdf)]

    vec_a = [0.0 for iDim in range(nDim)] 
    vec_b = [0.0 for iDim in range(nDim)]

    center = [[0.0 for iDim in range(nDim)] for iElem_bdf in range(nElem_bdf)]

    # Normal and Area Elem

    for iElem_bdf in range(nElem_bdf):
        iPoint_0 = elem_bdf[iElem_bdf][0]-1
        iPoint_1 = elem_bdf[iElem_bdf][1]-1
        iPoint_2 = elem_bdf[iElem_bdf][2]-1
        if (nNode > 3): iPoint_3 = elem_bdf[iElem_bdf][3]-1
        for iDim in range(nDim):
            vec_a[iDim] = coord_bdf[iPoint_0][iDim]-coord_bdf[iPoint_1][iDim]
            vec_b[iDim] = coord_bdf[iPoint_2][iDim]-coord_bdf[iPoint_1][iDim]
            if (nNode > 3):
                center[iElem_bdf][iDim] = 0.25*(coord_bdf[iPoint_0][iDim]+coord_bdf[iPoint_1][iDim]+coord_bdf[iPoint_2][iDim]+coord_bdf[iPoint_3][iDim])
            else:
                center[iElem_bdf][iDim] = 0.333333333*(coord_bdf[iPoint_0][iDim]+coord_bdf[iPoint_1][iDim]+coord_bdf[iPoint_2][iDim])
        normal[iElem_bdf][0] += 0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])
        normal[iElem_bdf][1] += -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0])
        normal[iElem_bdf][2] += 0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])
        if (nNode > 3): 
            for iDim in range(nDim):
                vec_a[iDim] = coord_bdf[iPoint_2][iDim]-coord_bdf[iPoint_3][iDim]
                vec_b[iDim] = coord_bdf[iPoint_0][iDim]-coord_bdf[iPoint_3][iDim]
            normal[iElem_bdf][0] += 0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])
            normal[iElem_bdf][1] += -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0])
            normal[iElem_bdf][2] += 0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])

        area[iElem_bdf] += np.sqrt(normal[iElem_bdf][0]*normal[iElem_bdf][0] + normal[iElem_bdf][1]*normal[iElem_bdf][1] + normal[iElem_bdf][2]*normal[iElem_bdf][2])

    return area, normal, center

#: def areas_elem()

def areas_normals_voronoi(nDim, nNode, coord, elem):

    nPoint = len(coord)
    nElem = len(elem)

    normal = [[0.0 for iDim in range(nDim)] for iPoint in range(nPoint)]
    area = [0.0 for iPoint in range(nPoint)]

    coordElemCG = [[0.0 for iDim in range(nDim)] for iElem in range(nElem)]
    coordEdgeCG = [0.0 for iDim in range(nDim)]
    vec_a = [0.0 for iDim in range(nDim)] 
    vec_b = [0.0 for iDim in range(nDim)] 

    # coordElemCG

    for iElem in range(nElem):
        for iNode in range(nNode):
            iPoint = elem[iElem][iNode]
            for iDim in range(nDim):
                coordElemCG[iElem][iDim] += coord[iPoint][iDim];

    for iElem in range(nElem):
        for iDim in range(nDim):
            coordElemCG[iElem][iDim] /= nNode * 1.0

    # Normal and Area Voronoi

    nNeighnour = 2
    for iElem in range(nElem):
        for iNode in range(nNode):
            iPoint = elem[iElem][iNode]
            for iNeighbour in range(nNeighnour):
                jNode = (iNode + 1 - nNeighnour * iNeighbour) % nNode
                jPoint = elem[iElem][jNode]
                for iDim in range(nDim):
                    coordEdgeCG[iDim] = 0.5 * (coord[iPoint][iDim] + coord[jPoint][iDim])
                if (iNeighbour == 0):
                    for iDim in range(nDim):
                        vec_a[iDim] = coord[iPoint][iDim]-coordElemCG[iElem][iDim]
                        vec_b[iDim] = coordEdgeCG[iDim]-coordElemCG[iElem][iDim]
                else:
                    for iDim in range(nDim):
                        vec_a[iDim] = coord[iPoint][iDim]-coordEdgeCG[iDim]
                        vec_b[iDim] = coordElemCG[iElem][iDim]-coordEdgeCG[iDim]
                normal[iPoint][0] += 0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])
                normal[iPoint][1] += -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0])
                normal[iPoint][2] += 0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])

    for iPoint in range(nPoint):
        area[iPoint] = (normal[iPoint][0]**2.0+normal[iPoint][1]**2.0+normal[iPoint][2]**2.0)**0.5

    return area, normal

#: def areas_normals()



