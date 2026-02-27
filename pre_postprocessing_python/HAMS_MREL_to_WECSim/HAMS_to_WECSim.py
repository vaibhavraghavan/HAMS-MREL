# This script can be used to convert the results of HAMS-MREL to the .out file of WAMIT, which can then be utilized as
# input in WEC-Sim, currently only implemented for multi-bodies
# The script is created by Vaibhav Raghavan, MSc TU Delft; Copyright TU Delft 2024
import math
import cmath
import os

import numpy as np

from postprocessing_functions_Mat import (hydrodynamic_coefficients_exciting_forces_all_bodies,
                                          obtain_specific_exciting_force_matrix,
                                          obtain_specific_coefficient_matrix,
                                          hydrodynamic_coefficients_exciting_forces_all_bodies_zero_inf)
from datetime import date, datetime

class HAMS_to_WAMIT():
    """
    This class is for converting the results from HAMS-MREL to the .out for WAMIT. This can be utilized for calculations
    in WEC-Sim
    """

    def __init__(self, filename, hams_result_folder, hams_input_folder, rho=1000):
        self.filename = filename
        self.hams_result_folder = hams_result_folder
        self.hams_input_folder = hams_input_folder
        self.rho = rho

    def create_initial_wamit_output(self):
        today = date.today()
        now = datetime.now()
        self.date = today.strftime("%d %b %Y")
        self.time = now.strftime("%H:%M:%S")
        self.output = open(f'{self.filename}.out','w')
        self.output.write(' -----------------------------------------------------------------------\n')
        self.output.write('\n')
        self.output.write('                        WAMIT  Version x.x\n')
        self.output.write('\n')
        self.output.write('     Copyright (c) 1999-2016 WAMIT Incorporated\n')
        self.output.write('     Copyright (c) 1998 Massachusetts Institute of Technology\n')
        self.output.write(' -----------------------------------------------------------------------\n')
        self.output.write('\n')
        self.output.write('     The WAMIT software performs computations of wave interactions with\n')
        self.output.write('     floating or submerged vessels.  WAMIT is a registered trademark of\n')
        self.output.write('     WAMIT Incorporated.  The results of HAMS-MREL are converted to WAMIT by\n')
        self.output.write('\n')
        self.output.write('            Marine Renewable Energies Lab\n')
        self.output.write('            CiTG, Stevinweg 1\n')
        self.output.write('            2628 CN, Delft, Netherlands\n')
        self.output.write(f'     for use in WEC-Sim.                      Release date: {self.date}')
        self.output.write('\n')
        self.output.write(' -----------------------------------------------------------------------\n')
        self.output.write('\n')
        self.output.write('\n')
    def get_number_of_bodies(self):
        f = open(self.hams_input_folder + '/ControlFile.in')
        lines = [line.strip().split() for line in f]
        self.n_bodies = int([line[1] for line in lines if 'Number_of_bodies' in line][0])
        f.close()

    def obtain_wavenumber_shallow_water(self, period, depth, g = 9.81):
        """
        This function calculates the wave length of regular waves in shallow water

        :param period: Wave period in s
        :param depth: Water depth in m
        :param g: gravitational constant in m/s^2
        :return: Wavelength in m
        """

        L_0 = (g / (2 * math.pi)) * period ** 2
        k = 2 * math.pi / L_0
        w = 2 * math.pi / period
        # print(f"Iterative process started from L_0 = {L_0}")

        w_2 = k * g * math.tanh(k * depth)
        error = w ** 2 - w_2

        L_0_final = L_0 # Assume this to be the default value, and then the iteration process can take over
        while error > 0.001:
            w_2 = k * g * math.tanh(k * depth)
            L_0_final = L_0
            error = w ** 2 - w_2
            if error > 0:
                L_0 = L_0 - L_0 / 2
            else:
                L_0 = L_0 + L_0 / 2
            k = 2 * math.pi / L_0
            error = abs(error)
        return k

    def get_frequency_information(self, g = 9.81):
        f = open(self.hams_input_folder + '/ControlFile.in')
        lines = [line.strip().split() for line in f]
        self.period_print = []
        self.wavenumber_print = []
        self.omega_print = []
        check_zero_infinite = int([line[1] for line in lines if '0_inf_frequency_limits' in line][0])
        if check_zero_infinite == 1:
            self.period_print.append(-1)
            self.period_print.append(0)
            self.wavenumber_print.append('zero')
            self.wavenumber_print.append('infinite')
            self.omega_print.append(0)
            self.omega_print.append(-1)
        self.input_frequency_type = int([line[1] for line in lines if 'Input_frequency_type' in line][0])
        self.output_frequency_type = int([line[1] for line in lines if 'Output_frequency_type' in line][0])
        self.n_frequencies = abs(int([line[1] for line in lines if 'Number_of_frequencies' in line][0]))
        self.minimum_frequency = float([line[1] for line in lines if 'Minimum_frequency_Wmin' in line][0].split('D0')[0])
        self.incidence_angle = int([lines[i+1][0] for i in range(len(lines)) if 'Number_of_headings' in lines[i]][0].split('.D0')[0])
        self.frequency_step = float([line[1] for line in lines if 'Frequency_step' in line][0].split('D0')[0])
        self.length_scale = float([line[1] for line in lines if 'Reference_body_length' in line][0].split('.D0')[0])
        self.water_depth = float([line[1] for line in lines if 'Waterdepth' in line][0].split('.D0')[0])
        self.body_center = []
        self.rotation_center = []
        for i in range(self.n_bodies):
            x = float([line[1] for line in lines if f'LCS_{i + 1}' in line][0].split('.D0')[0])
            y = float([line[2] for line in lines if f'LCS_{i + 1}' in line][0].split('.D0')[0])
            z = float([line[3] for line in lines if f'LCS_{i + 1}' in line][0].split('.D0')[0])
            self.body_center.append([x,y,z])
            xb = float([line[1] for line in lines if f'Reference_body_center_{i + 1}' in line][0].split('.D0')[0])
            yb = float([line[2] for line in lines if f'Reference_body_center_{i + 1}' in line][0].split('.D0')[0])
            zb = float([line[3] for line in lines if f'Reference_body_center_{i + 1}' in line][0].split('.D0')[0])
            self.rotation_center.append([xb, yb, zb])
        f.close()
        for i in range(self.n_frequencies):
            if self.input_frequency_type == 1:
                wavenumber = self.minimum_frequency + i * self.frequency_step
                omega = math.sqrt(wavenumber*g)
                period = 2*math.pi/(omega)
            elif self.input_frequency_type == 2:
                wavenumber = self.minimum_frequency + i * self.frequency_step
                omega = math.sqrt(wavenumber*g*np.tanh(wavenumber*self.water_depth))
                period = 2 * math.pi / (omega)
            elif self.input_frequency_type == 3:
                omega = self.minimum_frequency + i * self.frequency_step
                period = 2*math.pi/(omega)
                if self.water_depth > 0:
                    wavenumber = self.obtain_wavenumber_shallow_water(period=period, depth=self.water_depth)
                else:
                    wavenumber = omega**2 / g
            elif self.input_frequency_type == 4: # This is only for deepwaters at the moment
                period = self.minimum_frequency + i*self.frequency_step
                omega = 2*math.pi/(period)
                if self.water_depth > 0:
                    wavenumber = self.obtain_wavenumber_shallow_water(period=period, depth=self.water_depth)
                else:
                    wavenumber = omega**2 / g
            else:
                raise ValueError('This output type in currently not implemented. Please use one of the other inputs')
            self.omega_print.append(omega)
            self.period_print.append(period)
            self.wavenumber_print.append(wavenumber)

        self.number_of_panels = []
        for j in range(self.n_bodies):
            f=open(self.hams_input_folder + f'/HullMesh_{j+1}.pnl') #
            lines = [line.strip().split() for line in f]
            self.number_of_panels.append(int([lines[i+1][0] for i in range(len(lines)) if 'X-Symmetry' in lines[i]][0]))
            f.close()

    def hams_to_wamit_coefficient_conversion(self,body_number_1, body_number_2, dof_1, dof_2):
        return 6*(body_number_1-1) + dof_1, 6*(body_number_2-1) + dof_2

    def hams_to_wamit_coefficient_conversion_excitation(self,body_number_1, dof_1):
        return 6*(body_number_1-1) + dof_1

    def create_hams_compiled_results_for_wamit(self, hams_output_location, run):
        if run:
            results_full_dir = 'Results_full_HAMS'
            if not os.path.exists(results_full_dir):
                os.makedirs(results_full_dir)

            # Hydrodynamic coefficients
            for dof_1 in range(6):
                for dof_2 in range(6):
                    added_mass_file_location = hams_output_location + f'\OAMASS{dof_1 + 1}.txt'
                    added_damping_file_location = hams_output_location + f'\ODAMPING{dof_1 + 1}.txt'
                    external_force_file_location = hams_output_location + f'\OEXFOR{dof_1 + 1}.txt'
                    n_bodies, omega_final, added_mass_compiled, radiation_damping_compiled, exciting_force_compiled = (
                        hydrodynamic_coefficients_exciting_forces_all_bodies(
                            added_mass_file_location=added_mass_file_location,
                            added_damping_file_location=added_damping_file_location,#
                            external_force_file_location=external_force_file_location,
                            dof=dof_2 + 1, zero_inf_frequency=True))
                    for k in range(self.n_bodies):
                        for l in range(self.n_bodies):
                            a_list, r_list = obtain_specific_coefficient_matrix(added_mass_compiled=added_mass_compiled,
                                                                                radiation_damping_compiled=radiation_damping_compiled,
                                                                                n_bodies=n_bodies,
                                                                                x=k + 1, y=l + 1)
                            dof_final_1, dof_final_2 = self.hams_to_wamit_coefficient_conversion(k + 1, l + 1, dof_1 + 1,
                                                                                            dof_2 + 1)
                            f = open(results_full_dir + f'/hydrodynamic_coefficients_{dof_final_1}_{dof_final_2}.out', # The values for WAMIT need to be fixed, check conversion factor
                                     'w')  # This means - Body number 1, Body 2, DOF of Body 1, DOF of Body 2
                            for m in range(len(a_list)):
                                f.write(f'{a_list[m]} {r_list[m]}\n')
                            f.close()
                    self.omega = omega_final

                for n in range(self.n_bodies):
                    e_list = obtain_specific_exciting_force_matrix(exciting_force_compiled=exciting_force_compiled,
                                                                   x=n+1)
                    dof_final = self.hams_to_wamit_coefficient_conversion_excitation(n+1, dof_1+1)
                    f = open(results_full_dir + f'/excitation_forces_{dof_final}.out',
                             'w')
                    for p in range(len(e_list)):
                        f.write(f'{abs(e_list[p])} {cmath.phase(e_list[p])}\n') # The values for WAMIT need to be fixed, check conversion factor
                    f.close()

    def create_hams_compiled_results_for_wamit_zero_inf(self, hams_output_location, run):
        if run:
            results_full_dir = 'Results_full_HAMS_inf_zero'
            if not os.path.exists(results_full_dir):
                os.makedirs(results_full_dir)

            # Hydrodynamic coefficients
            for dof_1 in range(6):
                for dof_2 in range(6):
                    added_mass_file_location = hams_output_location + f'\OAMASS{dof_1 + 1}.txt'
                    added_damping_file_location = hams_output_location + f'\ODAMPING{dof_1 + 1}.txt'
                    external_force_file_location = hams_output_location + f'\OEXFOR{dof_1 + 1}.txt'
                    n_bodies, omega_final, added_mass_compiled, radiation_damping_compiled, exciting_force_compiled = (
                        hydrodynamic_coefficients_exciting_forces_all_bodies_zero_inf(
                            added_mass_file_location=added_mass_file_location,
                            added_damping_file_location=added_damping_file_location,
                            external_force_file_location=external_force_file_location,
                            dof=dof_2 + 1))
                    for k in range(self.n_bodies):
                        for l in range(self.n_bodies):
                            a_list, r_list = obtain_specific_coefficient_matrix(added_mass_compiled=added_mass_compiled,
                                                                                radiation_damping_compiled=radiation_damping_compiled,
                                                                                n_bodies=n_bodies,
                                                                                x=k + 1, y=l + 1)
                            dof_final_1, dof_final_2 = self.hams_to_wamit_coefficient_conversion(k + 1, l + 1, dof_1 + 1,
                                                                                            dof_2 + 1)
                            f = open(results_full_dir + f'/hydrodynamic_coefficients_{dof_final_1}_{dof_final_2}.out', # The values for WAMIT need to be fixed, check conversion factor
                                     'w')  # This means - Body number 1, Body 2, DOF of Body 1, DOF of Body 2
                            for m in range(len(a_list)):
                                f.write(f'{a_list[m]}\n')
                            f.close()

    def compile_hams_results_per_frequency_for_wamit(self, run):
        if run:
            results_full_dir = 'Results_frequency_HAMS'
            if not os.path.exists(results_full_dir):
                os.makedirs(results_full_dir)

            for i in range(len(self.omega)):
                print(f'Writing data for frequency {i+1}/{len(self.omega)}')
                m = open(results_full_dir + f'/Coefficients_frequency_{i + 1}', 'w')
                n = open(results_full_dir + f'/Excitation_frequency_{i + 1}', 'w')
                for j in range(self.n_bodies*6):
                    for k in range(self.n_bodies*6):
                        f = open( f'Results_full_HAMS/hydrodynamic_coefficients_{j+1}_{k+1}.out','r')
                        data = [line.strip().split() for line in f]
                        m.write(f'{j+1} {k+1} {data[i][0]} {data[i][1]}\n')
                        f.close()
                for l in range(self.n_bodies*6):
                    f = open(f'Results_full_HAMS/excitation_forces_{l + 1}.out', 'r')
                    data = [line.strip().split() for line in f]
                    n.write(f'{l + 1} {data[i][0]} {data[i][1]}\n')
                    f.close()
                print(f'Finished writing data for frequency {i+1}/{len(self.omega)}')
                m.close()
                n.close()

    def compile_hams_results_per_frequency_for_wamit_zero_inf(self, run):
        if run:
            results_full_dir = 'Results_frequency_HAMS_inf_zero'
            if not os.path.exists(results_full_dir):
                os.makedirs(results_full_dir)

            for i in range(2):
                print(f'Writing data for frequency {i+1}/{2}')
                m = open(results_full_dir + f'/Coefficients_frequency_{i + 1}', 'w')
                for j in range(self.n_bodies*6):
                    for k in range(self.n_bodies*6):
                        f = open( f'Results_full_HAMS/hydrodynamic_coefficients_{j+1}_{k+1}.out','r')
                        data = [line.strip().split() for line in f]
                        m.write(f'{j+1} {k+1} {data[i][0]}\n')
                        f.close()
                print(f'Finished writing data for frequency {i+1}/{2}')
                m.close()
    def add_input_file_description(self):
        self.output.write(' Low-order panel method  (ILOWHI=0)\n')
        self.output.write('\n')
        self.output.write(' Input from Geometric Data Files:\n')
        for i in range(self.n_bodies):
            self.output.write(f'                               N=  {i+1}     HullMesh_{i+1}.pnl\n')
            self.output.write(' Rhino->WAMIT file export (mesh)\n')
        self.output.write('\n')
        self.output.write(' Input from Potential Control File: no pot file used in HAMS-MREL\n')
        self.output.write(' No pot -- file type .pnl, ILOWHI=0, IRR=0\n')
        self.output.write('\n')
        self.output.write('\n')

    def add_frequency_information(self):
        self.output.write(f' POTEN run date and starting time:        {self.date}  --  {self.time} \n')
        self.output.write(f'   Period       Time           RAD      DIFF  (max iterations)\n')
        for i in range(len(self.period_print)):
            if self.period_print[i] == -1:
                now = datetime.now()
                self.time = now.strftime("%H:%M:%S")
                self.output.write(f'   {self.period_print[i]:.4f}    {self.time}          -1          \n')
            elif self.period_print[i] == 0:
                now = datetime.now()
                self.time = now.strftime("%H:%M:%S")
                self.output.write(f'    {self.period_print[i]:.4f}    {self.time}          -1          \n')
            else:
                now = datetime.now()
                self.time = now.strftime("%H:%M:%S")
                self.output.write(f'    {self.period_print[i]:.4f}    {self.time}          -1      -1 \n')
        self.output.write(f'\n')
        self.output.write(f' Gravity:     9.80665                Length scale:        {self.length_scale:.5f}\n')
        self.output.write(f' Water depth:        {self.water_depth:.5f}\n')
        self.output.write(f' Logarithmic singularity index:              ILOG =     1\n')
        self.output.write(f' Source formulation index:                   ISOR =     0\n')
        self.output.write(f' Diffraction/scattering formulation index: ISCATT =     0\n')
        self.output.write(f' Number of blocks used in linear system:   ISOLVE =     1\n')
        self.output.write(f' Number of unknowns in linear system:        NEQN =  {sum(self.number_of_panels)}\n')
        self.output.write(f' Logarithmic singularity index:              ILOG =     1\n')
        self.output.write(f' Global Symmetries: none\n')

    def add_hydrostatic_data(self, volumes, hydrostatic, center_of_gravity, center_of_bouyancy, g = 9.81):
        self.output.write(f'\n')
        self.output.write(f' BODY PARAMETERS:\n')
        self.output.write(f'\n')
        for i in range(self.n_bodies):
            self.output.write(f' Body number: N= {i+1}       Total panels:   {self.number_of_panels[i]}    Waterline panels:   32\n') # The waterline panels are currently a random number
            self.output.write(
                f'               Irregular frequency index: IRR =0      Symmetries: none\n')
            self.output.write(f'\n')
            self.output.write(f' XBODY =   {self.body_center[i][0]:.4f} YBODY =    {self.body_center[i][1]:.4f} ZBODY =    {self.body_center[i][2]:.4f} PHIBODY =   0.0\n')
            self.output.write(
                f' Volumes (VOLX,VOLY,VOLZ):           {volumes[i][0]:.3f}      {volumes[i][1]:.3f}      {volumes[i][2]:.3f}\n')
            self.output.write(f' Center of Buoyancy (Xb,Yb,Zb):     {center_of_bouyancy[i][0]:.6f}     {center_of_bouyancy[i][1]:.6f}    {center_of_bouyancy[i][2]:.6f}\n') # Currently rotation center is added. This can be modified later. This is also based on LCS?
            self.output.write(f' Hydrostatic and gravitational restoring coefficients:\n')
            self.output.write(f' C(3,3),C(3,4),C(3,5):   {hydrostatic[i][2][2]/(self.rho * g):.3f}     {hydrostatic[i][2][3]/(self.rho * g):.3f}     {hydrostatic[i][2][4]/(self.rho * g):.3f}\n')
            self.output.write(
                f' C(4,4),C(4,5),C(4,6):                {hydrostatic[i][3][3]/(self.rho * g):.3f}     {hydrostatic[i][3][4]/(self.rho * g):.3f}     {hydrostatic[i][3][4]/(self.rho * g):.3f}\n')
            self.output.write(
                f'        C(5,5),C(5,6):                             {hydrostatic[i][4][4]/(self.rho * g):.3f}     {hydrostatic[i][4][5]/(self.rho * g):.3f}\n')
            self.output.write(f' Center of Gravity  (Xg,Yg,Zg):     {center_of_gravity[i][0]:.6f}     {center_of_gravity[i][1]:.6f}     {center_of_gravity[i][2]:.6f}\n') # This should be the local COG?
            self.output.write(' Radii of gyration:     0.000000     0.000000     0.000000\n')
            self.output.write('                        0.000000     0.000000     0.000000\n')
            self.output.write('                        0.000000     0.000000     0.000000\n')
            self.output.write(f'\n')
        self.output.write(f'\n')
        self.output.write(f' ------------------------------------------------------------------------\n')

    def add_hydrodynamic_coefficients_per_frequency(self, omega_number):
        results_full_dir = 'Results_frequency_HAMS'
        self.output.write(f' ************************************************************************\n')
        self.output.write(f'\n')
        if type(self.period_print[omega_number]) == float:
            self.output.write(
                f' Wave period (sec) =  {float(self.period_print[omega_number]):.6E}        Wavenumber (kL) =   {float(self.wavenumber_print[omega_number]):.6E}\n')
        else:
            self.output.write(
                f' Wave period (sec) =  {self.period_print[omega_number]}        Wavenumber (kL) =   {self.wavenumber_print[omega_number]}\n')
        self.output.write(f' ------------------------------------------------------------------------\n')
        self.output.write(f'\n')
        self.output.write(f'\n')
        self.output.write(f'    ADDED-MASS AND DAMPING COEFFICIENTS\n')
        self.output.write(f'     I     J         A(I,J)         B(I,J)\n')
        self.output.write(f'\n')
        f = open(results_full_dir + f'/Coefficients_frequency_{omega_number-1}', 'r')
        data = [line.strip().split() for line in f]
        for i in range(len(data)):
            self.output.write(f'     {data[i][0]}     {data[i][1]}   {float(data[i][2])/self.rho:.6E}   {float(data[i][3])/(self.rho*self.omega_print[omega_number]):.6E}\n') # The quantities are made dimensionless here
        self.output.write(f'\n')
        self.output.write(f'\n')
        self.output.write(f'\n')

    def add_initial_output_info(self):
        self.output.write(f'                            Output from  WAMIT\n') # Actually from HAMS
        self.output.write(f' ------------------------------------------------------------------------\n')
        self.output.write(f' FORCE run date and starting time:                {self.date} -- {self.time}\n')
        self.output.write(f' ------------------------------------------------------------------------\n')
        self.output.write(f' I/O Files:              wamit_hams_output.out\n') #Only the output file is displayed for now
        self.output.write(f'  wamit_hams_output.out -- file type .pnl, ILOWHI=0, IRR=0\n')
        self.output.write(f'\n')
        self.output.write(f'\n')

    def add_added_mass_zero_inf(self):
        results_full_dir = 'Results_frequency_HAMS_inf_zero'
        ################For infinite period################################################################
        self.output.write(f' ************************************************************************\n')
        self.output.write(f'\n')
        self.output.write(
                f' Wave period (sec) =  infinite        Wavenumber (kL) =   zero\n')
        self.output.write(f' ------------------------------------------------------------------------\n')
        self.output.write(f'\n')
        self.output.write(f'\n')
        self.output.write(f'    ADDED-MASS COEFFICIENTS\n')
        self.output.write(f'     I     J         A(I,J)\n')
        self.output.write(f'\n')
        f = open(results_full_dir + f'/Coefficients_frequency_1', 'r')
        data = [line.strip().split() for line in f]
        for i in range(len(data)):
            self.output.write(f'     {data[i][0]}     {data[i][1]}   {float(data[i][2])/self.rho:.6E}\n') # The quantities are made dimensionless here
        self.output.write(f'\n')
        self.output.write(f'\n')
        self.output.write(f'\n')

        ################For zero period################################################################
        self.output.write(f' ************************************************************************\n')
        self.output.write(f'\n')
        self.output.write(
            f' Wave period (sec) =  zero        Wavenumber (kL) =   infinite\n')
        self.output.write(f' ------------------------------------------------------------------------\n')
        self.output.write(f'\n')
        self.output.write(f'\n')
        self.output.write(f'    ADDED-MASS COEFFICIENTS\n')
        self.output.write(f'     I     J         A(I,J)\n')
        self.output.write(f'\n')
        f = open(results_full_dir + f'/Coefficients_frequency_2', 'r')
        data = [line.strip().split() for line in f]
        for i in range(len(data)):
            self.output.write(
                f'     {data[i][0]}     {data[i][1]}   {float(data[i][2]) / self.rho:.6E}\n')  # The quantities are made dimensionless here
        self.output.write(f'\n')
        self.output.write(f'\n')
        self.output.write(f'\n')

    def add_excitation_forces_per_frequency(self, omega_number, g = 9.81):
        results_full_dir = 'Results_frequency_HAMS'
        f = open(results_full_dir + f'/Excitation_frequency_{omega_number - 1}', 'r')
        data = [line.strip().split() for line in f]
        self.output.write(f'\n')
        self.output.write(f'    DIFFRACTION EXCITING FORCES AND MOMENTS\n')
        self.output.write(f'\n')
        self.output.write(f'  Wave Heading (deg) :      {self.incidence_angle}\n')
        self.output.write(f'\n')
        self.output.write(f'     I     Mod[Xh(I)]     Pha[Xh(I)]\n')
        self.output.write(f'\n')
        for i in range(len(data)):
            self.output.write(f'     {data[i][0]}     {float(data[i][1])/(self.rho * g):.6E}   {int(float(data[i][2])*180/math.pi)}\n')
        self.output.write(f'\n')
        self.output.write(f'\n')
    def add_hydrodynamic_coefficients_excitation_forces(self):
        for i in range(2, len(self.period_print)):
            self.add_hydrodynamic_coefficients_per_frequency(omega_number=i)
            self.add_excitation_forces_per_frequency(omega_number=i)

# Add the information about the zero and infinite frequency added mass
# Addition of the excitation force coefficients - Added
# Making the quantities dimensionless - Added
# Obtaining the correct hydrostatic quantities (Is this where it is taken from by WEC Sim?) - Added
# The right wave number and time period needs to be added - Added
# Check the center of gravity input - Added
# Are the volumes relevant for WECSim?

########################################################################################################################
# Running the script
# NOTE: Reference length is always 1 here

create_data_from_HAMS = False
create_final_WECSIM_files = True

if create_data_from_HAMS:
    hams_input_folder = 'HAMS_case/Input'
    output_file = HAMS_to_WAMIT(filename='wamit_hams_output', hams_input_folder=hams_input_folder,
                                hams_result_folder='HAMS_case')
    output_file.get_number_of_bodies()
    output_file.get_frequency_information()
    output_file.create_initial_wamit_output()
    output_file.add_input_file_description()
    output_file.add_frequency_information()
    # Only needs to be run once for creating the results
    # Getting the input from the HAMS-MREL output
    hams_output_location = r'UGroningen_OceanGrazer_full\Output\Hams_format'

    # Compiling the added mass, radiation damping matrices and excitation force, and creating the results for all DOFs
    output_file.create_hams_compiled_results_for_wamit(hams_output_location=hams_output_location, run=True)

    # Convert the compiled files into a format can be used to easily create the output for WAMIT
    output_file.compile_hams_results_per_frequency_for_wamit(run=True)

    # Compiling the added mass, and creating the results for all DOFs for zero and inf frequency
    output_file.create_hams_compiled_results_for_wamit_zero_inf(hams_output_location=hams_output_location, run=True)

    # Convert the compiled files into a format can be used to easily create the output for WAMIT
    output_file.compile_hams_results_per_frequency_for_wamit_zero_inf(run=True)


if create_final_WECSIM_files:
    hams_input_folder = 'HAMS_case/Input'
    output_file = HAMS_to_WAMIT(filename='wamit_hams_output', hams_input_folder = hams_input_folder,
                                hams_result_folder='HAMS_case')
    output_file.get_number_of_bodies()
    output_file.get_frequency_information()
    output_file.create_initial_wamit_output()
    output_file.add_input_file_description()
    output_file.add_frequency_information()
    volumes = [[1238.37,1238.37,1238.37],[1238.37,1238.37,1238.37],[1238.37,1238.37,1238.37]]
    hydrostatic = [[0,0,0,0,0,0],
                   [0,0,0,0,0,0],
                   [0,0,0.27275E+06,0,0,0],
                   [0,0,0,-0.15754E+08,0,0],
                   [0,0,0,0,-0.15754E+08,0],
                   [0,0,0,0,0,0]]
    hydrostatic_final = [hydrostatic, hydrostatic, hydrostatic]
    center_of_gravity = [[0,0,0],[0,0,0],[0,0,0]]
    center_of_bouyancy = [[0, 0, -5.475], [0, 0, -5.475], [0, 0, -5.475]]
    output_file.add_hydrostatic_data(volumes=volumes, hydrostatic=hydrostatic_final, center_of_gravity=center_of_gravity,
                                     center_of_bouyancy=center_of_bouyancy)
    output_file.add_initial_output_info()
    output_file.add_added_mass_zero_inf()
    output_file.add_hydrodynamic_coefficients_excitation_forces()







