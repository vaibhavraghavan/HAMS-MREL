import math
import cmath
import os

import numpy as np

from post_processing_HAMS_multibodies import (hydrodynamic_coefficients_exciting_forces_all_bodies,
                                                        obtain_specific_exciting_force_matrix,
                                                        obtain_specific_coefficient_matrix,
                                                        hydrodynamic_coefficients_exciting_forces_all_bodies_zero_inf)
from datetime import date, datetime

class Multibodies_HAMS():
    """
    This class is for post-processing the results of HAMS-MREL specifically for multiple bodies
    """

    def __init__(self, filename, hams_result_folder, hams_input_folder, rho=1025):
        self.filename = filename
        self.hams_result_folder = hams_result_folder
        self.hams_input_folder = hams_input_folder
        self.rho = rho
    def get_number_of_bodies(self):
        """
        This function obtains the number of bodies in the simulation
        """
        f = open(self.hams_input_folder + '/ControlFile.in')
        lines = [line.strip().split() for line in f]
        self.n_bodies = int([line[1] for line in lines if 'Number_of_bodies' in line][0])
        f.close()

    def obtain_wavenumber_shallow_water(self, period, depth, g = 9.81):
        """
        This function calculates the wave length of regular waves in shallow water
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
        """
        This function obtains the all the information from the Control.in file for postprocessing
        """
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
        self.minimum_frequency = float([line[1] for line in lines if 'Minimum_frequency_Wmin' in line][0].split('.D0')[0])
        self.n_directions = int([lines[i][1] for i in range(len(lines)) if 'Number_of_headings' in lines[i]][0].split('.D0')[0])
        temp_list = [lines[i+1] for i in range(len(lines)) if 'Number_of_headings' in lines[i]][0]
        self.incidence_angle_list = [int(ele.split('.D0')[0]) for ele in temp_list]
        self.frequency_step = float([line[1] for line in lines if 'Frequency_step' in line][0].split('D0')[0])
        self.length_scale = float([line[1] for line in lines if 'Reference_body_length' in line][0].split('D0')[0])
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

    def hams_coefficient_conversion(self, body_number_1, body_number_2, dof_1, dof_2):
        return 6*(body_number_1-1) + dof_1, 6*(body_number_2-1) + dof_2

    def hams_coefficient_conversion_excitation(self, body_number_1, dof_1):
        return 6*(body_number_1-1) + dof_1

    def create_hams_compiled_results(self, hams_output_location, run):
        """
        This function writes the compiled results from HAMS-MREL for the hydrodynamic coefficients and exciting forces
        into text files. Can be used if the user wants this data for postprocessing
        """
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
                    n_bodies, omega_final, added_mass_compiled, radiation_damping_compiled, exciting_force_compiled,\
                        n_directions= (
                        hydrodynamic_coefficients_exciting_forces_all_bodies(
                            added_mass_file_location=added_mass_file_location,
                            added_damping_file_location=added_damping_file_location,#
                            external_force_file_location=external_force_file_location,
                            dof=dof_2 + 1, zero_inf_frequency=True))
                    for k in range(self.n_bodies):
                        for l in range(self.n_bodies):
                            a_list, r_list = obtain_specific_coefficient_matrix(added_mass_compiled=added_mass_compiled,
                                                                                radiation_damping_compiled=radiation_damping_compiled,
                                                                                n_bodies=self.n_bodies,
                                                                                x=k + 1, y=l + 1)
                            dof_final_1, dof_final_2 = self.hams_coefficient_conversion(k + 1, l + 1, dof_1 + 1,
                                                                                        dof_2 + 1)
                            f = open(results_full_dir + f'/hydrodynamic_coefficients_{dof_final_1}_{dof_final_2}.out', # The values for WAMIT need to be fixed, check conversion factor
                                     'w')  # This means - Body number 1, Body 2, DOF of Body 1, DOF of Body 2
                            for m in range(len(a_list)):
                                f.write(f'{a_list[m]} {r_list[m]}\n')
                            f.close()
                    self.omega = omega_final

                for n in range(self.n_bodies):
                    for p in range(self.n_directions):
                        e_list = obtain_specific_exciting_force_matrix(exciting_force_compiled=exciting_force_compiled,
                                                                       x=n+1, j=p+1, n_directions=self.n_directions)
                        dof_final = self.hams_coefficient_conversion_excitation(n + 1, dof_1 + 1)
                        f = open(results_full_dir + f'/excitation_forces_{dof_final}_dir_{p+1}.out',
                                 'w')
                        for p in range(len(e_list)):
                            if math.isnan(e_list[p].real) or math.isnan(e_list[p].imag):
                                f.write(f'{0.0} {0.0}\n')
                            else:
                                f.write(f'{-e_list[p].imag} {-e_list[p].real}\n')
                        f.close()

    def create_hams_compiled_results_zero_inf(self, hams_output_location, run):
        """
        This function writes the compiled results from HAMS-MREL for the hydrodynamic coefficients and exciting forces
        into text files for zero and infinite frequencies only. Can be used if the used wants this data for postprocessing.
        """
        if run:
            results_full_dir = 'Results_full_HAMS_inf_zero'
            if not os.path.exists(results_full_dir):
                os.makedirs(results_full_dir)

            # Hydrodynamic coefficients
            for dof_1 in range(6):
                for dof_2 in range(6):
                    added_mass_file_location = hams_output_location + f'\OAMASS{dof_1 + 1}.txt'
                    added_damping_file_location = hams_output_location + f'\ODAMPING{dof_1 + 1}.txt'
                    n_bodies, omega_final, added_mass_compiled, radiation_damping_compiled= (
                        hydrodynamic_coefficients_exciting_forces_all_bodies_zero_inf(
                            added_mass_file_location=added_mass_file_location,
                            added_damping_file_location=added_damping_file_location,
                            dof=dof_2 + 1))
                    for k in range(self.n_bodies):
                        for l in range(self.n_bodies):
                            a_list, r_list = obtain_specific_coefficient_matrix(added_mass_compiled=added_mass_compiled,
                                                                                radiation_damping_compiled=radiation_damping_compiled,
                                                                                n_bodies=self.n_bodies,
                                                                                x=k + 1, y=l + 1)
                            dof_final_1, dof_final_2 = self.hams_coefficient_conversion(k + 1, l + 1, dof_1 + 1,
                                                                                        dof_2 + 1)
                            f = open(results_full_dir + f'/hydrodynamic_coefficients_{dof_final_1}_{dof_final_2}.out', # The values for WAMIT need to be fixed, check conversion factor
                                     'w')  # This means - Body number 1, Body 2, DOF of Body 1, DOF of Body 2
                            for m in range(len(a_list)):
                                f.write(f'{a_list[m]}\n')
                            f.close()
                            
    def compile_hams_results_per_frequency(self, run):
        """
        This function writes the compiled results from HAMS-MREL per frequency for the hydrodynamic coefficients and exciting forces
        into text files. Can be used if the used wants this data for postprocessing.
        """
        if run:
            results_full_dir = 'Results_frequency_HAMS'
            if not os.path.exists(results_full_dir):
                os.makedirs(results_full_dir)

            for i in range(len(self.omega)):
                print(f'Writing data for frequency {i+1}/{len(self.omega)}')
                m = open(results_full_dir + f'/Coefficients_frequency_{i + 1}.out', 'w')
                n = open(results_full_dir + f'/Excitation_frequency_{i + 1}.out', 'w')
                for j in range(self.n_bodies*6):
                    for k in range(self.n_bodies*6):
                        f = open( f'Results_full_HAMS/hydrodynamic_coefficients_{j+1}_{k+1}.out','r')
                        data = [line.strip().split() for line in f]
                        m.write(f'{j+1} {k+1} {data[i][0]} {data[i][1]}\n')
                        f.close()
                for l in range(self.n_bodies*6):
                    for p in range(self.n_directions):
                        f = open(f'Results_full_HAMS/excitation_forces_{l + 1}_dir_{p+1}.out', 'r')
                        data = [line.strip().split() for line in f]
                        n.write(f'{l + 1} {self.incidence_angle_list[p]} {data[i][0]} {data[i][1]}\n')
                        f.close()
                print(f'Finished writing data for frequency {i+1}/{len(self.omega)}')
                m.close()
                n.close()

    def compile_hams_results_per_frequency_zero_inf(self, run):
        """
        This function writes the compiled results from HAMS-MREL per frequency for the hydrodynamic coefficients and exciting forces
        into text files for only the zero and infinite frequencies. Can be used if the used wants this data for postprocessing.
        """
        if run:
            results_full_dir = 'Results_frequency_HAMS_inf_zero'
            if not os.path.exists(results_full_dir):
                os.makedirs(results_full_dir)

            for i in range(2):
                print(f'Writing data for frequency {i+1}/{2}')
                m = open(results_full_dir + f'/Coefficients_frequency_{i + 1}.out', 'w')
                for j in range(self.n_bodies*6):
                    for k in range(self.n_bodies*6):
                        f = open( f'Results_full_HAMS/hydrodynamic_coefficients_{j+1}_{k+1}.out','r')
                        data = [line.strip().split() for line in f]
                        m.write(f'{j+1} {k+1} {data[i][0]}\n')
                        f.close()
                print(f'Finished writing data for frequency {i+1}/{2}')
                m.close()