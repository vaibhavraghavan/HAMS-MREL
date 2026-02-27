"""
    Created on Fri Oct 4 11:57:48 2024

    This script presents the postprocessing for an array of point absorbers for calculation of the RAOs. This is added
    as an example

    Original author: Vaibhav Raghavan
"""
import cmath
import math
import os

import numpy as np
from Multibodies_HAMS import Multibodies_HAMS

# Input folder for the multi-body simulation in HAMS-MREL
hams_input_folder = r'Results_multibodies_PA\Input'
output_file = Multibodies_HAMS(filename='hams_output', hams_input_folder=hams_input_folder,
                               hams_result_folder='HAMS_case')

# Switch for creating coefficients per frequency
create_coefficients_per_frequency = True
# Switch to calculate the RAOs. Note this is specific to this example.
obtain_RAO = False

output_file.get_number_of_bodies()
output_file.get_frequency_information()
number_of_bodies = output_file.n_bodies

###################################### Mass and stiffness matrix#######################################################
inertia_rotational =  304537
mass_matrix_single = np.array([[0.30454E+06,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00],
                              [0.00000E+00,0.30454E+06,0.00000E+00,0.00000E+00, 0.00000E+00,0.00000E+00],
                               [0.00000E+00,0.00000E+00,0.30454E+06,0.00000E+00,0.00000E+00, 0.00000E+00],
                               [0.00000E+00, 0.00000E+00, 0.00000E+00, inertia_rotational, 0.00000E+00, 0.00000E+00],
                               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, inertia_rotational, 0.00000E+00],
                               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, inertia_rotational]])
stiffness_matrix_single = np.array([[0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00],
                              [0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00, 0.00000E+00,0.00000E+00],
                               [0.00000E+00,0.00000E+00,0.27275E+06,0.00000E+00,0.00000E+00, 0.00000E+00],
                               [0.00000E+00, 0.00000E+00, 0.00000E+00, -0.15754E+08, 0.00000E+00, 0.00000E+00],
                               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, -0.15754E+08, 0.00000E+00],
                               [0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00]])
mass_matrix_compiled = np.zeros([6*number_of_bodies, 6*number_of_bodies])
stiffness_matrix_compiled = np.zeros([6*number_of_bodies, 6*number_of_bodies])
for i in range(number_of_bodies):
    mass_matrix_compiled[6*i:6*(i+1), 6*i:6*(i+1)] = mass_matrix_single
    stiffness_matrix_compiled[6*i:6*(i+1),6*i:6*(i+1)] = stiffness_matrix_single

###################################### Added mass, Radiation damping and excitation force matrices######################
if create_coefficients_per_frequency:
    # Output folder for
    hams_output_location = r'Results_multibodies_PA\Output\Hams_format'

    # Compiling the added mass, radiation damping matrices and excitation force, and creating the results for all DOFs
    output_file.create_hams_compiled_results(hams_output_location=hams_output_location, run=True)

    # Convert the compiled files into a format can be used to easily for postprocessing
    output_file.compile_hams_results_per_frequency(run=True)

####################################### Obtaining the RAOs #############################################################
if obtain_RAO:
    # RAOs
    omega_list = output_file.omega_print[2:]
    RAO_matrix_compiled = []

    for i in range(len(omega_list)):
        f = open(os.getcwd() + fr'/Results_frequency_HAMS/Coefficients_frequency_{i+1}')
        lines = [line.strip().split() for line in f]

        # Added mass and radiation damping

        added_mass_total = [float(ele[2]) for ele in lines]
        added_mass_compiled = np.array(added_mass_total).reshape((6 * number_of_bodies, 6 * number_of_bodies))
        radiation_damping_total = [float(ele[3]) for ele in lines]
        radiation_damping_compiled = np.array(radiation_damping_total).reshape((6 * number_of_bodies, 6 * number_of_bodies))
        f.close()
        mass_mod = -omega_list[i]**2 * (mass_matrix_compiled + added_mass_compiled)
        damping_mod = 1j * omega_list[i] * radiation_damping_compiled

        # Excitation force
        f = open(os.getcwd() + fr'/Results_frequency_HAMS/Excitation_frequency_{i + 1}')
        lines = [line.strip().split() for line in f]
        excitation_compiled = np.array([float(ele[1]) + 1j * float(ele[2]) for ele in lines]) # Correction based on the RAOs from WAMIT
        f.close()

        # RAO matrix
        RAO_matrix = np.matmul(np.linalg.inv(mass_mod + damping_mod + stiffness_matrix_compiled), excitation_compiled)
        RAO_matrix_compiled.append(RAO_matrix)

    # Creating the output file in WAMIT format to run the HAMS plotting code

    text_file = open('raos_hams.txt','w')
    period_list = output_file.period_print[2:]

    for i in range(len(period_list)):
        for j in range(6 * number_of_bodies):
            text_file.write(f'  {period_list[i]:.6E}  {output_file.incidence_angle:.6E}     {j+1}  {abs(RAO_matrix_compiled[i][j]):.6E} {cmath.phase(RAO_matrix_compiled[i][j]) * 180/math.pi:.6E} {RAO_matrix_compiled[i][j].real:.6E} {RAO_matrix_compiled[i][j].imag:.6E}\n')


