"""
    Created on Wed May 1 11:57:48 2024

    This script presents the postprocessing using a point absorber in Heave with and without PTO as a case study. 
    The properties of the PA can be seen in the paper RENEW_2022_final. The input and output for HAMS can be found in the
    folder 'PA_no_waterplane_results'.

    Original author: Vaibhav Raghavan
"""

# Packages required for running the script
from matplotlib import pyplot as plt
from shapely import LineString

from postprocessing_functions_hams_single_bodies import get_mass_and_damping_matrices_hams, get_external_force_hams
import numpy as np

# Parameters to be input by user

# dof_1 and dof_2 - If the user wants to retrieve added mass coefficient A_43, then dof_1 = 4 and dof_2 = 3.
# This is also for the radiation damping coefficient. For the excitation force, only dof_1 is sufficient.
dof_1 = 3
dof_2 = 3
result_folder = r'PB_no_waterplane_results\Hams_format'

rho = 1025 # Density of water
g = 9.81 # Acceleration due to gravity
S = 28.27 # Submerged area of Bouy
c_33 = rho * g * S # Restoring coefficient (Bouyancy force) in heaving. This can also be obtained from the KH.dat file
# (Hydrostatic stiffness matrix)
displacement_volume = 42.29 # Volume of water that is displaced. This is essentially based on the submerged height
# of the bouy
mass = displacement_volume * rho

# PTO system parameters
B_pto = 25000 # N/m, damping from the PTO system

# Plotting the coefficients, displacements and average power
plot_coeff = True
displacement_plot_without_pto = True
displacement_plot_with_pto = True
average_power_plot_with_pto = True

###########################################POSTPROCESSING###############################################################

# File locations
added_mass_file_location = result_folder + f'\OAMASS{dof_1}.txt'
added_damping_file_location = result_folder + f'\ODAMPING{dof_1}.txt'
external_force_file_location = result_folder + f'\OEXFOR{dof_1}.txt'

# Retrieving relevant information for postprocessing
omega, mass_all_frequencies, damping_all_frequencies =\
    get_mass_and_damping_matrices_hams(added_damping_file_location=added_damping_file_location,
                                  added_mass_file_location=added_mass_file_location,
                                  dof=dof_2)
# Wave excitation force is calculated assuming unit amplitude of wave excitation (In this case since there are two
# reference points indicated, the excitation force is calculated for (0,0,0) and (0,0,1.5)
external_force_all_frequencies_1 = get_external_force_hams(external_force_file_location=external_force_file_location,
                                                           column=1)
external_force_all_frequencies_2 = get_external_force_hams(external_force_file_location=external_force_file_location,
                                                           column=2)
exciting_forces_complex = [external_force_all_frequencies_1[i] + complex(0, 1) * external_force_all_frequencies_2[i]
                               for i in range(len(external_force_all_frequencies_1))]

# NOTE: For this calculation, the added mass at 0 and infinite frequency is also requested. Hence, in the plots,
# [2:] is utilized

if plot_coeff:
    # Plotting added mass coefficients
    plt.subplot(4, 1, 1)
    plt.plot(omega[2:], mass_all_frequencies[2:])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("$A_{33}$ (kg)")
    plt.grid()

    # Plotting added damping coefficients
    plt.subplot(4, 1, 2)
    plt.plot(omega[2:], damping_all_frequencies[2:])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("$B_{33}$ (Ns/m)")
    plt.grid()
    # plt.savefig(fname='A_33,B_33')

    # Plotting external force values
    plt.subplot(4, 1, 3)
    plt.plot(omega[2:], [-1*ele for ele in external_force_all_frequencies_1[2:]])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("$F_{e} (real part)$ (N)")
    plt.grid()

    plt.subplot(4, 1, 4)
    plt.plot(omega[2:], [-1*ele for ele in external_force_all_frequencies_2[2:]])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("$F_{e}(imaginary part)$ (N)")
    plt.grid()
    plt.show()

# Calculation of the natural frequency of the system

## Solving the dynamic equation in the frequency domain

# Dynamic equation in the frequency domain without PTO - [-omega^2*(mass+a_33) + i*omega*b_33 + c_33]eta_33 = F_e
displacement_all_frequencies = [exciting_forces_complex[i] / (-omega[i]**2 * (mass + mass_all_frequencies[i])
                                                                     + complex(0, 1) * omega[i] *
                                                                     damping_all_frequencies[i] + c_33) for i in
                                                                     range(len(omega))]

# Calculating the natural frequency of the system (wherever the following equation passes 0 is the natural frequency of the PA)
omega_natural_squared = [(omega[i]**2 - (c_33/(mass + mass_all_frequencies[i])))/50 for i in range(2,len(mass_all_frequencies))]
first_line = LineString(np.column_stack((omega[2:], [abs(ele) for ele in displacement_all_frequencies[2:]])))
second_line = LineString(np.column_stack((omega[2:], omega_natural_squared)))
natural_frequency = first_line.intersection(second_line)

if displacement_plot_without_pto:
    # Plotting the real part of the displacements
    plt.subplot(3, 1, 1)
    plt.plot(omega[2:], [ele.real for ele in displacement_all_frequencies][2:])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("Real($\zeta_{33}$) (m/m)")
    plt.title('Displacements of the WEC at the center of gravity')
    plt.grid()

    plt.subplot(3, 1, 2)
    plt.plot(omega[2:], [ele.imag for ele in displacement_all_frequencies][2:])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("Imag($\zeta_{33}$) (m/m)")
    plt.grid()

    plt.subplot(3, 1, 3)
    plt.plot(omega[2:], [abs(ele) for ele in displacement_all_frequencies][2:])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("Abs($\zeta_{33}$) (m/m)")
    plt.grid()
    plt.show()
    plt.clf()

    plt.plot(omega[2:], [abs(ele) for ele in displacement_all_frequencies[2:]])
    plt.plot(omega[2:], omega_natural_squared)
    plt.plot(*natural_frequency.xy, 'o')
    plt.xlim([0, 4])
    plt.xlabel("$\omega$ (rad/s)")
    plt.legend(['Abs($z_{33}/\zeta_{33}$) (m/m)', '$x50 (\omega^{2} -(c_{33}/(m + a_{33})))  (rad^{2}/s^{2})$'])
    plt.grid()
    plt.show()
    plt.clf()

    print(f'Natural frequency: {natural_frequency.x} rad/s') # Printing the natural frequency of the PA

# Dynamic equation in the frequency domain with linear PTO -
# [-omega^2*(mass+a_33) + i*omega*(b_33 + B_PTO) + c_33]eta_33 = F_e
displacement_all_frequencies = [exciting_forces_complex[i] / (-omega[i]**2 * (mass + mass_all_frequencies[i])
                                                                     + complex(0, 1) * omega[i] *
                                                                     (damping_all_frequencies[i] + B_pto) + c_33) for i
                                                                     in range(len(omega))]

if displacement_plot_with_pto:
    # Plotting the real part of the displacements
    plt.subplot(3, 1, 1)
    plt.plot(omega[2:], [ele.real for ele in displacement_all_frequencies][2:])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("Real($\zeta_{33}$) (m/m)")
    plt.title('Displacements of the WEC at the center of gravity')
    plt.grid()

    plt.subplot(3, 1, 2)
    plt.plot(omega[2:], [ele.imag for ele in displacement_all_frequencies][2:])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("Imag($\zeta_{33}$) (m/m)")
    plt.grid()

    plt.subplot(3, 1, 3)
    plt.plot(omega[2:], [abs(ele) for ele in displacement_all_frequencies][2:])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("Abs($\zeta_{33}$) (m/m)")
    plt.grid()
    plt.show()

# Average power calculated
average_power_all_frequencies = [(0.5 * omega[i]**2 * B_pto * abs(displacement_all_frequencies[i])**2)/1000
                                 for i in range(len(omega))]
if average_power_plot_with_pto:
    plt.plot(omega[2:], average_power_all_frequencies[2:])
    plt.xlabel("$\omega$ (rad/s)")
    plt.ylabel("$P_{a}$ (kW)")
    plt.title('Average power generated by the WEC with')
    plt.grid()
    plt.show()













