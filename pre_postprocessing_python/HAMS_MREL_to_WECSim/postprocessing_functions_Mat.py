import math

def get_mass_and_damping_matrices_hams_multibodies(added_mass_file_location: str, added_damping_file_location: str, dof: int):

    # Mass matrix
    # Importing the files
    with open(added_mass_file_location) as f:
        lines = [line.strip().split() for line in f]

    # Obtaining the data
    omega = []
    mass_matrix_for_all_frequencies = []
    for i in range(len(lines)):
        omega.append(float(lines[i][1]))
        mass_matrix_for_all_frequencies.append([float(lines[i][2]), float(lines[i][3]), float(lines[i][4]), float(lines[i][5]), float(lines[i][6]),
                                                         float(lines[i][7])])
    mass_all_frequencies = [mass_matrix_for_all_frequencies[i][dof-1] for i in
                                    range(len(mass_matrix_for_all_frequencies))]

    # Damping matrix
    # Importing the files
    with open(added_damping_file_location) as f:
        lines = [line.strip().split() for line in f]

    # Obtaining the elements
    damping_matrix_for_all_frequencies = []
    for i in range(len(lines)):
        damping_matrix_for_all_frequencies.append(
            [float(lines[i][2]), float(lines[i][3]), float(lines[i][4]), float(lines[i][5]), float(lines[i][6]),
             float(lines[i][7])])
    damping_all_frequencies = [damping_matrix_for_all_frequencies[i][dof-1] for i in
                                    range(len(damping_matrix_for_all_frequencies))]

    if len(omega) != len(mass_all_frequencies) or len(omega) != len(damping_all_frequencies):
        raise ValueError('The lengths of the mass/damping list and different from the frequency lists. Please check.')

    return omega, mass_all_frequencies, damping_all_frequencies

def get_external_force_hams_multibodies(external_force_file_location: str, column: int):

    # Importing the files
    with open(external_force_file_location) as f:
        lines = [line.strip().split() for line in f]

    # Getting the data
    omega = []
    external_force_all_frequencies = []
    for i in range(len(lines)):
        omega.append(float(lines[i][1]))
        external_force_all_frequencies.append(float(lines[i][column + 1]))

    if len(omega) != len(external_force_all_frequencies):
        raise ValueError('The lengths of the external force list and different from the frequency lists. Please check.')

    return external_force_all_frequencies

def compile_mass_damping_matrices_multiple_bodies(n_bodies: int, length_lists: int, added_mass_final: list,
                                                  radiation_damping_final: list):
    added_mass_compiled = []
    radiation_damping_compiled = []


    for j in range(n_bodies**2):
        added_mass_compiled.append([added_mass_final[n_bodies**2 * k + j] for k in range(length_lists)])
        radiation_damping_compiled.append([radiation_damping_final[n_bodies ** 2 * k + j] for k in range(length_lists)])

    return added_mass_compiled, radiation_damping_compiled

def compile_excitation_force_matrix_multiple_bodies(n_bodies: int, length_lists: int, exciting_force_final: list):

    exciting_force_compiled = []

    for j in range(n_bodies):
        exciting_force_compiled.append([exciting_force_final[n_bodies * k + j] for k in range(length_lists)])

    return exciting_force_compiled

def hydrodynamic_coefficients_exciting_forces_all_bodies(added_damping_file_location: str,
                                                         added_mass_file_location: str,
                                                         external_force_file_location: str,
                                                         dof: int, zero_inf_frequency: bool):
    omega, mass_all_frequencies, damping_all_frequencies = \
        get_mass_and_damping_matrices_hams_multibodies(added_damping_file_location=added_damping_file_location,
                                                       added_mass_file_location=added_mass_file_location,
                                                       dof=dof)
    external_force_all_frequencies_1 = get_external_force_hams_multibodies(
        external_force_file_location=external_force_file_location,
        column=1)
    external_force_all_frequencies_2 = get_external_force_hams_multibodies(
        external_force_file_location=external_force_file_location,
        column=2)

    # Getting unique list of frequencies not considering the 0 and infinite frequencies
    omega_final = list(set(omega))
    n_bodies = int(math.sqrt(len(omega) / len(omega_final))) # No of bodies
    omega_final.sort()
    if zero_inf_frequency:
        omega_final = omega_final[2:] # Excluding the 0 and infinite frequencies
        # Sorting the added mass and damping matrices
        added_mass_final = mass_all_frequencies[n_bodies**2 * 2:]
        radiation_damping_final = damping_all_frequencies[n_bodies**2 * 2:]
        exciting_forces_final = [external_force_all_frequencies_1[i] + complex(0, 1) * external_force_all_frequencies_2[i]
                               for i in range(len(external_force_all_frequencies_1))][2*n_bodies:]
    else:
        # Sorting the added mass and damping matrices
        added_mass_final = mass_all_frequencies
        radiation_damping_final = damping_all_frequencies
        exciting_forces_final = [external_force_all_frequencies_1[i] + complex(0, 1) * external_force_all_frequencies_2[
            i] for i in range(len(external_force_all_frequencies_1))]
    list_size = len(omega_final)

    # Creating the mass, damping and excitation force compiled matrices
    added_mass_compiled, radiation_damping_compiled = (
        compile_mass_damping_matrices_multiple_bodies(n_bodies=n_bodies, length_lists=list_size,
                                                      added_mass_final=added_mass_final,
                                                      radiation_damping_final=radiation_damping_final))

    exciting_force_compiled = compile_excitation_force_matrix_multiple_bodies(n_bodies=n_bodies, length_lists=list_size,
                                                                              exciting_force_final=exciting_forces_final)

    return n_bodies, omega_final, added_mass_compiled, radiation_damping_compiled, exciting_force_compiled

def hydrodynamic_coefficients_exciting_forces_all_bodies_zero_inf(added_damping_file_location: str,
                                                         added_mass_file_location: str,
                                                         external_force_file_location: str,
                                                         dof: int):
    omega, mass_all_frequencies, damping_all_frequencies = \
        get_mass_and_damping_matrices_hams_multibodies(added_damping_file_location=added_damping_file_location,
                                                       added_mass_file_location=added_mass_file_location,
                                                       dof=dof)
    external_force_all_frequencies_1 = get_external_force_hams_multibodies(
        external_force_file_location=external_force_file_location,
        column=1)
    external_force_all_frequencies_2 = get_external_force_hams_multibodies(
        external_force_file_location=external_force_file_location,
        column=2)

    # Getting unique list of frequencies not considering the 0 and infinite frequencies
    omega_final = list(set(omega))
    n_bodies = int(math.sqrt(len(omega) / len(omega_final))) # No of bodies
    omega_final.sort()
    omega_final = omega_final[:2] # Only zero and infinite frequencies
    # Sorting the added mass and damping matrices
    added_mass_final = mass_all_frequencies[:n_bodies**2 * 2]
    radiation_damping_final = damping_all_frequencies[:n_bodies**2 * 2]
    exciting_forces_final = [external_force_all_frequencies_1[i] + complex(0, 1) * external_force_all_frequencies_2[i]
                           for i in range(len(external_force_all_frequencies_1))][:2*n_bodies]
    list_size = len(omega_final)

    # Creating the mass, damping and excitation force compiled matrices
    added_mass_compiled, radiation_damping_compiled = (
        compile_mass_damping_matrices_multiple_bodies(n_bodies=n_bodies, length_lists=list_size,
                                                      added_mass_final=added_mass_final,
                                                      radiation_damping_final=radiation_damping_final))

    exciting_force_compiled = compile_excitation_force_matrix_multiple_bodies(n_bodies=n_bodies, length_lists=list_size,
                                                                              exciting_force_final=exciting_forces_final)

    return n_bodies, omega_final, added_mass_compiled, radiation_damping_compiled, exciting_force_compiled

def obtain_specific_coefficient_matrix(added_mass_compiled: list, radiation_damping_compiled: list, n_bodies: int,
                                       x: int, y: int):
    index = n_bodies * (x - 1) + (y - 1)
    return added_mass_compiled[index], radiation_damping_compiled[index]


def obtain_specific_exciting_force_matrix(exciting_force_compiled: list, x: int):
    return exciting_force_compiled[x - 1]