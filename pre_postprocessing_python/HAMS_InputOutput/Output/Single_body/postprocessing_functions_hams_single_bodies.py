import numpy as np

def get_mass_and_damping_matrices_hams(added_mass_file_location: str, added_damping_file_location: str, dof: int):
    """
        The function obtains the added mass and radiation damping coefficients from the HAMS output files

        Input:
        added_mass_file_location: Location of the added mass output file
        added_damping_file_location: Location of the radiation damping output file
        dof: The degree of freedom for which the output should be retrieved

        Output:
        The frequency list (this depends on the input format in HAMS), added mass list and radiation damping list are returned
    """

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

def get_external_force_hams(external_force_file_location: str, column: int):
    """
        The function obtains the excitation force from the HAMS output files

        Input:
        external_force_file_location: Location of the excitation force output file
        dof: The degree of freedom for which the output should be retrieved

        Output:
        The excitation force list is returned
        """

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
