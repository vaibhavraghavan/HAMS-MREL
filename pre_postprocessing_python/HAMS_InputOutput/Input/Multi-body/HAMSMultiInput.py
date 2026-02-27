# Packages to be imported
import os
import numpy as np
from HAMS_multi_bodies import create_control_file_hams

##############################################CREATING CONTROL FILE#####################################################

# Parameters to be filled in
file_location = os.getcwd() # File location of analysis
depth = 10.0                # Water depth (+ve value indicates finite depth while -ve value indicates infinite depth)
zero_inf_limits = 1         # If zero and infinite frequencies should be evaluated - 1: Evaluate, 0: Do not evaluate
input_frequency_type = 1    # Input frequency type - 1:deepwater wave number; 2:finite-depth wave number; 3:wave frequency; 4:wave period; 5:wave length
output_frequency_type = 4   # Output frequency type - Same as input frequency type
number_frequencies = -30    # Number of frequencies to be evaluated - If negative value is input, a minimum frequency and frequency step should be provided; If positive, explicitly provide the frequencies to be evaluated
minimum_frequency = 0.1     # In case of negative number_frequencies, this is the minimum frequency from which the diffraction radiation calculation will start
frequency_step = 0.1        # In case of negative number_frequencies, this is the frequency step
number_headings = 1         # Number of headings for the input wave
value_heading = [0.0]       # Direction in degrees of the input wave
reference_body_length = 1.0 # Reference body length, default is 1
wave_diffraction_solution = 2 # Wave diffraction solution type, default is 2 which will evaluate the diffraction+incident wave potential
remove_irr_freq = 0 # Removal of irregular frequencies; If 1, these will be removed
number_of_threads = 20 # Number of threads for OpenMP parallelization
number_of_bodies = 10 # Number of bodies for the multi-body simulation
local_body_coordinates = [[1, 1, -2, 0], [2, 1, -2, 0], [3, 1, -2, 0], [4, 1, -2, 0], [5, 1, -2, 0], [6, 1, -2, 0], [7, 1, -2, 0], [8, 1, -2, 0], [9, 1, -2, 0],
                          [10,1,-2,0]] # Parameters of the Local Coordinate System (LCS) for each body; The first three elements in each list represent the x,y,z coordinate of the LCS while the fourth element represent the orientation of the body w.r.t x-axis
reference_body_center = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
                        # Center of rotation for each body. This is the point about which the rotation degrees of freedom are evaluated.

# # Pressure Elevation Input (Currently unavailable, will be added in the next release)
number_pressure_elevation_points = None
pressure_elevation_input = None

# number_pressure_elevation_points = 3
# pressure_elevation_input = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]]

# pressure_elevation_input = []
# x_grid,y_grid = np.linspace(-100.0, 100.0, 200), np.linspace(-100.0, 100.0, 200)
# for i in range(len(x_grid)):
#     for j in range(len(y_grid)):
#         pressure_elevation_input.append([x_grid[i], y_grid[j], 0.0])
# number_pressure_elevation_points = len(pressure_elevation_input)

create_control_file_hams(file_location=file_location,
                         depth=depth,
                         zero_inf_limits=zero_inf_limits,
                         input_frequency_type=input_frequency_type,
                         output_frequency_type=output_frequency_type,
                         number_frequencies=number_frequencies,
                         minimum_frequency=minimum_frequency,
                         frequency_step=frequency_step,
                         number_headings=number_headings,
                         value_heading=value_heading,
                         reference_body_center=reference_body_center,
                         reference_body_length=reference_body_length,
                         wave_diffraction_solution=wave_diffraction_solution,
                         remove_irr_freq=remove_irr_freq,
                         number_of_threads=number_of_threads,
                         number_pressure_elevation_points=number_pressure_elevation_points,
                         pressure_elevation_input=pressure_elevation_input,
                         number_of_bodies=number_of_bodies,
                         coordinates_body_center=local_body_coordinates)

##############################################CREATING THE REMAINING FOLDERS AND FILES##################################
if not os.path.isdir('Input'):
    os.makedirs('Input')
if not os.path.isdir('Output'):
    os.makedirs('Output')
if os.path.isfile('ControlFile.in'):
    os.replace('ControlFile.in','Input/ControlFile.in')
else:
    print('The ControlFile is not found. Please check.')
if os.path.isfile('Hydrostatic.in'):
    os.replace('Hydrostatic.in','Input/Hydrostatic.in')
else:
    print('The Hydrostatic file is not found. Please check.')

if os.path.isdir('Output'):
    if not os.path.isdir('Output/Hams_format'):
        os.makedirs('Output/Hams_format')
    if not os.path.isdir('Output/Hydrostar_format'):
        os.makedirs('Output/Hydrostar_format')
    if not os.path.isdir('Output/Wamit_format'):
        os.makedirs('Output/Wamit_format')
    if not os.path.isfile('Output/ErrorCheck.txt'):
        f=open("Output/ErrorCheck.txt", "w")
        f.close()


