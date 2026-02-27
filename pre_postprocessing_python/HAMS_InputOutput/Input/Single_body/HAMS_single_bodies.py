class HAMSControlInput():
    """
        This class is used to create the control file for HAMS (single bodies)
    """
    def __init__(self, file_location: str, depth: float, zero_inf_limits: int, input_frequency_type: int, output_frequency_type: int,
                 number_frequencies: int, minimum_frequency: float, frequency_step: float, number_headings: int,
                 value_heading: list, reference_body_center: list, reference_body_length: float,
                 wave_diffraction_solution: int, remove_irr_freq: int, number_of_threads: int,
                 number_pressure_elevation_points: int, pressure_elevation_input: list, frequency_list: list = None,
                 minimum_heading: float = None, heading_step: float = None):

        self.file_location = file_location
        self.depth = depth
        self.input_frequency_type = input_frequency_type
        self.output_frequency_type = output_frequency_type
        self.number_frequencies = number_frequencies
        self.minimum_frequency = minimum_frequency
        self.number_headings = number_headings
        self.reference_body_center = reference_body_center
        self.reference_body_length = reference_body_length
        self.wave_diffraction_solution = wave_diffraction_solution
        self.remove_irr_freq = remove_irr_freq
        self.value_heading = value_heading
        self.number_of_threads = number_of_threads
        self.frequency_step = frequency_step
        self.zero_inf_limits = zero_inf_limits
        self.number_pressure_elevation_points = number_pressure_elevation_points
        self.pressure_elevation_input = pressure_elevation_input
        self.frequency_list = frequency_list
        self.minimum_heading = minimum_heading
        self.heading_step = heading_step
    def create_control_file(self):
        self.control_file = open("ControlFile.in", "w")
        self.control_file.write('   --------------HAMS Control file---------------\n')
        self.control_file.write('\n')
    def control_file_write_water_depth(self):
        self.control_file.write(f'   Waterdepth  {round(self.depth,1)}D0\n')
        self.control_file.write('\n')

    def control_file_write_wave_frequencies(self):
        self.control_file.write('   #Start Definition of Wave Frequencies\n')
        self.control_file.write(f'    0_inf_frequency_limits      {self.zero_inf_limits}\n')
        self.control_file.write(f'    Input_frequency_type        {self.input_frequency_type}\n')
        self.control_file.write(f'    Output_frequency_type       {self.output_frequency_type}\n')
        if self.number_frequencies < 0:
            self.control_file.write(f'    Number_of_frequencies      {self.number_frequencies}\n')
            self.control_file.write(f'    Minimum_frequency_Wmin    {round(self.minimum_frequency,1)}D0\n')
            self.control_file.write(f'    Frequency_step            {round(self.minimum_frequency,1)}D0\n')
        else:
            self.control_file.write(f'    Number_of_frequencies       {self.number_frequencies}\n')
            list_fre = ' '.join([str(float(ele)) for ele in self.frequency_list])
            self.control_file.write(f'    {list_fre}\n')
        self.control_file.write(f'   #End Definition of Wave Frequencies\n')
        self.control_file.write('\n')

    def control_file_write_wave_headings(self, number_headings):
        self.control_file.write('   #Start Definition of Wave Headings\n')
        self.control_file.write(f'    Number_of_headings          {number_headings}\n')
        print_headings = []
        if self.number_headings > 0: #TODO Add the possibility to provide negative number of wave headings
            if self.number_headings != len(self.value_heading):
                raise ValueError('The number of headings do not match the length of the list of headings. Please check.')
            for i in range(self.number_headings):
                print_headings.append(f'{round(self.value_heading[i],1)}' + ' ')
            separator = ''
            self.control_file.write(f'    {separator.join(print_headings)}\n')
        else:
            self.control_file.write(f'    Minimum_heading           {round(self.minimum_heading,1)}D0\n')
            self.control_file.write(f'    Heading_step             {round(self.heading_step,1)}D0\n')
        self.control_file.write(f'   #End Definition of Wave Headings\n')
        self.control_file.write('\n')

    def control_file_write_additional_terms(self):
        # This method writes the reference body center, reference body length, wave diffraction solution type,
        # if irregular frequencies should be nullified, and number of threads for parallelization

        self.control_file.write(f'    Reference_body_center       {self.reference_body_center[0]:.3f}       '
                                f'{self.reference_body_center[1]:.3f}       '
                                f'{self.reference_body_center[2]:.3f}\n')
        self.control_file.write(f'    Reference_body_length   {self.reference_body_length}D0\n')
        self.control_file.write(f'    Wave_diffrac_solution    {self.wave_diffraction_solution}\n')
        self.control_file.write(f'    If_remove_irr_freq       {self.remove_irr_freq}\n')
        self.control_file.write(f'    Number of threads       {self.number_of_threads}\n')
        self.control_file.write('\n')

    def control_file_write_pressure_elevation_input(self,number_pressure_elevation_points):
        self.control_file.write('   #Start Definition of Pressure and/or Elevation (PE)\n')
        self.control_file.write(f'    Number_of_field_points     {self.number_pressure_elevation_points}\n')
        if self.number_pressure_elevation_points != len(self.pressure_elevation_input):
          raise ValueError('The number of points for pressure/elevation calculation do not match the length of the'
                           ' list of pressure/elevation points. Please check.')
        for i in range(number_pressure_elevation_points):
            self.control_file.write(f'    {self.pressure_elevation_input[i][0]:6f}    '
                                    f'{self.pressure_elevation_input[i][1]:6f}    '
                                    f'{self.pressure_elevation_input[i][2]:6f}    Global_coords_point_{i+1}\n')
        self.control_file.write('   #End Definition of Pressure and/or Elevation\n')
        self.control_file.write('\n')

    def control_file_write_end_hams(self):
        self.control_file.write('   ----------End HAMS Control file---------------\n')
        self.control_file.write('\n')
        self.control_file.write('   Note: \n')
        self.control_file.write('   1. A negative value for depth indicated infinite depth/deep water condition \n')
        self.control_file.write('   2. Input_frequency_type options:\n')
        self.control_file.write('      1--deepwater wave number; 2--finite-depth wave number; 3--wave frequency; '
                                '4--wave period; 5--wave length\n')
        self.control_file.write('   3. Output_frequency_type options: same as Input_frequency_type options\n')
        self.control_file.write('   4. If_remove_irr_freq=1 indicates that the irregular frequencies will be '
                                'nullified\n')

class HAMSHydrostaticFile():
      """
      This class is used to create the Hydrostatic file for HAMS (single bodies)
      """

      def __init__(self,file_location: str, center_of_gravity: list, body_mass_matrix: list,
                   external_linear_damping_matrix: list, external_quadratic_damping_matrix: list,
                   hydrostatic_restoring_matrix: list, external_restoring_matrix: list):

          self.file_location = file_location
          self.center_of_gravity = center_of_gravity
          self.body_mass_matrix = body_mass_matrix
          self.external_linear_damping_matrix = external_linear_damping_matrix
          self.external_quadratic_damping_matrix = external_quadratic_damping_matrix
          self.hydrostatic_restoring_matrix = hydrostatic_restoring_matrix
          self.external_restoring_matrix = external_restoring_matrix

      def create_hydrostatic_file(self):
          self.hydrostatic_file = open("Hydrostatic.in", "w")

      def add_hydrostatic_information(self):
          self.hydrostatic_file.write(' Center of Gravity:\n')
          self.hydrostatic_file.write(f'  {float(self.center_of_gravity[0]):.15E}0  {float(self.center_of_gravity[1]):.15E}0  {float(self.center_of_gravity[2]):.15E}0\n')

          self.hydrostatic_file.write(' Body Mass Matrix:\n')
          for i in range(len(self.body_mass_matrix)):
              self.hydrostatic_file.write(
                  f'  {   float(self.body_mass_matrix[i][0]):.5E}   {float(self.body_mass_matrix[i][1]):.5E}   {float(self.body_mass_matrix[i][2]):.5E}   {float(self.body_mass_matrix[i][3]):.5E}   {float(self.body_mass_matrix[i][4]):.5E}   {float(self.body_mass_matrix[i][5]):.5E}\n')

          self.hydrostatic_file.write(' External Linear Damping Matrix:\n')
          for i in range(len(self.external_linear_damping_matrix)):
              self.hydrostatic_file.write(
                  f'  {float(self.external_linear_damping_matrix[i][0]):.5E}   {float(self.external_linear_damping_matrix[i][1]):.5E}   {float(self.external_linear_damping_matrix[i][2]):.5E}   {float(self.external_linear_damping_matrix[i][3]):.5E}   {float(self.external_linear_damping_matrix[i][4]):.5E}   {float(self.external_linear_damping_matrix[i][5]):.5E}\n')

          self.hydrostatic_file.write(' External Quadratic Damping Matrix:\n')
          for i in range(len(self.external_quadratic_damping_matrix)):
              self.hydrostatic_file.write(
                  f'  {float(self.external_quadratic_damping_matrix[i][0]):.5E}   {float(self.external_quadratic_damping_matrix[i][1]):.5E}   {float(self.external_quadratic_damping_matrix[i][2]):.5E}   {float(self.external_quadratic_damping_matrix[i][3]):.5E}   {float(self.external_quadratic_damping_matrix[i][4]):.5E}   {float(self.external_quadratic_damping_matrix[i][5]):.5E}\n')

          self.hydrostatic_file.write(' Hydrostatic Restoring Matrix:\n')
          for i in range(len(self.hydrostatic_restoring_matrix)):
              self.hydrostatic_file.write(
                  f'  {float(self.hydrostatic_restoring_matrix[i][0]):.5E}   {float(self.hydrostatic_restoring_matrix[i][1]):.5E}   {float(self.hydrostatic_restoring_matrix[i][2]):.5E}   {float(self.hydrostatic_restoring_matrix[i][3]):.5E}   {float(self.hydrostatic_restoring_matrix[i][4]):.5E}   {float(self.hydrostatic_restoring_matrix[i][5]):.5E}\n')

          self.hydrostatic_file.write(' External Restoring Matrix:\n')
          for i in range(len(self.external_restoring_matrix)):
              self.hydrostatic_file.write(
                  f'  {float(self.external_restoring_matrix[i][0]):.5E}   {float(self.external_restoring_matrix[i][1]):.5E}   {float(self.external_restoring_matrix[i][2]):.5E}   {float(self.external_restoring_matrix[i][3]):.5E}   {float(self.external_restoring_matrix[i][4]):.5E}   {float(self.external_restoring_matrix[i][5]):.5E}\n')

def create_control_file_hams(file_location: str, depth: float, zero_inf_limits: int, input_frequency_type: int, output_frequency_type: int,
                 number_frequencies: int, minimum_frequency: float, frequency_step: float, number_headings: int
                 , reference_body_center: list, reference_body_length: float,
                 wave_diffraction_solution: int, remove_irr_freq: int, number_of_threads: int,
                 number_pressure_elevation_points: int, pressure_elevation_input: list,
                             frequency_list: list = None, value_heading: list = None, minimum_heading: float = None,
                             heading_step: float = None):

    """
    This function creates the ControlFile.in input file for the radiation/diffraction analysis for HAMS
    """

    input_file = HAMSControlInput(file_location=file_location,
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
                                  frequency_list=frequency_list,
                                  heading_step=heading_step,minimum_heading=minimum_heading)
    input_file.create_control_file()
    input_file.control_file_write_water_depth()
    input_file.control_file_write_wave_frequencies()
    input_file.control_file_write_wave_headings(number_headings=number_headings)
    input_file.control_file_write_additional_terms()
    input_file.control_file_write_pressure_elevation_input(number_pressure_elevation_points=
                                                           number_pressure_elevation_points)
    input_file.control_file_write_end_hams()

def create_hydrostatic_file(file_location: str, center_of_gravity: list, body_mass_matrix: list,
                   external_linear_damping_matrix: list, external_quadratic_damping_matrix: list,
                   hydrostatic_restoring_matrix: list, external_restoring_matrix: list):
    """
    This function creates the Hydrostatic.in input file for the radiation/diffraction analysis for HAMS
    """

    input_file = HAMSHydrostaticFile(file_location=file_location,
                                     center_of_gravity=center_of_gravity, body_mass_matrix=body_mass_matrix,
                                     external_linear_damping_matrix=external_linear_damping_matrix,
                                     external_quadratic_damping_matrix=external_quadratic_damping_matrix,
                                     hydrostatic_restoring_matrix = hydrostatic_restoring_matrix,
                                     external_restoring_matrix = external_restoring_matrix)
    input_file.create_hydrostatic_file()
    input_file.add_hydrostatic_information()