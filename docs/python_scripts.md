# Python Scripts for Pre and Post Processing

HAMS-MREL includes a set of Python scripts for generating input files, postprocessing output, and computing generalized mode shapes. The scripts are organized as follows:

```
HAMS_InputOutput/
├── Input/
│   ├── Single_body/
│   │   ├── HAMSInput.py
│   │   └── HAMS_single_bodies.py
│   └── Multi-body/
│       ├── HAMSMultiInput.py
│       └── HAMS_multi_bodies.py
└── Output/
    ├── Single_body/
    │   ├── Example_single_body.py
    │   └── postprocessing_functions_hams_single_bodies.py
    └── Multi-body/
        ├── Multibodies_HAMS.py
        ├── Obtaining_RAOs_from_HAMS_PA_array.py
        └── post_processing_HAMS_multibodies.py
HAMS_to_WECSim/
├── HAMS_to_WECSim.py
└── postprocessing_functions_Mat.py
generalized_modes.py
```

Dependencies: `numpy`, `matplotlib`, `shapely`. The WEC-Sim conversion scripts additionally use `cmath`, `os`, and `datetime` (all standard library).

## Input Generation

### Single Body

**`HAMS_single_bodies.py`** contains two classes and two wrapper functions for generating the `ControlFile.in` and `Hydrostatic.in` files for a single-body simulation.

- `HAMSControlInput` — class that constructs the control file line by line. Methods write each section: water depth, wave frequencies, wave headings, solver options, pressure/elevation field points, and closing notes.
- `HAMSHydrostaticFile` — class that constructs the hydrostatic file. Writes the center of gravity, 6×6 mass matrix, external linear damping matrix, external quadratic damping matrix, hydrostatic restoring matrix, and external restoring matrix.
- `create_control_file_hams(...)` — wrapper function that instantiates `HAMSControlInput` and calls all write methods in the correct order.
- `create_hydrostatic_file(...)` — wrapper function that instantiates `HAMSHydrostaticFile` and writes the hydrostatic data.

**`HAMSInput.py`** is an example script that demonstrates how to use `HAMS_single_bodies.py`. The user sets simulation parameters as variables at the top of the script (water depth, frequency range, headings, rotation center, solver options, mass and stiffness matrices, etc.), then calls `create_control_file_hams` and `create_hydrostatic_file`. The script also creates the `Input/` and `Output/` directory structure and moves the generated files into `Input/`.

### Multi-Body

**`HAMS_multi_bodies.py`** contains the `HAMSInput` class and a `create_control_file_hams` wrapper function for generating the `ControlFile.in` for multi-body simulations. The class extends the single-body version with:

- `add_body_info()` — writes the number of bodies and the LCS (Local Coordinate System) coordinates and orientation for each body.
- Per-body rotation centers (`Reference_body_center_i`) in the solver options section.

**`HAMSMultiInput.py`** is an example script that demonstrates how to use `HAMS_multi_bodies.py`. It sets up a 10-body simulation with LCS coordinates and rotation centers per body. Note that the `Hydrostatic[_i].in` files must be created separately for multi-body cases.


## Output Postprocessing

### Single Body

**`postprocessing_functions_hams_single_bodies.py`** provides two functions for reading HAMS output files in Hams_format:

- `get_mass_and_damping_matrices_hams(added_mass_file_location, added_damping_file_location, dof)` — reads `OAMASS{i}.txt` and `ODAMPING{i}.txt` and returns the frequency list, added mass coefficients, and radiation damping coefficients for a given DOF.
- `get_external_force_hams(external_force_file_location, column)` — reads `OEXFOR{i}.txt` and returns the excitation force values for a given column (real or imaginary part).

**`Example_single_body.py`** is an example script that demonstrates postprocessing for a single point absorber in heave. It reads the HAMS output, plots added mass, radiation damping, and excitation force coefficients, solves the frequency-domain equation of motion with and without a PTO system, computes the natural frequency, and plots the RAOs and average absorbed power. The dynamic equation solved is:

```
[-ω²(m + a₃₃) + iω(b₃₃ + B_PTO) + c₃₃] · η₃₃ = Fₑ
```

### Multi-Body

**`post_processing_HAMS_multibodies.py`** provides functions for reading and reorganizing HAMS output for multi-body simulations:

- `get_mass_and_damping_matrices_hams_multibodies(...)` — reads added mass and radiation damping from HAMS output files for a given DOF. Same interface as the single-body version.
- `get_external_force_hams_multibodies(...)` — reads excitation forces from HAMS output files.
- `compile_mass_damping_matrices_multiple_bodies(n_bodies, length_lists, added_mass_final, radiation_damping_final)` — reorganizes the raw coefficient lists into a compiled format indexed by body pair. Returns lists of size `n_bodies²`, each containing the coefficient values across all frequencies.
- `compile_excitation_force_matrix_multiple_bodies(n_bodies, length_lists, exciting_force_final, n_directions)` — reorganizes excitation force lists indexed by body and wave direction. Returns a list of size `n_bodies × n_directions`.
- `hydrodynamic_coefficients_exciting_forces_all_bodies(...)` — main function that reads all HAMS output files for a DOF, optionally strips zero/infinite frequency entries, and returns the compiled added mass, radiation damping, and excitation force matrices along with the number of bodies and wave directions.
- `hydrodynamic_coefficients_exciting_forces_all_bodies_zero_inf(...)` — same as above but returns only the zero and infinite frequency coefficients.
- `obtain_specific_coefficient_matrix(added_mass_compiled, radiation_damping_compiled, n_bodies, x, y)` — retrieves added mass and radiation damping for a specific body pair (`x`, `y`).
- `obtain_specific_exciting_force_matrix(exciting_force_compiled, x, j, n_directions)` — retrieves excitation force for a specific body `x` and direction `j`.

**`Multibodies_HAMS.py`** contains the `Multibodies_HAMS` class that provides higher-level postprocessing for multi-body simulations. Key methods:

- `get_number_of_bodies()` — reads the number of bodies from `ControlFile.in`.
- `get_frequency_information()` — parses all simulation parameters from `ControlFile.in`, including frequencies, headings, body positions, water depth, and number of panels per body. Converts between frequency types (wave number, frequency, period) including a shallow-water dispersion solver.
- `hams_coefficient_conversion(body_number_1, body_number_2, dof_1, dof_2)` — converts body/DOF pairs to global DOF indices in the combined `6×NBODY` system.
- `create_hams_compiled_results(hams_output_location, run)` — reads all HAMS output files, compiles hydrodynamic coefficients and excitation forces across all DOF combinations, and writes them to individual text files in `Results_full_HAMS/`.
- `create_hams_compiled_results_zero_inf(hams_output_location, run)` — same as above but for zero and infinite frequency only.
- `compile_hams_results_per_frequency(run)` — reads the compiled results from `Results_full_HAMS/` and reorganizes them into per-frequency output files in `Results_frequency_HAMS/`, with one file per frequency containing all DOF combinations.
- `compile_hams_results_per_frequency_zero_inf(run)` — same as above for zero/infinite frequencies.

**`Obtaining_RAOs_from_HAMS_PA_array.py`** is an example script that demonstrates RAO computation for an array of point absorbers. It uses the `Multibodies_HAMS` class to compile the per-frequency coefficient and excitation matrices, then solves the coupled equation of motion:

```
[-ω²(M + A) + iωB + K] · ξ = Fₑ
```

where `M`, `A`, `B`, `K` are `6N × 6N` block-diagonal matrices (N = number of bodies) and `Fₑ` is the `6N × 1` excitation force vector. The resulting RAOs are written to a text file in WAMIT format.


### Multi-Body (Matrix-Based Variant)

**`postprocessing_functions_Mat.py`** provides an alternative set of multi-body postprocessing functions with a simplified interface tailored for matrix-based workflows (e.g. WEC-Sim conversion). The functions mirror those in `post_processing_HAMS_multibodies.py` but omit wave-direction indexing, assuming a single wave heading. Key functions:

- `get_mass_and_damping_matrices_hams_multibodies(added_mass_file_location, added_damping_file_location, dof)` — reads HAMS output files and returns the frequency list, added mass coefficients, and radiation damping coefficients for a given DOF.
- `get_external_force_hams_multibodies(external_force_file_location, column)` — reads excitation force values for a given column (real or imaginary part).
- `compile_mass_damping_matrices_multiple_bodies(n_bodies, length_lists, added_mass_final, radiation_damping_final)` — reorganizes raw coefficient lists into compiled format indexed by body pair. Returns lists of size `n_bodies²`.
- `compile_excitation_force_matrix_multiple_bodies(n_bodies, length_lists, exciting_force_final)` — reorganizes excitation force lists indexed by body. Returns a list of size `n_bodies`.
- `hydrodynamic_coefficients_exciting_forces_all_bodies(...)` — main function that reads all HAMS output files for a DOF, optionally strips zero/infinite frequency entries, constructs complex excitation forces from real and imaginary parts, and returns the compiled matrices.
- `hydrodynamic_coefficients_exciting_forces_all_bodies_zero_inf(...)` — same as above but returns only the zero and infinite frequency coefficients.
- `obtain_specific_coefficient_matrix(added_mass_compiled, radiation_damping_compiled, n_bodies, x, y)` — retrieves added mass and radiation damping for a specific body pair (`x`, `y`).
- `obtain_specific_exciting_force_matrix(exciting_force_compiled, x)` — retrieves excitation force for a specific body `x`.


## Generalized Modes

**`generalized_modes.py`** computes generalized mode shape functions for use with the generalized modes feature in HAMS-MREL. The script produces `gen_mod_{body}_{dof}.txt` files that can be placed in the input directory to override standard rigid-body normals. See the [Input Files](input_files.md#gen_mod-optional) documentation for details on how these files are read by the Fortran code.

The script contains the following functions:

- `read_pnl_mesh(filepath)` — reads a HAMS `.pnl` hull mesh file and returns node coordinates and panel connectivity.
- `extract_body_number(mesh_filepath)` — extracts the body number from the mesh filename (e.g., `HullMesh_2.pnl` → 2).
- `read_lcs_from_controlfile(control_filepath, body_number)` — reads the LCS translation offsets and z-axis rotation angle for a given body from `ControlFile.in`.
- `local_to_global(coords_local, x0, y0, z0, rot_deg)` — transforms coordinates from the local body frame to the global frame by applying z-axis rotation followed by translation.
- `compute_centroids(nodes, panels)` — computes the centroid of each panel as the arithmetic mean of its vertex coordinates.
- `compute_generalized_normal(centroids_global, nx, ny, nz)` — evaluates the generalized mode shape function at each panel centroid as `n_j = nx * X + ny * Y + nz * Z`.

The user sets the mesh file, control file, DOF, and shape function components (`nx`, `ny`, `nz`) as variables in `main()`, then runs the script to produce the output file.


## HAMS to WEC-Sim Conversion

**`HAMS_to_WECSim.py`** converts HAMS-MREL multi-body results into a WAMIT-format `.out` file that can be used as input for WEC-Sim. The script imports postprocessing functions from `postprocessing_functions_Mat.py` and contains the `HAMS_to_WAMIT` class with the following methods:

- `__init__(filename, hams_result_folder, hams_input_folder, rho)` — initializes the converter with output filename, HAMS directory paths, and water density (default 1000 kg/m³).
- `create_initial_wamit_output()` — writes the WAMIT-format header block to the output file, including MREL attribution.
- `get_number_of_bodies()` — reads the number of bodies from `ControlFile.in`.
- `get_frequency_information(g)` — parses all simulation parameters from `ControlFile.in` including frequencies, headings, body positions (LCS), rotation centers, water depth, and panel counts. Converts between frequency types (wave number, frequency, period) for all four HAMS input frequency types. Stores period, wavenumber, and angular frequency lists for WAMIT output formatting.
- `obtain_wavenumber_shallow_water(period, depth, g)` — iterative shallow-water dispersion solver that computes the wavenumber from wave period and water depth.
- `hams_to_wamit_coefficient_conversion(body_number_1, body_number_2, dof_1, dof_2)` — maps body/DOF pairs to global indices in the combined `6N` system used by WAMIT.
- `hams_to_wamit_coefficient_conversion_excitation(body_number_1, dof_1)` — maps body/DOF to a global excitation force index.
- `create_hams_compiled_results_for_wamit(hams_output_location, run)` — reads all HAMS output files, compiles hydrodynamic coefficients (added mass, radiation damping) and excitation forces across all DOF combinations, and writes them to individual text files in `Results_full_HAMS/`. Excitation forces are stored as magnitude and phase.
- `create_hams_compiled_results_for_wamit_zero_inf(hams_output_location, run)` — same as above but for zero and infinite frequency added mass only, writing to `Results_full_HAMS_inf_zero/`.
- `compile_hams_results_per_frequency_for_wamit(run)` — reads the compiled results from `Results_full_HAMS/` and reorganizes them into per-frequency files in `Results_frequency_HAMS/`, with one coefficient file and one excitation file per frequency.
- `compile_hams_results_per_frequency_for_wamit_zero_inf(run)` — same as above for zero/infinite frequencies, writing to `Results_frequency_HAMS_inf_zero/`.
- `add_input_file_description()` — writes the WAMIT input file description block listing hull mesh files.
- `add_frequency_information()` — writes the frequency table, gravity, water depth, solver parameters, and symmetry information in WAMIT format.
- `add_hydrostatic_data(volumes, hydrostatic, center_of_gravity, center_of_bouyancy, g)` — writes per-body hydrostatic data including displaced volumes, center of buoyancy, hydrostatic restoring coefficients (non-dimensionalized by `ρg`), center of gravity, and radii of gyration.
- `add_added_mass_zero_inf()` — writes the zero and infinite frequency added mass coefficient blocks, non-dimensionalized by `ρ`.
- `add_hydrodynamic_coefficients_per_frequency(omega_number)` — writes the added mass (non-dimensionalized by `ρ`) and radiation damping (non-dimensionalized by `ρω`) for a given frequency.
- `add_excitation_forces_per_frequency(omega_number, g)` — writes the diffraction exciting force magnitude (non-dimensionalized by `ρg`) and phase for a given frequency and wave heading.
- `add_hydrodynamic_coefficients_excitation_forces()` — loops over all frequencies and writes the combined coefficient and excitation force blocks.
- `add_initial_output_info()` — writes the WAMIT FORCE run header.

The script includes a two-stage execution workflow controlled by flags at the bottom:

1. **`create_data_from_HAMS`** — reads raw HAMS output, compiles all DOF combinations, and produces intermediate per-frequency result files.
2. **`create_final_WECSIM_files`** — reads the intermediate files and assembles the final WAMIT `.out` file with user-specified hydrostatic data, volumes, and centers of gravity/buoyancy.
