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
generalized_modes.py
```

Dependencies: `numpy`, `matplotlib`, `shapely`.

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
