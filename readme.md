# HAMS-MREL

HAMS-MREL (Hydrodynamic Analysis of Marine Structures – Marine Renewable Energies Lab) is an open-source numerical tool for studying wave–structure interactions involving multiple floating bodies. It uses a Boundary Integral Equation Method to compute hydrodynamic coefficients, wave excitation forces, and related quantities for diffraction and radiation problems. The software includes OpenMP-based parallelization to improve computational efficiency and has been validated against WAMIT for a range of multi-body test cases.

HAMS-MREL is developed as an extension of the original [HAMS solver](https://github.com/YingyiLiu/HAMS). While HAMS focuses on single-body hydrodynamic analysis, HAMS-MREL expands the framework to handle interactions among multiple bodies.


## Installation and Usage

Installation instructions are provided for:
- [Windows](docs/installation_windows.md)
- [Linux (Ubuntu)](docs/installation_linux.md)
- [High Performance Computing Clusters](docs/installation_linux.md)

HAMS-MREL requires input files to be placed in an input directory before execution. The program can be run in two ways, depending on whether the input and output directory locations are specified explicitly.

### 1. Default mode (Input/Output directories next to the executable)

In this mode, the program expects two directories named **Input** and **Output** to be located in the same folder as the executable.

**Windows (Command Prompt or Power Shell)**
```
hams-mrel.exe
```

**Linux (Terminal)**
```
./hams-mrel
```

### 2. Specifying input and output directory paths

The user may also provide paths to the input and output directories. These directories may be located anywhere on the system. In this mode, the program expects exactly two arguments: the input directory and the output directory. Providing only one argument is not allowed and will result in an error.

**Windows**
```
hams-mrel.exe <path-to-input-dir> <path-to-output-dir>
```

**Linux**
```
./hams-mrel <path-to-input-dir> <path-to-output-dir>
```

## Input and Output Files

### Input

HAMS-MREL requires a set of input files to run. These files contain the general setting for the simulation (`ControlFile.in`), hydrostatic data (`Hydrostatic[_i].in`) and meshes for the submerged bodies (`HullMesh[_i].pnl`) and inner water plane for each body (`WaterPlaneMesh[_i].pnl`). The exact number of files depends on the number of bodies and on whether irregular frequencies are analyzed.  

**A full description of each inout file is provided in [docs/input_files.md](docs/input_files.md)**.


### Output

HAMS-MREL produces output files in multiple established formats to support post-processing and compatibility with external tools. Every simulation outputs data in the HAMS and WAMIT formats (plus Hydrostar for single-body cases), covering added mass, damping, excitation forces, motions, pressure/elevation, and hydrostatics, along with an _ErrorCheck.txt_ file summarizing key system parameters. Output formats are fixed and always generated.

**A full description of each file is provided in [docs/output_files.md](docs/output_files.md)**.

## Developer Documentation

Additional developer-focused documentation is available in [`docs/developer_documentation.md`](docs/developer_documentation.md). This includes details on the test suite, the GitHub Actions continuous integration pipeline and results from memory profiling analysis.

## Citation

If you use HAMS‑MREL in your research, please cite:

> Raghavan, V., Loukogeorgaki, E., Mantadakis, N., Metrikine, A. V., & Lavidas, G. (2024). _HAMS‑MREL, a new open‑source multiple body solver for marine renewable energies: Model description, application and validation._ Renewable Energy, 237, 127857. [https://doi.org/10.1016/j.renene.2024.121577](https://doi.org/10.1016/j.renene.2024.121577)

> Raghavan, V., Metrikine, A.V. & Lavidas, G. Theory and validation of the new features in BIEM solver HAMS-MREL. J. Ocean Eng. Mar. Energy 12, 53–71 (2026). [https://doi.org/10.1007/s40722-025-00431-8]

## License

HAMS-MREL is available under the Apache 2.0 license. 

## Contact

HAMS‑MREL was developed by Vaibhav Raghavan. For questions, please contact the author at <email>. Users are encouraged to open a new issue in this repository to report bugs, request new features, or suggest improvements.

This work is part of the TU Delft Marine Renewable Energies Lab (MREL). Learn more about the laboratory [here](https://www.tudelft.nl/en/ceg/about-faculty/departments/hydraulic-engineering/sections/offshore-engineering/research/marine-renewable-energies-lab-mrel). 

HAMS-MREL was supported by Yasel Quintero Lares from the Digital Competence Centre, Delft University of Technology.  

