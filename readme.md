# HAMS-MREL

HAMS-MREL (Hydrodynamic Analysis of Marine Structures – Marine Renewable Energies Lab) is an open-source numerical tool for studying wave–structure interactions involving multiple floating bodies. It uses a Boundary Integral Equation Method to compute hydrodynamic coefficients, wave excitation forces, and related quantities for diffraction and radiation problems. The software includes OpenMP-based parallelization to improve computational efficiency and has been validated against WAMIT for a range of multi-body test cases.

HAMS-MREL is developed as an extension of the original [HAMS solver](https://github.com/YingyiLiu/HAMS). While HAMS focuses on single-body hydrodynamic analysis, HAMS-MREL expands the framework to handle interactions among multiple bodies.


## Usage

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

> [!IMPORTANT]
Detailed descriptions of each input file are provided in the project Wiki.

### Input

HAMS-MREL requires a set of input files to run. The exact files depend on whether you are analyzing a single body or multiple bodies:
- `HullMesh[_i].pnl`
- `Hydrostatic[_i].in`
- `WaterPlaneMesh[_i].pnl` (only if analyzing irregular frequencies)

For multi-body cases, `_i` denotes the body number (`i = 1 … numbodies`). For a single body, no index is used (`HullMesh.pnl`, `Hydrostatic.in` and `WaterPlaneMesh.pnl`).

The program will check that all required files exist before running. If any file is missing, execution will stop.

### Output

HAMS-MREL produces output files in multiple established formats to support post-processing and compatibility with external tools. Every simulation outputs data in the HAMS and WAMIT formats (plus Hydrostar for single-body cases), covering added mass, damping, excitation forces, motions, pressure/elevation, and hydrostatics, along with an _ErrorCheck.txt_ file summarizing key system parameters. Output formats are fixed and always generated. **A full description of each file is provided in [docs/output_files.md](docs/output_files.md)**.

## Citation

If you use HAMS‑MREL in your research, please cite:

> Raghavan, V., Loukogeorgaki, E., Mantadakis, N., Metrikine, A. V., & Lavidas, G. (2024). _HAMS‑MREL, a new open‑source multiple body solver for marine renewable energies: Model description, application and validation._ Renewable Energy, 237, 127857. [https://doi.org/10.1016/j.renene.2024.121577](https://doi.org/10.1016/j.renene.2024.121577)

> Yingyi Liu (2019). _HAMS: A Frequency-Domain Preprocessor for Wave-Structure Interactions—Theory, Development, and Application._ Journal of Marine Science and Engineering, 7: 81. [https://doi.org/10.3390/jmse7030081](https://doi.org/10.3390/jmse7030081)

## License

HAMS-MREL is available under the <> license. 

## Contact

HAMS‑MREL was developed by Vaibhav Raghavan. For questions, please contact the author at <email>. Users are encouraged to open a new issue in this repository to report bugs, request new features, or suggest improvements.

Tis work is part of the TU Delft Marine Renewable Energies Lab (MREL). Learn more about the laboratory [here](https://www.tudelft.nl/en/ceg/about-faculty/departments/hydraulic-engineering/sections/offshore-engineering/research/marine-renewable-energies-lab-mrel).