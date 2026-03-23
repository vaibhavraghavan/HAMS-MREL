# Linux Installation Instructions

This document describes how to set up a development environment that allows you to compile and run HAMS-MREL on Linux (Ubuntu) systems, including the DelftBlue High Performance Computing Cluster.

## Step 1: Install Dependencies

HAMS-MREL requires the following software:

### Intel Fortran Compiler and oneMK Library

Download and install the [Intel Fortran Essentials](https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2025-0/intel-fortran-essentials.html) bundle from the [Intel website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?packages=fortran-essentials&fortran-essentials-os=linux&fortran-essentials-lin=offline).

There are two types of types of installers: offline and online. Both can be installed with the terminal or via a GUI. For detailed installation instructions see the official documentation [here](https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2025-0/online-offline-installer-001.html).

**Important**
After installation you must run the `setvars.sh` script to add the compilers and libraries to your system's PATH:
```bash
# Open installation directory (typically /home/<user>/intel/oneapi)
cd <installation-directory> 

# Run the script
source setvars.sh

# Test the installation
ifx --version
```

### Fortran Package Manager (optional)

HAMS-MREL can optionally be built using the [Fortran Package Manager](https://fpm.fortran-lang.org). Using fpm provides a convenient and reproducible way to build the application, handle dependencies, and manage compiler options.

## Step 2: Download HAMS-MREL

You may obtain the HAMS-MREL source code using one of the following methods.

- **Clone the repository using Git**: if Git is installed, open a Command Prompt or PowerShell and run:
```
git clone https://github.com/vaibhavraghavan/HAMS-MREL.git
```

- **Download a ZIP archive**: Open the [GitHub repository](https://github.com/vaibhavraghavan/HAMS-MREL) in a web browser. Click "Code -> Download ZIP". Extract the contents of the ZIP file to a local directory of your choice.


## Step 3: Build and Execute

Run the following commands from the terminal:

```bash
# Navigate to the root directory of the HAMS-MREL repository
cd <path-to-HAMS-MREL>

# Ensure the Intel oneAPI environment is loaded
source ~/intel/oneapi/setvars.h     # path may be different on your system 

# Open the src directory
cd src

# Compile with the IFX Intel compiler
ifx -qmkl -o HAMS_MREL WavDynMods.f90 InputFiles.f90 WavDynSubs.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 AssbMatx_irr.f90 AssbMatx.f90 AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 SolveMotion.f90 SolveMotionMulti.f90 PrintOutput.f90 HAMS_Prog.f90 -qopenmp

# Run the appplication
./hamsmrel <path-to-input-files> <path-to-output-directory>
```

This will create an executable called `hamsmrel` (specified by the `-o` option) inside the HAMS-MREL/src directory. The `-qmkl` and `-qopenmp` options indicate that the MKL and OpenMP libraries must be linked.

> [!IMPORTANT]
> By default, HAMS-MREL looks for input files in the executable’s directory. You can also specify input and output paths as command-line arguments. Output files will be written to the specified directory or the executable’s directory if none is provided. See [input_files.md](input_files.md) and [output_files.md](output_files.md) for a full description of the input and output files.

### Build and Execute Using the Fortran Package Manager

Alternatively, you may build, install and run HAMS-MREL using the Fortran Package Manager. This option makes use of the [fpm.toml](../fpm.toml) file which defines source files, dependencies, etc.

To build and run using FPM:

```bash
# Navigate to the root directory of the HAMS-MREL repository
cd <path-to-HAMS-MREL>

# Ensure the Intel oneAPI environment is loaded
source ~/intel/oneapi/setvars.h     # path may be different on your system

# Set the FPM compiler and flags
export FPM_FC=IFX                   # compiler
export FPM_FFLAGS="-qmkl -qopenmp"  # flags

# Compile
fpm build

# Run the application via fpm
fpm run -- <path-to-input-files> <path-to-output-directory>

# Or copy the executable to a specified directory and run manually
# fpm install --prefix <path-to-install-dir>
# cd <path-to-install-dir>/bin
# ./hamsmrel <path-to-input-files> <path-to-output-directory>
```

## Building in a High Performance Computer Cluster

This section provides instructions for building HAMS-MREL on the [DelftBlue HPC cluster](https://doc.dhpc.tudelft.nl/delftblue/). The steps are general and can be adapted to other HPC environments that use Slurm and provide Intel compilers with the oneAPI MKL library.

### Step 1: Compile HAMS-MREL

```bash
# Load required modules
module load 2025
module load intel/oneapi-all

# Clone the repository
git clone https://github.com/vaibhavraghavan/HAMS-MREL.git

# Navigate to the src directory
cd HAMS-MREL-DCC/src/

# Compile with the IFX compiler
ifx -qmkl -o HAMS_MREL WavDynMods.f90 InputFiles.f90 WavDynSubs.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 AssbMatx_irr.f90 AssbMatx.f90 AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 SolveMotion.f90 SolveMotionMulti.f90 PrintOutput.f90 HAMS_Prog.f90 -qopenmp
```

This produces the executable `hamsmrel` in the `src` directory. The `-qmkl` and `-qopenmp` flags link the MKL library and enable OpenMP parallelization.

### Step 2: Submit a Slurm job

Create a Slurm job script, for example `slurmjob.sh`:

```bash
#!/bin/sh
#SBATCH --job-name=hamsmulti
#SBATCH --partition=compute
#SBATCH --time=00:10:00         # job duration in hh:mm:ss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32      # number of OpenMP threads
#SBATCH --mem-per-cpu=1G

# Load Intel compilers
module load intel/oneapi-all

# Set the number of OpenMP threads 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the application
srun <path-to-hams-repository>/src/hamsmrel <input-files-path> <output-directory-path>
```

Submit the job with:
```bash
sbatch slurmjob.sh
```

This example runs HAMS-MREL on 1 node, 1 task, and 32 threads. Adjust `--cpus-per-task` depending on your simulation size and cluster resource allocation. This number must match the number of threads specified in the input file `ControlFile.in`.

Note that DelftBlue imposes strict memory quotas on the `/home` directory (approx. 30 GBs). Run jobs in `/scratch` to avoid exceeding memory limits.


