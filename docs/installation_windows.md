# Windows Installation Instructions

This document describes how to set up a development environment that allows you to compile and run HAMS-MREL on Windows.

## Step 1: Install Dependencies

HAMS-MREL requires the following software:

### Intel Fortran Compiler and oneMK Library

Download and install the [Intel Fortran Essentials](https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-windows/2025-0/intel-fortran-essentials.html) bundle (includes the oneAPI MKL library) bundle from the [Intel website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?packages=fortran-essentials&fortran-essentials-os=windows&fortran-essentials-win=offline).

There are two types of types of installers: offline and online. Both can be installed with the command line (Command Prompt or Windows PowerShell) or via a GUI. If installing with the command line, make sure you activate the integration with Visual Studio by using the `NEED_VS<yyy>_INTEGRATION` property. For example, if Working With Visual Studio 2022:

```shell
./intel-fortran-essentials-2025.1.0.555_offline.exe -a --silent --eula accept -p=NEED_VS2022_INTEGRATION=1
```

**Important**
After installation you must run the `setvars.bat` script to add the libraries to your system PATH. This ensures the Intel compiler and libraries are discoverable by Visual Studio. The script is typically located in `C:\Program Files (x86)\Intel\oneAPI\setvars.bat`.

### Visual Studio

Visual Studio (Community, Professional, or Enterprise) must be installed with the "Desktop development for C++" workload. For more information see the [installation docs](https://learn.microsoft.com/en-us/visualstudio/install/install-visual-studio?view=vs-2022).

When using Visual Studio 2022 and 2026, do ensure that 'Properties/libraries/Runtime Library -> Multithread DLL' and 'Properties/libraries/Use Intel Math Kernel Library -> Parallel (/Qmkl:parallel)' are turned on.

## Step 2: Download HAMS-MREL

You may obtain the HAMS-MREL source code using one of the following methods.

- **Clone the repository using Git**: if Git is installed, open a Command Prompt or PowerShell and run:
```
git clone https://github.com/vaibhavraghavan/HAMS-MREL.git
```

- **Download a ZIP archive**: Open the [GitHub repository](https://github.com/vaibhavraghavan/HAMS-MREL) in a web browser. Click "Code -> Download ZIP". Extract the contents of the ZIP file to a local directory of your choice.


## Step 3: Build and Execute HAMS-MREL using Visual Studio

Launch Visual Studio. On the welcome screen, select "Open a project or solution", then open the Visual Studio solution file [HAMS_MREL.sln](../HAMS_MREL.sln) located in the root directory of the HAMS-MREL repository.

In the top toolbar, set:
- **Configuration** to `Release`
- **Platform** to `x64`

Compile the application by selecting "Build -> Build Solution" from the top menu. The resulting executable `HAMS_MREL.exe` will be placed in the `HAMS-MREL/x64/Release` folder.

**Running**

To execute HAMS-MREL, click the **Start** button (green triangle) at the top, or select "Debug -> Start Without Debugging". A command prompt window will open and display the runtime output of the application.

> [!IMPORTANT]
> HAMS-MREL expects a set of Input Files to be located next to the executable. A full description of the files can be found [here](input_files.md). Output files will also be generated next to the executable.

---

## Alternative: Build Using gfortran (Command Line)

If you do not have access to the Intel oneAPI toolkit, HAMS-MREL can also be compiled using the free and open-source **GFortran** compiler together with the **LAPACK** and **BLAS** linear algebra libraries (as a substitute for Intel MKL).

### Step 1: Install GFortran

Download and install [MSYS2](https://www.msys2.org), which provides GFortran and associated libraries for Windows. After installation, open the **MSYS2 UCRT64** terminal and run:

```bash
pacman -Syu
pacman -S mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-lapack mingw-w64-ucrt-x86_64-openblas
```

This installs GFortran, LAPACK, and OpenBLAS (which provides BLAS). After installation, add the MSYS2 `ucrt64\bin` directory to your Windows PATH (typically `C:\msys64\ucrt64\bin`) so that `gfortran` is accessible from the Command Prompt or PowerShell.

To verify the installation, open a new Command Prompt and run:

```shell
gfortran --version
```

### Step 2: Compile HAMS-MREL

Open a Command Prompt or PowerShell, navigate to the `src` directory of the repository, and run:

```shell
cd <path-to-HAMS-MREL>\src

gfortran -O3 -fopenmp -ffree-line-length-none -o HAMS_MREL.exe ^
  WavDynMods.f90 InputFiles.f90 ^
  WavDynSubs.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 ^
  FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 ^
  CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 ^
  BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 ^
  HMatrix.f90 IterativeSolvers.f90 ^
  AssbMatx_irr.f90 AssbMatx.f90 ^
  AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 ^
  NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 ^
  PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 ^
  PressureElevationMulti.f90 SolveMotion.f90 SolveMotionMulti.f90 ^
  PrintOutput.f90 HAMS_Prog.f90 ^
  -llapack -lopenblas
```

The key flags are:

| Flag | Purpose |
|---|---|
| `-O3` | Optimisation level 3 (equivalent to ifort `/O3`) |
| `-fopenmp` | Enable OpenMP parallelisation (equivalent to ifort `-openmp`) |
| `-ffree-line-length-none` | Allow source lines longer than 132 characters |
| `-llapack -lopenblas` | Link LAPACK and BLAS (substitute for Intel MKL `/Qmkl`) |

> [!NOTE]
> `InputFiles.f90` must appear **before** the files that depend on `module IO`. The order shown above ensures all modules are compiled before they are used.

This produces the executable `hamsmrel.exe` inside the `src` directory.

### Step 3: Run HAMS-MREL

Copy `hamsmrel.exe` to the directory containing your `Input\` folder and run:

```shell
hamsmrel.exe
```

> [!IMPORTANT]
> The required runtime DLLs (`libgfortran`, `libgomp`, `libopenblas`, etc.) must be on your PATH or copied alongside the executable. If the application fails to launch with a missing DLL error, copy the relevant `.dll` files from `C:\msys64\ucrt64\bin\` into the same folder as `hamsmrel.exe`. You can identify the required DLLs by running `dumpbin /dependents hamsmrel.exe` in a Visual Studio Developer Command Prompt.
