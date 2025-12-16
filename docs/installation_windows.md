# Windows Installation Instructions

This document describes how to set up a development environment that allows you to compile and run HAMS-MREL on the Windows.

## Step 1: Install Dependencies

HAMS-MREL requires the following software:

### Intel Fortran Compiler and oneMK Library

Download and install the [Intel Fortran Essentials](https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-windows/2025-0/intel-fortran-essentials.html) bundle (includes the oneAPI MKL library) from the [Intel website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?packages=fortran-essentials&fortran-essentials-os=windows&fortran-essentials-win=offline).

There are two types of types of installers: offline and online. Both can be installed with the command line (Command Prompt or Windows PowerShell) or via a GUI. If installing with the command line, make sure you activate the integration with Visual Studio by using the `NEED_VS<yyy>_INTEGRATION` property. For example, if Working With Visual Studio 2022:

```shell
./intel-fortran-essentials-2025.1.0.555_offline.exe -a --silent --eula accept -p=NEED_VS2022_INTEGRATION=1
```

**Important**
After installation you must run the `setvars.bat` script to add the libraries to your system PATH. This ensures the Intel compiler and libraries are discoverable by Visual Studio. The script is typically located in `C:\Program Files (x86)\Intel\oneAPI\setvars.bat`.

### Visual Studio

Visual Studio (Community, Professional, or Enterprise) must be installed with the "Desktop development for C++" workload. For more information see the [installation docs](https://learn.microsoft.com/en-us/visualstudio/install/install-visual-studio?view=vs-2022).

## Step 2: Download HAMS-MREL

You may obtain the HAMS-MREL source code using one of the following methods.

- **Clone the repository using Git**: if Git is installed, open a Command Prompt or PowerShell and run:
```
git clone https://github.com/vaibhavraghavan/HAMS-MREL-DCC.git
```

- **Download a ZIP archive**: Open the [GitHub repository](https://github.com/vaibhavraghavan/HAMS-MREL-DCC) in a web browser. Click "Code -> Download ZIP". Extract the contents of the ZIP file to a local directory of your choice.


## Step 3: Build and Execute HAMS-MREL using Visual Studio

Launch Visual Studio. On the welcome screen, select "Open a project or solution", then open the Visual Studio solution file [HAMS_original.sln](../HAMS_original.sln) located in the root directory of the HAMS-MREL repository.

In the top toolbar, set:
- **Configuration** to `Release`
- **Platform** to `x64`

Compile the application by selecting "Build -> Build Solution" from the top menu. The resulting executable `HAMS_original.exe` will be placed in the `HAMS-MREL-DCC/x64/Release` folder.

**Running**

To execute HAMS-MREL, click the **Start** button (green triangle) at the top, or select "Debug -> Start Without Debugging". A command prompt window will open and display the runtime output of the application.

> [!IMPORTANT]
> HAMS-MREL expects a set of Input Files to be located next to the executable. A full description of the files can be found [here](input_files.md). Output files will also be generated next to the executable.