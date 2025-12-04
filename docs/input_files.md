# Input Files

HAMS-MREL requires a set of input files to run. The number files depend on whether you are analyzing a single body or multiple bodies:
- `ControlFile.in`
- `Hydrostatic[_i].in`
- `HullMesh[_i].pnl`
- `WaterPlaneMesh[_i].pnl` (only if analyzing irregular frequencies)

For multi-body cases, `_i` denotes the body number (`i = 1 … numberofbodies`). For a single body, no index is used (`HullMesh.pnl`, `Hydrostatic.in` and `WaterPlaneMesh.pnl`).

The program will check that all required files exist before running. If any file is missing, execution will stop.

### Mesh Files

As described by the [HAMS user manual](https://github.com/YingyiLiu/HAMS/blob/master/Manual/A%20brief%20manual.pdf), to obtain the water and body meshes in a HAMS-compatible format, you can start by exporting the hydrodynamic CAD model from Rhinoceros to WAMIT’s *.gdf format. Next, use the built-in tool WAMIT_MeshTran.exe to convert the *.gdf mesh into the HAMS mesh format. To do this, run RunWAMIT_MeshTran.bat and enter the filename of the *.gdf file when prompted. Ensure that the *.gdf file includes meshes for both the waterplane and the submerged bodies. The tool will automatically separate them into two files: WaterplaneMesh.pnl and HullMesh.pnl. Note that the WAMIT_MeshTran tool is currently only available on Windows.

### Fortran Format Descriptors

HAMS-MREL uses Fortran format descriptors to read values from the file. These are indicated in the Format column in the tables below:
- `x`: skip spaces. For example, `14x` means skip 14 characters.
- `iN`: integer with width N. Example: `i16` reads a 16-character wide integer.
- `fW.D`: floating-point number with width W and D digits after the decimal point. Example: `f30.15` reads a number in 30 characters with 15 digits after the decimal.
- `A100`: alphanumeric string of length 100.

**Example:** `20x,f30.15` means skip the first 20 characters of the line, then read a floating-point number using 30-character width and 15 digits after the decimal point.

**Example:** `FORMAT(6(2x,E12.5))` reads 6 numbers per line, each in 12-character wide scientific notation with 5 digits after the decimal, preceded by 2 spaces.

## ControlFile.in

The general settings of the simulation are specified in `ControlFile.in`. The file is read by the subroutine `ReadControlFile` in [InputFiles.f90](../src/InputFiles.f90). The format for the file is as follows.

**Water Depth and Wave Frequencies**

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| 1 | free | - | Comment or blank line (ignored) |
| 2 | free | - | Comment or blank line (ignored) |
| 3 | `14x,f30.15` | `H` | Water depth (m). <0 for infinite depth, >0 for finite depth |
| 4 | free | - | Comment or blank line |
| 5 | free | - | Comment or blank line |
| 6 | `27x,i16` | `SYB0` | Switch: 0 or 1. 1 to include zero and infinite frequency added mass calculations; 0 to exclude them.  |
| 7 | `25x,i16` | `INFT` | Input frequency type: 1 for deepwater wave number, 2 for finite depth wave number, 3 for wave frequency, 4 for wave period and 5 for wave length. 25 spaces before the integer |
| 8 | `25x,i17` | `OUFT` | Output frequency type; same options as input frequency type |
| 9 | `26x,i16` | `NPET` | Number of frequencies |
| if `NPET > 0` | space-separated list of length `NPET` (single line) | `WVNB` | Discrete wave frequencies, periods, numbers or lengths, depending on `Input_frequency_type` (see [example](https://github.com/YingyiLiu/HAMS/blob/master/CertTest/Moonpool/Input/ControlFile.in))
| if `NPET < 0` | `27x,f30.15` | `WK1` | Minimum frequency
| if `NPET < 0` | `19x,f30.15` | `DWK` | Frequency step

**Body Information**

As shown in the table below, the first two lines after the last value from the previous section are ignored.

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| next | `A100` | `BODY_CHECK` | String indicating whether multiple bodies are used. Include `"multi"` if multiple bodies. |
| if `BODY_CHECK` contains "multi": | `24x,i16` | `NBODY` | Number of bodies; must be ≥1 |
| next `NBODY` lines | `26x,4f12.3` |`LCS_MULTI` | origin coordinates of the local coordinate system (LCS), per body mesh

**Wave Headings**

Wave incident angles with respect to the x-direction. All wave heading input values are expressed in degrees.

As shown in the table below, the first three lines after the last value from the previous section are ignored.

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| next | `23x,i16` | `NBETA` | Number of wave headings; a negative sign indicates values are generated from a minimum and a step size rather than listed explicitly
| if `NBETA > 0` | space-separated list of length `NBETA` | `WVHD` | `NBETA` wave headings |
| if `NBETA < 0` | `20x,f30.15` | `BETA1` | minimum heading
| if `NBETA < 0` | `17x,f30.15` | `DBETA` | heading step

**Rotation Centers**

As shown in the table below, the first two lines after the last value from the previous section are ignored.

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| if `NBODY=1` | `28x,3f12.3` | `XR(1:3)` | Rotation center coordinates |
| if `NBODY>1` | `28x,3f12.3` | `XR_MULTI(NBODY,1:3)` | Each body's rotation center coordinates; one set of coordinates per line |

**Solver Options**

This section begines on the next line after the previous one. 

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | `26x,f30.15` | `REFL` | Reference body length |
| next | `26x,i16` | `ISOL` | Wave diffraction solution. Switch: 1 or 2. 1 if incident and scattering potentials are combined; 2 if only scattering.  |
| next | `23x,i16` | `IRSP` | Irregular frequencies: 1 for to remove them, 0 to keep them |
| next | `23x,i16` | `NTHREAD` | Number of OpenMP threads; setting it to 1 removes parallelization |

**Pressure and Elevation**

As shown in the table below, the first two lines after the last value from the previous section are ignored.

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| next | `27x,i16` | `NFP` | Number of field points to output field pressure and elevation. |
| next `NFP` lines | free | `XFP(I,1:3)` | Coordinates of field points (1 set per line) |

## Hydrostatic.in

Hydrostatic and damping properties are specified in the `Hydrostatic[_i].in` files. They are read by the subroutine `ReadHydroStaticMulti` in [HydroStaticMulti.f90](../src/HydroStaticMulti.f90). There is one file per body.


| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| 1 | free | - | Comment or blank line (ignored) |
| 2 | `6(2x,E12.5)` | `XG_MULTI(BODY_N,1:3)` | Center of gravity coordinates (x, y, z) |
| 3 | free | - | Comment or blank line (ignored) |
| 4-9 | `6(2x,E12.5)` | `MATX_MULTI(BODY_N,1:6,1:6)` | 6×6 mass matrix, one row per line |
| 10 | free | - | Comment or blank line (ignored) |
| 11-16 | `6(2x,E12.5)` | `BLNR_MULTI(BODY_N,1:6,1:6)` | 6×6 external linear damping matrix, one row per line |
| 17 | free | - | Comment or blank line (ignored) |
| 18-23 | `6(2x,E12.5)` | `BQDR_MULTI(BODY_N,1:6,1:6)` | 6×6 external quadratic damping matrix, one row per line |
| 24 | free | - | Comment or blank line (ignored) |
| 25-30 | `6(2x,E12.5)` | `CRS_MULTI(NBODY,1:6,1:6)` | 6×6 hydrostatic restoring matrix, one row per line |
| 31 | free | - | Comment or blank line (ignored) |
| 32-37 | `6(2x,E12.5)` | `KSTF_MULTI(NBODY,1:6,1:6)` | 6×6 external restoring matrix, one row per line |


## HullMesh.pnl

Node and element data on the surface of each body are specified in the `HullMesh[_i].pnl` files. There is one file per body.

The files are read in [HAMS_Prog.f90](../src/HAMS_Prog.f90) and in the subroutine `ReadBodyMeshMulti` [ReadPanelMeshMulti.f90](../src/ReadPanelMeshMulti.f90). The expected file structure is described below.

A few notes:
- Numbers in each line can be separated by whitespaces but **not commas**.
- Any trailing text on each line must be preceded by `!`.
- The number of vertex indices must match `NTND` exactly.

| Line(s) | Format | Variable | Description |
| :------ | :------ | :-------- | :----------- |
| 1–3 | free | — | Comment or blank lines (ignored) |
| 4 | free | `NELEM, NTND, ISX, ISY` | Number of panels, number of nodes, X-symmetry flag, Y-symmetry flag |
| 5-6 | free | — | Comment or blank lines (ignored) |
| next `NTND` | free | `M` and `XYZ_LOCAL_MULTI(bodyID,nodeID,1:3)` | Node coordinates. Each line contains: `node_ID  x_local  y_local  z_local` |
| next 3 lines | free | — | Comment or blank (ignored) |
| next `NTND` lines | free | `M`, `NCN_MULTI(bodyID,nodeID)` and `NCON_MULTI(bodyID,nodeID,1:numvertices)`| Panel definitions. Each line contains: panel number, number of vertices and (vertex1, vertex2, ...)


## WaterPlaneMesh.pnl

Node and element data on the inner water plane for multiple bodies are specified in the `WaterPlaneMesh[_i].pnl` files. There is one file per body.

The files are read in [HAMS_Prog.f90](../src/HAMS_Prog.f90) and in the subroutine `ReadWTPLMeshMulti` in [ReadPanelMeshMulti.f90](../src/ReadPanelMeshMulti.f90). The file structure is exactly the same as the one for the HullMesh.pnl files, which is described above.

The node coordinates of the waterplane mesh are stored in `iXYZ_LOCAL_MULTI(bodyID,nodeID,1:3)`. The number of vertices in each panel is stored in `iNCN_MULTI(bodyID,nodeID)` and the vertices (vertex1, vertex2, ...) are stored in `iNCON_MULTI(bodyID,nodeID,1:numvertices)`.