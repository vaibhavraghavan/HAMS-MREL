%ComSpec% /E:ON /K ""C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019"
rem  HAMS-MREL build flags
rem    /Qmkl          - link Intel MKL (provides multi-threaded ZGETRF/ZGETRS)
rem    /O3            - high-level optimisation
rem    /heap-arrays0  - REQUIRED. Forces all automatic arrays and compiler-generated array temporaries
rem                     onto the heap. Without this, multi-body cases with TNELEM > ~3000 abort with
rem                     `forrtl: severe (170): Program Exception - stack overflow` on the first
rem                     frequency iteration, because the default 1-MB Windows thread stack is too
rem                     small for the array sections created inside the inlined CALGREEN_*/SGLINTBD_*
rem                     call chain. The `0` argument means "no threshold - every automatic array".
rem    -openmp        - enable OpenMP parallelisation
rem  Optional, opt-in for advanced users:
rem    /QxHost            - generate the best instruction set for the build machine
rem                          (AVX/AVX2/AVX-512). Ties the binary to the build machine's ISA -
rem                          for distribution builds use /arch:AVX2 instead.
rem    /align:array64byte - align arrays to 64 bytes for AVX-512 cache lines.
rem  Do NOT enable /fp:fast - it relaxes complex-arithmetic IEEE rules and breaks bit-identity on
rem    the Green's-function kernels.
ifort /Qmkl /O3 /heap-arrays0 /exe:HAMS_ifort_Win.exe WavDynMods.f90 InputFiles.f90 WavDynSubs.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 AssbMatx_irr.f90 AssbMatx.f90 AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 SolveMotion.f90 SolveMotionMulti.f90 PrintOutput.f90 HAMS_Prog.f90 -openmp
pause