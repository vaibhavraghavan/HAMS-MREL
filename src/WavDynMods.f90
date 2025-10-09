!  ------------------------------------------------------------------------------------------------------
!                                                               
!    Program HAMS for the diffraction and radiation of waves 
!    by 3D structures.
! 
!  License:
! 
!    This routine is part of HAMS.
!
!    HAMS is a free software framework: you can redistribute it and/or modify it 
!    under the terms of the Apache License, Version 2.0 (the "License"); you may 
!    not use these subroutines except in compliance with the License. The software
!    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
!
!    You should have received a copy of the Apache License, Version 2.0, along with 
!    HAMS. If not, see <http://www.apache.org/licenses/LICENSE-2.0>.
!
!  Code Original Author:
!
!    Yingyi Liu
!
!  ------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!        Data module for declaring input variables in HAMS (This is for variables for all the modules which are implemented in HAMS. These are hence mostly public.)
!---------------------------------------------------------------------------------------------					
!
        MODULE HAMS_mod
!
        INTEGER,PUBLIC:: IRSP,INFT,OUFT,SYBO,ISOL
        INTEGER,PUBLIC:: NTHREAD                                 ! The number of threads to be used in openmpi (parallel)

        REAL*8,PUBLIC::  WK1,DWK,BETA1,DBETA

        END MODULE HAMS_mod

        
!---------------------------------------------------------------------------------------------
!        Data module for declaring wave relevant variables in Morison_stick numerical model: Is this for the calculation of forces?
!---------------------------------------------------------------------------------------------	
!
        MODULE WaveDyn_mod
!
        INTEGER,PUBLIC:: NPER,NBETA,NK                                  
        PARAMETER(NK=200)
        REAL*8,PUBLIC:: A,H,AMP,V,BETA,BLNR(6,6),BQDR(6,6)                           ! BLNR - External linear damping matrix, BQDR - External quadratic damping matrix, V is the wave number in deep water, H is the dimensional finite water depth,
        REAL*8,PUBLIC:: WK,W1,WL,TP,WVN(NK),INFR,OUFR                                ! WVN - Real, an array storing roots of the dispersion equation, WK - Positive root of the dispersion equation, W1 - Current wave frequency, WL - Wave length corresponding to WK, Time period corresponding to the current wave frequency, INFR - Input frequency format, OUFR - Output frequency format
!
        REAL*8,ALLOCATABLE,PUBLIC:: WVNB(:),WVFQ(:),WVHD(:)                          
        REAL*8,ALLOCATABLE,PUBLIC:: AMAS(:,:,:),BDMP(:,:,:)                          ! Mass and damping matrix respectively

        COMPLEX*16,ALLOCATABLE,PUBLIC:: EXFC(:,:,:),DSPL(:,:,:)                      ! Excitation force matrix and Displacement matrix respectively
        
        ! Variables for multi-body interactions
        
        INTEGER,PUBLIC:: NBODY                                                       ! Number of bodies to be considered for multi-body interaction
        REAL*8,ALLOCATABLE,PUBLIC:: LCS_MULTI(:,:)                                   ! Array of location of the Local Coordinate System (LCS) along with the rotation of the LCS w.r.t the Global Coordinate System (GCS) 
        REAL*8,ALLOCATABLE,PUBLIC:: AMAS_MULTI(:,:,:,:),BDMP_MULTI(:,:,:,:),AMAS_MULTI_COMB(:,:,:),BDMP_MULTI_COMB(:,:,:)          ! Mass and damping matrix respectively
        COMPLEX*16,ALLOCATABLE,PUBLIC:: EXFC_MULTI(:,:,:,:), EXFC_MULTI_COMB(:,:,:),DSPL_MULTI(:,:,:,:),DSPL_MULTI_COMB(:,:,:)      ! Excitation force matrix and Displacement matrix respectively
        REAL*8,ALLOCATABLE,PUBLIC:: BLNR_MULTI(:,:,:),BQDR_MULTI(:,:,:)                          ! BLNR - External linear damping matrix, BQDR - External quadratic damping matrix
        END MODULE WaveDyn_mod

        
!---------------------------------------------------------------------------------------------
!        Data module for declaring body relevant variables in Morison_stick numerical model
!---------------------------------------------------------------------------------------------	
!
        MODULE Body_mod
        
        REAL*8,PUBLIC:: XR(3),XG(3),XB(3),XW(2)
        REAL*8,PUBLIC:: VOL,MASS,REFL
!
        REAL*8,PUBLIC:: IB(3,3),MATX(6,6),CRS(6,6),KSTF(6,6),RAO(6)
        
        ! Variables for multi-body interactions
        
        REAL*8,ALLOCATABLE,PUBLIC:: XR_MULTI(:,:),XG_MULTI(:,:),XB_MULTI(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: VOL_MULTI(:),MASS_MULTI(:)                                            ! As for the single bodies, the Reference length (REFL) is not added here, since that is expected to be the same for all bodies.
        REAL*8,ALLOCATABLE,PUBLIC:: IB_MULTI(:,:,:),MATX_MULTI(:,:,:),CRS_MULTI(:,:,:),KSTF_MULTI(:,:,:),RAO_MULTI(:,:)
!
        END MODULE Body_mod
        
        
!---------------------------------------------------------------------------------------------
!        Data module for declaring variables in panel method
!---------------------------------------------------------------------------------------------
!
        MODULE PanelMesh_mod
!
        INTEGER NELEM,NTND
        INTEGER ISYS,NSYS,ISX,ISY
!
        INTEGER,ALLOCATABLE,PUBLIC:: NCN(:),NCON(:,:),NCOND(:,:)
        INTEGER,ALLOCATABLE,PUBLIC:: IPIV(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: XYZ(:,:),DXYZ_P(:,:),XYZ_P(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: DS(:),PNSZ(:)
        
        ! Variables for multi-body interactions (Add the description of the variables)
        
        INTEGER,ALLOCATABLE:: NELEM_MULTI(:), NTND_MULTI(:)                         ! Number of elements and number of nodes 
        INTEGER,ALLOCATABLE,PUBLIC:: NCN_MULTI(:,:),NCON_MULTI(:,:,:), NCN_MULTI_COMB(:),NCON_MULTI_COMB(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: XYZ_LOCAL_MULTI(:,:,:),XYZ_GLOBAL_MULTI(:,:,:),DXYZ_MULTI_P(:,:,:),XYZ_MULTI_P(:,:,:),XYZ_GLOBAL_MULTI_COMB(:,:),XYZ_GLOBAL_MULTI_COMB_P(:,:),DXYZ_MULTI_COMB(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: DS_MULTI(:,:),PNSZ_MULTI(:,:), PNSZ_MULTI_COMB(:), DS_MULTI_COMB(:)
        INTEGER,ALLOCATABLE,PUBLIC:: IPIV_MULTI_COMB(:,:)
        
        ! Variables for removal of irregular frequencies (The variables with i and _COMB are for only water plane mesh panels)
        
        INTEGER,ALLOCATABLE,PUBLIC:: iNELEM_MULTI(:), iNTND_MULTI(:)                  ! Number of elements and number of nodes in the water plane mesh
        INTEGER,ALLOCATABLE,PUBLIC:: iNCN_MULTI(:,:),iNCON_MULTI(:,:,:),iNCN_MULTI_COMB(:),iNCON_MULTI_COMB(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: iXYZ_LOCAL_MULTI(:,:,:),iXYZ_GLOBAL_MULTI(:,:,:),iDXYZ_MULTI_P(:,:,:),iXYZ_MULTI_P(:,:,:),iXYZ_GLOBAL_MULTI_COMB(:,:),iXYZ_GLOBAL_MULTI_COMB_P(:,:),iDXYZ_MULTI_COMB(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: iDS_MULTI(:,:),iPNSZ_MULTI(:,:),iDS_MULTI_COMB(:),iPNSZ_MULTI_COMB(:)
        INTEGER,ALLOCATABLE,PUBLIC:: iNELEM_TOTAL,NELEM_TOTAL,iNTND_TOTAL,NTND_TOTAL
        
!
        END MODULE PanelMesh_mod
        
!---------------------------------------------------------------------------------------------
!        Data module for declaring variables of general computation
!---------------------------------------------------------------------------------------------						
!
        MODULE Const_mod
!
        REAL*8,PUBLIC:: G,RHO,PI
        REAL*8,PUBLIC:: RXY(2,2),RX(2,2),RY(2,2),SY(2,2),SX(2,2)
        REAL*8,PUBLIC:: RXY_MULTI(2,2),RX_MULTI(2,2),RY_MULTI(2,2),SY_MULTI(2,2),SX_MULTI(2,2)
!
        COMPLEX*16,PUBLIC:: CI
!
        DATA G,PI,RHO/9.80665D0,3.141592653589793D0, 1000.D0/                  !G is mostly acceleration due to gravity, PI is the mathematical pi, RHO/9.81 is conversion to N, Actual pi value, 1025 - Default density of sea water 
        DATA CI/(0.0D0, 1.0D0)/
        DATA RXY  /  1.D0,  1.D0,  1.D0, -1.D0  /
        DATA RX   /  1.D0, -1.D0,  1.D0,  1.D0  /
        DATA RY   /  1.D0,  1.D0,  1.D0, -1.D0  /
        DATA SX   /  1.D0,  1.D0,  1.D0,  1.D0  /
        DATA SY  /   1.D0, -1.D0, -1.D0,  1.D0  /
        
        ! For Multi-bodies, symmetry is implemented
        DATA RXY_MULTI  /  1.D0,  1.D0,  1.D0,  -1.D0  /
        DATA RX_MULTI   /  1.D0,  -1.D0,  1.D0,  1.D0  /
        DATA RY_MULTI   /  1.D0,  1.D0,  1.D0,  -1.D0  /
        DATA SX_MULTI   /  1.D0,  1.D0,  1.D0,  1.D0  /
        DATA SY_MULTI  /   1.D0,  -1.D0,  -1.D0,  1.D0  /

        END MODULE Const_mod
        

!---------------------------------------------------------------------------------------------
!        Data module for declaring variables for removing irregular frequencies
!---------------------------------------------------------------------------------------------			
!
        MODULE Inerfs_mod
!
        INTEGER,PUBLIC:: iNELEM,iNTND,tNELEM,tNTND,NELEM_PE
!     
        INTEGER,ALLOCATABLE,PUBLIC:: iNCN(:),iNCON(:,:),iNCOND(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: iXYZ(:,:),iDXYZ_P(:,:),iXYZ_P(:,:)
        REAL*8,ALLOCATABLE,PUBLIC:: iDS(:),iPNSZ(:)
!
        END MODULE Inerfs_mod
    
!---------------------------------------------------------------------------------------------
!        Data module for declaring variables in linear algeraic system - This is for solving the final linear algebraic equations
!---------------------------------------------------------------------------------------------						
!
        MODULE LinearMatrix_mod
!
        COMPLEX*16,ALLOCATABLE,PUBLIC:: AMAT(:,:,:),BRMAT(:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: BDMAT(:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: CMAT(:,:,:),DRMAT(:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: DDMAT(:,:)
        
        ! Variable to be used for multi-body interaction
        COMPLEX*16,ALLOCATABLE,PUBLIC:: AMAT_MULTI(:,:,:),BRMAT_MULTI(:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: BDMAT_MULTI(:,:)
        
        ! Variables for irregular frequency removal
        COMPLEX*16,ALLOCATABLE,PUBLIC:: CMAT_MULTI(:,:,:),DRMAT_MULTI(:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: DDMAT_MULTI(:,:)
        
!
        END MODULE LinearMatrix_mod

!---------------------------------------------------------------------------------------------
!        Data module for declaring variables for potentials
!---------------------------------------------------------------------------------------------						
!
        MODULE Potentials_mod
!
        COMPLEX*16,ALLOCATABLE,PUBLIC:: MXPOT(:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: CGRN(:,:,:,:),RKBN(:,:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: DGRN(:,:,:,:),PKBN(:,:,:,:)
        
        ! Varible for multi-body interaction
        
        COMPLEX*16,ALLOCATABLE,PUBLIC:: CGRN_MULTI_COMB(:,:,:,:),RKBN_MULTI_COMB(:,:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: MXPOT_MULTI_COMB(:,:,:)
        COMPLEX*16,ALLOCATABLE,PUBLIC:: DGRN_MULTI_COMB(:,:,:,:),PKBN_MULTI_COMB(:,:,:,:)
!
        END MODULE Potentials_mod


!---------------------------------------------------------------------------------------------
!        Data module for declaring variables for field points
!---------------------------------------------------------------------------------------------
!
        MODULE FieldOutput_mod
        !USE Precision
!
        INTEGER NFP
!
        REAL*8,ALLOCATABLE,PUBLIC:: XFP(:,:)
!
        END MODULE FieldOutput_mod
        
