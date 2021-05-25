
! ======================================================================
!
! This module contains data and parameters that are used by the GEM2
! evaporation and fission library.
!
! The module should, in good practice, only be used by the Evaporation
! class. Setup of initial data (2-D parameters) should only be called
! by the Evaporation module, however can be called by others.
!
! Written by CMJ, XCP-3, 8/2018 (Evap class creation)
!
! ======================================================================

module evaporationFissionData

  use, intrinsic :: iso_fortran_env, only: int32, real64
  use evaporationParams, only: zro, one, two

  implicit none
  private


  ! ('safety' to ensure arrays have been initialized)
  logical,        public,  protected :: evaporationDataEstablished = .FALSE.
  integer(int32), private, parameter :: noErrorsFound = 0_int32
  integer(int32), public,  parameter :: evaporationDataInitialized = 0_int32
  integer(int32), public,  parameter :: massDataFailed  =   1_int32
  integer(int32), public,  parameter :: shellDataFailed =  10_int32
  integer(int32), public,  parameter :: levelDataFailed = 100_int32


  ! Regarding level density parameterization sets [=0 (GCCI), =1 (a=A/8), >1 (a=a/alev)]:
  real(real64),   public,  parameter :: levelDensityGCCIFlag = 0.0_real64
  real(real64),   public,  parameter :: levelDensitySpecFlag = 1.0_real64


  ! Default variables for the data files
  character(LEN=*), public, parameter :: defaultMassFile  = "mass.tbl"
  character(LEN=*), public, parameter :: defaultLevelFile = "level.tbl"
  character(LEN=*), public, parameter :: defaultShellFile = "shell.tbl"


  ! For routine(s): gemdec, fprob
  real(real64),   public,  parameter    :: amu       = 931.4943d0
  integer(int32), public,  parameter    :: inn       = 150_int32
  integer(int32), public,  parameter    :: iiz       = 98_int32
  integer(int32), public,  parameter    :: aMax      = 300_int32
  integer(int32), public,  parameter    :: maxSize   = 70_int32

  ! Data from 'inigem' routine
  real(real64), protected, public, dimension(  3,   4)          :: t       = zro   ! Used by 'dost'
  real(real64), protected, public, dimension( :, :), allocatable:: &
       & exm,    & ! Used by 'gamma' (of size [maxSize, 200])
       & spin,   & ! Used by 'gamma'
       & width     ! Used by 'gamma'
  real(real64), protected, public, dimension(:, :), allocatable:: &
       & wapsm     ! Used by 'energy' (of size [0:150, 0:250])
  real(real64), protected, public, dimension(:, :), allocatable:: &
       & shellc    ! Used by 'efms' (of size [inn, 250])


  integer(int32) :: i
!     data tables of Cook et. al.  aaec/tm392, supplemented by 
!     Gilbert and Cameron.
  real(real64), public, parameter, dimension(iiz) :: sz = &
       & [   zro, zro, zro, zro, zro, zro, zro, zro, &
       &  -0.11d0,  -0.81d0,  -2.91d0,  -4.17d0,  -5.72d0, &
       &   -7.8d0,  -8.97d0,  -9.70d0, -10.10d0, -10.70d0, -11.38d0, &
       & -12.07d0, -12.55d0, -13.24d0, -13.93d0, -14.71d0, -15.53d0, &
       & -16.37d0, -17.36d0, -18.60d0, -18.70d0, -18.01d0, -17.87d0, &
       & -17.08d0, -16.60d0, -16.75d0, -16.50d0, -16.35d0, -16.22d0, &
       & -16.41d0, -16.89d0, -16.43d0, -16.68d0, -16.73d0, &
       & -17.45d0, -17.29d0, -17.44d0, -17.82d0, -18.62d0, -18.27d0, &
       & -19.39d0, -19.91d0, -19.14d0, -18.26d0, -17.40d0, -16.42d0, &
       & -15.77d0, -14.37d0, -13.91d0, -13.10d0, -13.11d0, -11.43d0, &
       & -10.89d0, -10.75d0, -10.62d0, -10.41d0, -10.21d0,  -9.85d0, &
       &  -9.47d0,  -9.03d0,  -8.61d0,  -8.13d0,  -7.46d0,  -7.48d0, &
       &  -7.20d0,  -7.13d0,  -7.06d0,  -6.78d0,  -6.64d0,  -6.64d0, &
       &  -7.68d0,  -7.89d0,  -8.41d0,  -8.49d0,  -7.88d0,  -6.30d0, &
       &  -5.47d0,  -4.78d0,  -4.37d0,  -4.17d0,  -4.13d0,  -4.32d0, &
       &  -4.55d0,  -5.04d0,  -5.28d0,  -6.06d0,  -6.28d0,  -6.87d0, &
       &  -7.20d0,  -7.74d0   ]

  real(real64), public, parameter, dimension(inn) :: sn = &
       ! Items 1 to 110
       & [          zro, zro, zro, zro, zro, zro, zro, zro, &
       &         10.3d0,   5.66d0,    6.8d0,   7.53d0, &
       &         7.55d0,   7.21d0,   7.44d0,   8.07d0,   8.94d0,   9.81d0, &
       &         10.6d0,  11.39d0,  12.54d0,  13.68d0,  14.34d0,  14.19d0, &
       &        13.83d0,  13.50d0,  13.00d0,  12.13d0,  12.60d0,  13.26d0, &
       &        14.13d0,  14.92d0,  15.52d0,  16.38d0,  17.16d0,  17.55d0, &
       &        18.03d0,  17.59d0,  19.03d0,  18.71d0,  18.80d0,  18.99d0, &
       &        18.46d0,  18.25d0,  17.76d0,  17.38d0,  16.72d0,  15.62d0, &
       &        14.38d0,  12.88d0,  13.23d0,  13.81d0,  14.90d0,  14.86d0, &
       &        15.76d0,  16.20d0,  17.62d0,  17.73d0,  18.16d0,  18.67d0, &
       &        19.69d0,  19.51d0,  20.17d0,  19.48d0,  19.98d0,  19.83d0, &
       &        20.20d0,  19.72d0,  19.87d0,  19.24d0,  18.44d0,  17.61d0, &
       &        17.10d0,  16.16d0,  15.90d0,  15.33d0,  14.76d0,  13.54d0, &
       &        12.63d0,  10.65d0,  10.10d0,   8.89d0,  10.25d0,   9.79d0, &
       &        11.39d0,  11.72d0,  12.43d0,  12.96d0,  13.43d0,  13.37d0, &
       &        12.96d0,  12.11d0,  11.92d0,  11.00d0,  10.80d0,  10.42d0, &
       &        10.39d0,   9.69d0,   9.27d0,   8.93d0,   8.57d0,   8.02d0, &
       &         7.59d0,   7.33d0,   7.23d0,   7.05d0,   7.42d0,   6.75d0, &
       &          6.6d0,    6.38d0, &
       ! Items 111 to 150
       &         6.36d0,   6.49d0,   6.25d0,   5.85d0, &
       &         5.48d0,   4.53d0,   4.3d0,    3.39d0,   2.35d0,   1.66d0, &
       &         0.81d0,   0.46d0,  -0.96d0,  -1.69d0,  -2.53d0,  -3.16d0, &
       &        -1.87d0,  -0.41d0,   0.71d0,   1.66d0,   2.62d0,   3.22d0, &
       &         3.76d0,   4.10d0,   4.46d0,   4.83d0,   5.09d0,   5.18d0, &
       &         5.17d0,   5.10d0,   5.01d0,   4.97d0,   5.09d0,   5.03d0, &
       &         4.93d0,   5.28d0,   5.49d0,   5.50d0,   5.37d0,   5.30d0   ]

  real(real64), public, parameter, dimension(iiz) :: pz = &
       & [        zro,  5.44d0, zro,  2.76d0, zro,  3.34d0, zro,  2.7d0, &
       &          zro,  1.90d0, zro,  2.12d0, zro,  2.09d0, zro,  1.62d0,  &
       &          zro,  1.62d0, zro,  1.83d0, zro,  1.73d0, zro,  1.35d0, &
       &          zro,  1.54d0,  zro,    1.28d0, 0.26d0, 0.88d0,  0.19d0, &
       &       1.35d0, -0.05d0, 1.52d0, -0.09d0, 1.17d0, 0.04d0,  1.24d0, &
       &       0.29d0,  1.09d0, 0.26d0,  1.17d0, 0.23d0, 1.15d0, -0.08d0, &
       &       1.35d0,  0.34d0, 1.05d0,  0.28d0, 1.27d0,  zro,    1.05d0, &
       &          zro, one,   0.09d0,  1.2d0,  0.2d0,  1.4d0, 0.93d0, one, &
       &       -0.2d0,  1.19d0, 0.09d0,  0.97d0, zro,    0.92d0, 0.11d0, &
       &       0.68d0,  0.05d0, 0.68d0, -0.22d0, 0.79d0, 0.09d0, 0.69d0, &
       &       0.01d0,  0.72d0,   zro,   0.4d0,  0.16d0, 0.73d0, zro, &
       &       0.46d0,  0.17d0,  0.89d0,   zro,  0.79d0,   zro,  0.89d0, &
       &          zro,  0.81d0, -0.06d0, 0.69d0, -0.2d0, 0.71d0, -0.12d0, &
       &       0.72d0,   zro,   0.77d0   ]

  real(real64), public, parameter, dimension(inn) :: pn = &
       ! Items 1 through 125
       & [        zro,  5.98d0,    zro,  2.77d0,    zro, &
       &       3.16d0,    zro,  3.01d0,    zro,  1.68d0,    zro,   1.73d0, &
       &         zro,   2.17d0,   zro,   1.67d0,    zro,  1.86d0,    zro, &
       &       2.04d0,    zro,  1.64d0,    zro,  1.44d0,    zro,   1.54d0, &
       &         zro,   1.3d0,    zro,   1.27d0,    zro,  1.29d0,  0.08d0, &
       &       1.41d0, -0.08d0,  1.5d0, -0.05d0, 2.24d0, -0.47d0,  1.43d0, &
       &      -0.15d0,  1.44d0, 0.06d0,  1.56d0, 0.25d0,  1.57d0, -0.16d0, &
       &       1.46d0,    zro,  0.93d0,  0.01d0, 0.62d0,  -0.5d0,  1.42d0, &
       &       0.13d0, 1.52d0, -0.65d0,  0.8d0, -0.08d0,  1.29d0, -0.47d0, &
       &       1.25d0, -0.44d0, 0.97d0,  0.08d0, 1.65d0, -0.11d0,  1.26d0, &
       &      -0.46d0,  1.06d0, 0.22d0,  1.55d0, -0.07d0, 1.37d0,  0.10d0, &
       &       1.20d0, -0.27d0, 0.92d0, -0.35d0, 1.19d0,    zro,   1.05d0, &
       &      -0.25d0,  1.61d0, -0.21d0, 0.9d0, -0.21d0,  0.74d0, -0.38d0, &
       &       0.72d0, -0.34d0, 0.92d0, -0.26d0, 0.94d0,  0.01d0,  0.65d0, &
       &      -0.36d0,  0.83d0, 0.11d0,  0.67d0, 0.05d0,  1.00d0,  0.51d0, &
       &       1.04d0,  0.33d0, 0.68d0, -0.27d0, 0.81d0,  0.09d0,  0.75d0, &
       &       0.17d0,  0.86d0, 0.14d0,  1.1d0, -0.22d0,  0.84d0, -0.47d0, &
       &       0.48d0,  0.02d0, 0.88d0,  0.24d0, 0.52d0,  0.27d0,  0.41d0, &
       &      -0.05d0, &
       ! Items 126 through 150
       &       0.38d0,  0.15d0, 0.67d0, zro, 0.61d0, zro, &
       &       0.78d0, zro, 0.67d0, zro, 0.67d0, zro, 0.79d0, zro, 0.6d0, &
       &       0.04d0, 0.64d0, -0.06d0, 0.45d0, 0.05d0, 0.26d0,-0.22d0, &
       &       0.39d0, zro, 0.39d0   ]

  real(real64), public, protected, dimension(:, :), allocatable:: &
       & st0, & ! (of sfize [iiz, inn])
       & paire0


  ! Julich mean value level density parameters
  real(real64), public, parameter, dimension(240) :: amean = &
       ! Items 1 through 100
       & [   0.125d0, 0.25d0, 0.375d0, 0.5d0, 0.625d0, &
       &     0.75d0,  0.875d0, 1.d0   , 1.125d0, 1.25d0,  1.375d0, 1.5d0,  &
       &    1.625d0,  1.75d0,  1.875d0, 2.d0,    2.125d0, 2.25d0, 2.375d0, &
       &     3.94d0,  2.63d0,  2.75d0,  2.88d0,  3.55d0,  4.35d0,  3.25d0, &
       &     3.38d0,  3.96d0,  3.63d0,  3.75d0,  3.88d0,  4.82d0,  4.44d0, &
       &     4.43d0,  4.43d0,  4.42d0,  4.63d0,  5.66d0,  5.81d0,  5.95d0, &
       &     5.49d0,  6.18d0,  7.11d0,  6.96d0,  7.2d0,   7.73d0,  6.41d0, &
       &     6.85d0,  6.77d0,  6.91d0,  7.26d0,  7.2d0,   6.86d0,  8.06d0, &
       &     7.81d0,  7.82d0,  8.41d0,  8.13d0,  7.19d0,  8.35d0,  8.13d0, &
       &     8.02d0,  8.93d0,  8.9d0,   9.69d0,  9.65d0, 10.55d0,  9.38d0, &
       &     9.72d0, 10.66d0, 11.98d0, 12.76d0, 12.10d0, 12.86d0, 13.03d0, &
       &    12.81d0, 12.54d0, 12.65d0, 12.00d0, 12.69d0, 14.05d0, 13.33d0, &
       &    13.28d0, 13.23d0, 13.17d0,  8.66d0, 11.09d0, 10.40d0, 13.47d0, &
       &    10.17d0, 12.22d0, 11.62d0, 12.95d0, 13.15d0, 13.57d0, 12.87d0, &
       &    16.16d0, 14.71d0, 15.69d0, 14.09d0, &
       ! Items 101 through 200
       &   18.56d0, 16.22d0, 16.67d0, 17.13d0, &
       &     17.d0,  16.86d0, 15.33d0, 15.61d0, 16.77d0, 17.93d0, 17.45d0, &
       &    16.97d0, 17.88d0, 17.58d0, 15.78d0, 16.83d0, 17.49d0, 16.03d0, &
       &    15.08d0, 16.74d0, 17.74d0, 17.43d0, 18.14d0, 17.06d0, 19.01d0, &
       &    17.02d0, 17.02d0, 17.02d0, 18.51d0, 17.20d0, 16.75d0, 16.97d0, &
       &    16.94d0, 16.91d0, 17.69d0, 15.55d0, 14.56d0, 14.35d0, 16.55d0, &
       &    18.29d0, 17.80d0, 17.05d0, 21.31d0, 19.15d0, 19.51d0, 19.87d0, &
       &    20.39d0, 20.90d0, 21.85d0, 22.89d0, 25.68d0, 24.64d0, 24.91d0, &
       &    23.24d0, 22.85d0, 22.46d0, 21.98d0, 21.64d0, 21.75d0, 21.85d0, &
       &    21.77d0, 21.69d0, 23.74d0, 21.35d0, 23.03d0, 20.66d0, 21.81d0, &
       &    20.77d0, 22.18d0, 22.58d0, 22.55d0, 21.45d0, 21.16d0, 21.02d0, &
       &    20.87d0, 22.09d0, 22.00d0, 21.28d0, 23.05d0, 21.70d0, 21.45d0, &
       &    22.28d0, 23.00d0, 22.11d0, 23.56d0, 22.83d0, 24.88d0, 22.64d0, &
       &    23.27d0, 23.89d0, 23.92d0, 23.94d0, 21.16d0, 22.3d0,  21.75d0, &
       &    21.19d0, 20.72d0, 20.24d0, 21.34d0, 19.00d0, &
       ! Items 201 through 240
       &    17.93d0, 17.85d0, 15.70d0, 13.54d0, &
       &    11.78d0, 10.02d0, 10.98d0, 10.28d0, 11.72d0, 13.81d0, 14.46d0, &
       &    15.30d0, 16.15d0, 16.99d0, 17.84d0, 18.68d0, 19.53d0, 20.37d0, &
       &    21.22d0, 22.06d0, 22.91d0, 23.75d0, 24.60d0, 25.44d0, 26.29d0, &
       &    27.13d0, 27.98d0, 28.82d0, 29.67d0, 30.71d0, 30.53d0, 31.45d0, &
       &    29.63d0, 30.15d0, 30.65d0, 30.27d0, 29.52d0, 30.08d0, 29.80d0, &
       &    29.87d0   ]

!   Flags for deformed nuclei:
  real(real64), public, parameter, dimension(2) :: con = [  0.142d0, 0.12d0  ]
  logical,      public, parameter, dimension(iiz) :: isz = [ &
       & (.false., i= 1, 53 ), &
       & (.true.,  i=54, 78 ), &
       & (.false., i=79, 85 ), &
       & (.true.,  i=86, 98 )   ]
  logical,      public, parameter, dimension(inn) :: isn = [ &
       & (.false., i=  1,  85 ), &
       & (.true.,  i= 86, 122 ), &
       & (.false., i=123, 129 ), &
       & (.true.,  i=130, 150 )   ]

! ======================================================================
! The following (omega, ifa, ifz) contain information on the 66 ejectiles
! that are considered for evaporation in the model.
! ======================================================================
!               n     p      d      t     He3   He4   He6   He8
!              Li6   Li7    Li8    Li9
!              Be7   Be9    Be10   Be11   Be12  
!              B8    B10    B11    B12    B13
!              C10   C11    C12    C13    C14   C15   C16
!              N12   N13    N14    N15    N16   N17
!              O14   O15    O16    O17    O18   O19   O20
!              F17   F18    F19    F20    F21
!             Ne18   Ne19   Ne20   Ne21   Ne22  Ne23  Ne24
!             Na21   Na22   Na23   Na24   Na24
!             Mg22   Mg23   Mg24   Mg25   Mg26  Mg27  Mg28
! ====================================================================== 
! **omega=spin
  real(real64), public, parameter, dimension(maxSize) :: omega = [ &
       &       0.5d0, 0.5d0, 1.0d0, 0.5d0, 0.5d0, 0.0d0, 0.d0, 0.d0, &
       &       1.0d0, 1.5d0, 2.0d0, 1.5d0, &
       &       1.5d0, 1.5d0, 0.0d0, 0.5d0, 0.0d0, &
       &       2.0d0, 3.0d0, 1.5d0, 1.0d0, 1.5d0, &
       &       0.0d0, 1.5d0, 0.0d0, 0.5d0, 0.0d0, 0.5d0, 0.d0, &
       &       1.0d0, 0.5d0, 1.0d0, 0.5d0, 2.0d0, 0.5d0, &
       &       0.0d0, 0.5d0, 0.0d0, 2.5d0, 0.0d0, 2.5d0, 0.d0, &
       &       2.5d0, 1.0d0, 0.5d0, 2.0d0, 2.5d0, &
       &       0.0d0, 0.5d0, 0.0d0, 1.5d0, 0.0d0, 2.5d0, 0.d0, &
       &       1.5d0, 3.0d0, 1.5d0, 4.0d0, 2.5d0, &
       &       0.0d0, 1.5d0, 0.0d0, 2.5d0, 0.0d0, 0.5d0, 0.d0, &
       &       0.d0,  0.d0,  0.d0,  0.d0   ]

  !  ifa = mass number
  integer(int32), public, parameter, dimension(maxSize) :: ifa = [ &
       &         1,     1,     2,     3,     3,     4,     6,     8,   &
       &         6,     7,     8,     9, &
       &         7,     9,    10,    11,    12, &
       &         8,    10,    11,    12,    13, &
       &        10,    11,    12,    13,    14,    15,    16, &
       &        12,    13,    14,    15,    16,    17,  &
       &        14,    15,    16,    17,    18,    19,    20, &
       &        17,    18,    19,    20,    21, &
       &        18,    19,    20,    21,    22,    23,    24, &
       &        21,    22,    23,    24,    25, &
       &        22,    23,    24,    25,    26,    27,    28,   &
       &         0,     0,     0,     0   ]

  !  ifz = atomic number (charge)
  integer(int32), public, parameter, dimension(maxSize) :: ifz = [ &
       &   0            , &
       & ( 1, i=  2,  4), &
       & ( 2, i=  5,  8), &
       & ( 3, i=  9, 12), &
       & ( 4, i= 13, 17), &
       & ( 5, i= 18, 22), &
       & ( 6, i= 23, 29), &
       & ( 7, i= 30, 35), &
       & ( 8, i= 36, 42), &
       & ( 9, i= 43, 47), &
       & (10, i= 48, 54), &
       & (11, i= 55, 59), &
       & (12, i= 60, 66), &
       & ( 0, i= 67, 70)   ]
! ======================================================================

  ! From 'auxeye' common block previously (not dependent on 'alev' calc. option)
  real(real64), public, parameter, dimension(maxSize)  :: omfct  = &
       & [ ( (two * omega(i) + one),   i = 1, maxSize) ]
  real(real64), public, parameter, dimension(aMax)     :: ux0    = &
       & [ ( (2.5d0 + 150.d0/dble(i)), i = 1, aMax) ]
  real(real64), public, parameter, dimension(aMax)     :: ux140  = &
       & [ ( ( sqrt(sqrt( ux0(i) )) ), i = 1, aMax) ]

  ! Data used by 'gamma' routine
  real(real64), public, protected, dimension(aMax) :: alpp  = zro
  real(real64), public, protected, dimension(aMax) :: betap = zro


  public :: initializeEvaporationData

contains

! ======================================================================

  include "initializeData.f90"

! ======================================================================
end module evaporationFissionData
