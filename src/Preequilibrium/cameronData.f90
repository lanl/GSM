
! ======================================================================
!
! ----------------------------------------------------------------------
!      Data for calculating nuclear masses and inverse cross sections  *
!  see: tabs. 96-97 and formulae (6.58a-6.58b), (6.79-6.80) and p. 422 *
!  in [1]                                                              *
!   Removed "OLD" binding energies for emitted particles dlmo; AJS 10/03
!   Removed Cameron shell corrections; 12/04/03.
!
! ----------------------------------------------------------------------
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by A. J. Sierk, LANL T-16, December, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
!
! LMK 07/2012
! Particle of type "j" for aj and zj:
! 1 = n
! 2 = p
! 3 = d
! 4 = t
! 5 = He3
! 6 = He4
! 7 = He6
! 8 = He8
! 9 = Li6
! 10 = Li7
! 11 = Li8
! 12 = Li9
! 13 = Be7
! 14 = Be9
! 15 = Be10
! 16 = Be11
! 17 = Be12
! 18 = B8
! 19 = B10
! 20 = B11
! 21 = B12
! 22 = B13
! 23 = C10
! 24 = C11
! 25 = C12
! 26 = C13
! 27 = C14
! 28 = C15
! 29 = C16
! 30 = N12
! 31 = N13
! 32 = N14
! 33 = N15
! 34 = N16
! 35 = N17
! 36 = O14
! 37 = O15
! 38 = O16
! 39 = O17
! 40 = O18
! 41 = O19
! 42 = O20
! 43 = F17
! 44 = F18
! 45 = F19
! 46 = F20
! 47 = F21
! 48 = Ne18
! 49 = Ne19
! 50 = Ne20
! 51 = Ne21
! 52 = Ne22
! 53 = Ne23
! 54 = Ne24
! 55 = Na21
! 56 = Na22
! 57 = Na23
! 58 = Na24
! 59 = Na25
! 60 = Mg22
! 61 = Mg23
! 62 = Mg24
! 63 = Mg25
! 64 = Mg26
! 65 = Mg27
! 66 = Mg28
!
! ======================================================================

     ! A of the allowed fragments
     real(real64), public, parameter, dimension(numAllowedFragments) :: aj = [ &
         & 1.d0,      1.d0,     2.d0,     3.d0,     3.d0,     4.d0, &
         & 6.d0,      8.d0,     6.d0,     7.d0,     8.d0,     9.d0, &
         & 7.d0,      9.d0,    10.d0,    11.d0,    12.d0,     8.d0, &
         & 10.d0,    11.d0,    12.d0,    13.d0,    10.d0,    11.d0, &
         & 12.d0,    13.d0,    14.d0,    15.d0,    16.d0,    12.d0, &
         & 13.d0,    14.d0,    15.d0,    16.d0,    17.d0,    14.d0, &
         & 15.d0,    16.d0,    17.d0,    18.d0,    19.d0,    20.d0, &
         & 17.d0,    18.d0,    19.d0,    20.d0,    21.d0,    18.d0, &
         & 19.d0,    20.d0,    21.d0,    22.d0,    23.d0,    24.d0, &
         & 21.d0,    22.d0,    23.d0,    24.d0,    25.d0,    22.d0, &
         & 23.d0,    24.d0,    25.d0,    26.d0,    27.d0,    28.d0   ]
 
    ! A^(1/3) for the allowed fragments
    integer(int32), private :: i
    real(real64), public, parameter, dimension(numAllowedFragments) :: ajthr = &
         & [ (ato3rd(nint(aj(i))), i = 1, numAllowedFragments) ]

    ! Z of the allowed fragments
    real(real64), public, parameter, dimension(numAllowedFragments) :: zj = [ &
         &  0.d0,    1.d0,   1.d0,   1.d0,   2.d0,   2.d0, &
         &  2.d0,    2.d0,    3.d0,    3.d0,    3.d0,    3.d0, &
         &  4.d0,    4.d0,    4.d0,    4.d0,    4.d0,    5.d0, &
         &  5.d0,    5.d0,    5.d0,    5.d0,    6.d0,    6.d0, &
         &  6.d0,    6.d0,    6.d0,    6.d0,    6.d0,    7.d0, &
         &  7.d0,    7.d0,    7.d0,    7.d0,    7.d0,    8.d0, &
         &  8.d0,    8.d0,    8.d0,    8.d0,    8.d0,    8.d0, &
         &  9.d0,    9.d0,    9.d0,    9.d0,    9.d0,    10.d0, &
         & 10.d0,    10.d0,    10.d0,    10.d0,    10.d0,    10.d0, &
         & 11.d0,    11.d0,    11.d0,    11.d0,    11.d0,    12.d0, &
         & 12.d0,    12.d0,    12.d0,    12.d0,    12.d0,    12.d0   ]

!    real(real64), public, parameter, dimension( 6) :: dlmo = [ &
!         & 8.368d0, 7.569d0, 13.835d0, 15.835d0, 15.817d0, 3.607d0   ]

    real(real64), public, parameter, dimension(numAllowedFragments) :: dlmn = [ &
         &  8.071d0,   7.289d0,  13.136d0,  14.950d0,  14.931d0,   2.425d0, &
         & 17.590d0,  31.600d0,  14.090d0,  14.910d0,  20.950d0,  24.950d0, &
         & 15.770d0,  11.350d0,  12.610d0,  20.170d0,  25.080d0,  22.920d0, &
         & 12.050d0,   8.670d0,  13.370d0,  16.560d0,  15.700d0,  10.650d0, &
         &  0.000d0,   3.130d0,   3.020d0,   9.870d0,  13.690d0,  17.340d0, &
         &  5.350d0,   2.860d0,   0.100d0,   5.680d0,   7.870d0,   8.010d0, &
         &  2.860d0,  -4.740d0,  -0.810d0,  -0.780d0,   3.330d0,   3.800d0, &
         &  1.950d0,   0.870d0,  -1.490d0,  -0.020d0,  -0.050d0,   5.320d0, &
         &  1.750d0,  -7.040d0,  -5.730d0,  -8.030d0,  -5.160d0,  -5.950d0, &
         & -2.190d0,  -5.180d0,  -9.530d0,  -8.420d0,  -9.360d0,  -0.400d0, &
         & -5.470d0, -13.930d0, -13.190d0, -16.210d0, -14.590d0, -15.020d0   ]

    ! NOTE: These values are NOT used anywhere
    real(real64), public, parameter, dimension( 5) :: &
         & z1 = [  10.d0,     20.d0,    30.d0,    50.d0,    70.d0    ], &
         & a1 = [   0.42d0,    0.58d0,   0.68d0,   0.77d0,   0.80d0  ], &
         & c1 = [   0.50d0,    0.28d0,   0.20d0,   0.15d0,   0.10d0  ], &
         & a2 = [   0.68d0,    0.82d0,   0.91d0,   0.97d0,   0.98d0  ], &
         & c2 = [   0.10d0,    0.10d0,   0.10d0,   0.08d0,   0.06d0  ]

! ======================================================================
