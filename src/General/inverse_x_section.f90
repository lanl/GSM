! Leslie M. Kerby, 04/2014, LANL

module inverse_x_section

! ==============================================================================
! Calculates inverse cross section the following ways:
!    1) Current/Old preequilibrium (Original Dostrovsky[1] with custom Vc)
!    2) Current GEM2 (Modified Dostrovsky[1])
!    3) NASA [2, 3, 4] - Kalbach [5]
!
! DEFAULT is ___NASA-Kalbach___ 
!
! [1] Dostrovsky et al., "Monte Carlo Calculations of Nuclear Evaporation Processes. III. Applications to Low-Energy Reactions",
!    Physical Review, 116 (1959) 683.
! [2] Tripathi et al., "Accurate universal parameterization of absorption cross sections," NIM B 117 (1996) 347.
! [3] Tripathi et al., "^... - neutron absorption cross sections," NIM B 129 (1997) 11.
! [4] Tripathi et al., "Universal Parameterization of Absorption Cross Sections: Light Systems," NASA/TP-1999-209726.
! [5] Kalbach, C., "Toward a global exciton model; lessons at 14 MeV," J. Phys. G: Nucl. Part. Phys. 24 (1998) 847.
! ==============================================================================

  use, intrinsic :: iso_fortran_env, only: real64
  use fund_data, only: r_rms
  use coulomb_barrier, only: cb_nasa

  implicit none
  private
! SUBROUTINE declarations
  public  :: orig_dost        ! Calculates inverse cross section as in current/old preequilibrium (orig. Dostrovsky)
  public  :: orig_dost_init   ! Initializes parameters for orig. Dostrovsky calculation
  public  :: mod_dost         ! Calculates inverse cross section as in current GEM2 evaporation (mod. Dostrovsky)
  public  :: mod_dost_init    ! Initializes parameters for mod. Dostrovsky calculation
  public  :: nasa             ! Calculates inverse cross section by NASA model
  private :: kalbach

! Dostrovsky
  real (real64) :: cj(66) = 1.0            ! Dostrovsky's original parametrization cj, 1=n, 2=p, 3=d...
  real (real64) :: aj(66) = 1.0            ! "alpha" for neutron channel and (1+cj) for charged particle channel


! ==============================================================================

contains

! ==============================================================================

  subroutine orig_dost(part_type, z_res_before,a_res_before,energy_frag, &
       & z_frag, a_frag, coul_bar, x_sec)

! ==============================================================================
! The inverse cross section calculation as found in the old/current
!    preequilibrium with the exception that "am" = 0.1
!
! Calculation is Original Dostrovsky except for Vc
! Coulomb Barrier Vc = vk(l)*(1.439976/radncl)*zj(l)*zfj(l)/(ajthr(l) +
!    afjthr(l))*(1 – u/(81*a*am))
!
! x_sec = 70.68375*(a_res_before)**(2./3.)*aj*(1-coul_bar/energy_frag)    (mb)
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    real (real64), intent(in)  :: z_res_before ! Z and A of residual nucleus BEFORE emission of fragment
    real (real64), intent(in)  :: a_res_before ! Z and A of residual nucleus BEFORE emission of fragment
    real (real64), intent(in)  :: z_frag       ! Z and A of the emitted particle
    real (real64), intent(in)  :: a_frag       ! Z and A of the emitted particle
    real (real64), intent(in)  :: coul_bar     ! Coulomb barrier (MeV)
    real (real64), intent(in)  :: energy_frag  ! Energy of emitted fragment (MeV)
    real (real64), intent(out) :: x_sec        ! Inverse cross section (mb)
 
    integer (int32) :: part_type   ! particle type (from 1-66)

! ==============================================================================

    ! Calculate inverse cross section
    if (energy_frag > coul_bar) then
       x_sec = ( 70.68375 * (a_res_before)**(2./3.) ) * aj(part_type) * &
            & ( 1. - coul_bar/energy_frag )   ! (mb) PiR^2 = 70.68375*(a_res_before)**(2./3.)
    else 
       x_sec = 0.0
    end if

    return
! ==============================================================================
  end subroutine orig_dost


  subroutine orig_dost_init(z_res_before,a_res_before)

! ==============================================================================
! Linear interpolation used in the old preequilibrium code. Lagrange interpolation used here.
! Original Dostrovsky parameters:
!    z    cp    calpha    
!    10    0.50    0.10
!    20    0.28    0.10
!    30    0.20    0.10
!    50    0.15    0.08
!    70    0.10    0.06
!    (100    0.08    0.05)    (used to have better extrapolation)
! Lagrange polynomial: 
! cj(2)= cp = −1.13095e−9*z^5 +3.57738e−7*z^4 −4.322e−5*z^3 +2.50e−3*z^2 −0.071781*z + 1.0075
! cj(6)= calpha = −3.439e−10*z^5 +8.690e−8*z^4 −7.6415e−6*z^3 +2.7218e−4*z^2 −4.01336e−3*z + 0.11972
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none
    ! Z and A of residual nucleus BEFORE emission of fragment
    real (real64), intent(in) :: z_res_before
    real (real64), intent(in) :: a_res_before

! ==============================================================================

    ! aj = "alpha" for neutron channel
    aj(1) = 0.76+2.2*a_res_before**(-1./3.)
    ! aj = (1+cj) for charged particle channels
    cj(2) = -1.13095e-9*(z_res_before**5) +3.57738e-7*(z_res_before**4) - &
         & 4.322e-5*(z_res_before**3) + 2.50e-3*(z_res_before**2) - &
         & 0.071781*z_res_before + 1.0075
    cj(6) = -3.439e-10*(z_res_before**5) +8.690e-8*(z_res_before**4) - &
         & 7.6415e-6*(z_res_before**3) + 2.7218e-4*(z_res_before**2) - &
         & 4.01336e-3*z_res_before + 0.11972
    cj(3) = cj(2)/2.
    cj(4) = cj(2)/3.
    cj(5) = (4./3.)*cj(6)
    cj(7:66) = 0.0
    aj(2:66) = (1. + cj(2:66))

    return
! ==============================================================================
  end subroutine orig_dost_init


  subroutine mod_dost(z_res_before,a_res_before,energy_frag,num_preeq_part, &
       & z_type, a_type, coul_bar, x_sec)

! ==============================================================================
! The inverse cross section calculation as found in GEM2
!
! x_sec = geom_x_sec*aj*(1-coul_bar/energy_frag)    (mb)
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    real (real64),   intent(in)  :: z_res_before   ! Z of residual nucleus BEFORE emission of fragment 
    real (real64),   intent(in)  :: a_res_before   ! A of residual nucleus BEFORE emission of fragment 
    real (real64),   intent(in)  :: energy_frag    ! Energy of emitted fragment (MeV)
    integer (int32), intent(in)  :: num_preeq_part ! Number of preequilibrium particles (1-66)
    real (real64),   intent(in)  :: z_type(66)     ! Z for the particle types 1-66, used for emitted particle
    real (real64),   intent(in)  :: a_type(66)     ! A for the particle types 1-66, used for emitted particle
    real (real64),   intent(in)  :: coul_bar(66)   ! Coulomb barrier (MeV)
    real (real64),   intent(out) :: x_sec(66)      ! Inverse cross section (mb)


    integer (int32) :: i
    real (real64) :: a_res_after(66)=1.0    ! A of residual nucleus AFTER emission of fragment
    real (real64) :: geom_x_sec(66)=1.0     ! Geometric cross section (mb), 10*Pi*radius^2 (1 b = 100 fm^2, 10 mb = 1 fm^2))
    real (real64) :: radius(66)=1.0         ! Radius for geom_x_sec (fm): 1.5A^(1/3) for n,p; 1.5(Ad^(1/3) + Aj^(1/3)) for d-He4, &
                                            ! (1.12*(A_d**(1./3.)-0.86*(A_d)**(-1./3.)) + (1.12*(A_f)**(1./3.)-0.86*(A_f)**(-1./3.)) + 2.85 for > He4

 ! ==============================================================================

    do i=1,num_preeq_part  
       if (z_res_before < z_type(i) .or. a_res_before <= a_type(i)) then 
          exit  
       end if
       a_res_after(i) = a_res_before - a_type(i)
    end do

    ! Calculate inverse cross section
    radius(1:2) = 1.5*a_res_after(1:2)**(1./3.)   ! (fm)
    radius(3:6) = 1.5*(a_res_after(3:6)**(1./3.) + a_type(3:6)**(1./3.))   !(fm)
    radius(7:66) = (1.12*(a_res_after(7:66))**(1./3.)-0.86 * &
         & (a_res_after(7:66))**(-1./3.)) + (1.12*(a_type(7:66))**(1./3.) - &
         & 0.86*(a_type(7:66))**(-1./3.)) + 2.85   !(fm)
    geom_x_sec = 10.0*3.141593*radius**2   !(mb)
    do i =1,num_preeq_part
       if (energy_frag > coul_bar(i)) then
          x_sec(i) = geom_x_sec(i)*aj(i)*(1-coul_bar(i)/energy_frag)   ! (mb) 
!          x_sec(i) = aj(i)*(1-coul_bar(i)/energy_frag)   ! Gives ratio => cross section/geom. cross section
       else 
          x_sec(i) = 0.0
       end if
    end do

    return
! ==============================================================================

  end subroutine mod_dost


! ==============================================================================

  subroutine mod_dost_init(z_res_before,a_res_before)

! ==============================================================================
! Linear interpolation used in the GEM2 code. Exponential regression used here.
! Modified Dostrovsky parameters:
!    z    cp    calpha    
!    20    0.00    0.0
!    30    -0.06    0.0
!    40    -0.10    0.0
!    50    -0.10    0.0
!
! Exponential regression: 
! cj(2)= cp = 0.295519*exp(-0.0613123*Z) - 0.11
! cj(6)= calpha = 0.0
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none
    real (real64), intent(in) :: z_res_before    ! Z of residual nucleus BEFORE emission of fragment
    real (real64), intent(in) :: a_res_before    ! A of residual nucleus BEFORE emission of fragment

! ==============================================================================

    ! aj = "alpha" for neutron channel
    aj(1) = 0.76+1.93*a_res_before**(-1./3.)
    ! aj = (1+cj) for charged particle channels
    cj(2) = 0.295519*exp(-0.0613123*z_res_before) - 0.11
    cj(6) = 0.0
    cj(3) = cj(2)/2.
    cj(4) = cj(2)/3.
    cj(5) = (4./3.)*cj(6)
    cj(7:66) = 0.0
    aj(2:66) = (1. + cj(2:66))

    return
! ==============================================================================
  end subroutine mod_dost_init


  subroutine nasa(z_res_after, a_res_after, z_frag, a_frag, t_a_frag, t_cm, coul_bar, x_sec)

! ==============================================================================
! Calculates inverse cross section according to NASA parameterization
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none
    real (real64), intent(in)  :: z_res_after ! Z and A of residual nucleus AFTER emission of fragment
    real (real64), intent(in)  :: a_res_after ! Z and A of residual nucleus AFTER emission of fragment
    real (real64), intent(in)  :: z_frag      ! Z and A of emitted fragment
    real (real64), intent(in)  :: a_frag      ! Z and A of emitted fragment
    real (real64), intent(in)  :: t_a_frag    ! Kinetic Energy of emitted fragment, in lab frame (A MeV)
    real (real64), intent(in)  :: t_cm        ! Kinetic energy in center-of-momentum frame (MeV)
    real (real64), intent(out) :: coul_bar    ! Coulomb barrier (MeV)
    real (real64), intent(out) :: x_sec      ! Inverse cross section (mb)   (10 mb = 1 fm^2)

    real (real64) :: delta   ! Delta_E, the energy dependence of the cross section (mainly due to transparency and Pauli blocking)
    real (real64) :: s       ! S, mass asymmetry term
    real (real64) :: ce      ! C_E, related to transparency and Pauli blocking
    real (real64) :: d       ! D, related to density dependence of colliding system - scaled to C+C density dependence
    real (real64) :: xm      ! X_m, modifies neutron cross section according to imaginary part of optical potential at surface
    real (real64) :: x1      ! X_1, used to calculate Xm
    real (real64) :: sl      ! S_L, used to calculate Xm
    real (real64) :: t1      ! T_1, parameter related to surface of colliding system
    real (real64) :: dg      ! D_g, used to calculate D for neutrons
    real (real64) :: rc      ! R_c, coulomb multiplier for light systems
    real (real64) :: g       !! G, parameter for alpha-nucleus systems

    ! Kalbach stuff
    real(real64), parameter, dimension(300) :: enmax = [   &
         & 14.55, 14.97, 55.64, 6.32, 6.30, 6.24, 6.46, 6.24, 6.36, 5.95, &  ! A_res_after=2...
         & 12.94, 8.86, 8.85, 8.70, 8.69, 8.63, 8.64, 9.03, 8.99, 9.01, & 
         & 8.97, 8.87, 8.80, 8.66, 8.72, 8.82, 8.81, 8.47, 8.47, 8.34, & 
         & 8.33, 8.21, 8.20, 8.08, 8.07, 7.96, 7.94, 7.83, 8.15, 8.75, & 
         & 8.61, 8.59, 8.47, 8.43, 8.41, 8.30, 8.27, 8.17, 8.14, 8.09, & 
         & 7.32, 7.29, 7.21, 7.17, 7.15, 7.07, 7.05, 6.97, 6.71, 6.67, & ! A_res=52...
         & 6.63, 6.59, 6.23, 6.20, 6.17, 6.46, 6.42, 6.40, 6.38, 6.34, & 
         & 6.32, 6.30, 6.28, 6.25, 6.24, 6.22, 6.21, 6.19, 6.17, 6.17, & 
         & 6.16, 6.15, 6.14, 6.14, 6.13, 6.12, 6.12, 6.12, 6.12, 6.11, & 
         & 6.12, 6.14, 6.14, 6.15, 6.16, 6.19, 6.18, 6.21, 6.23, 6.24, & 
         & 6.26, 6.25, 6.31, 6.33, 6.36, 6.38, 6.42, 6.45, 6.48, 6.52, & ! A_res=102
         & 6.55, 6.58, 6.62, 6.67, 6.71, 6.74, 6.80, 6.85, 6.89, 6.95, & 
         & 7.00, 7.05, 7.10, 7.16, 7.21, 7.28, 7.33, 7.40, 7.46, 7.53, & 
         & 7.60, 7.68, 7.74, 7.82, 7.89, 7.97, 8.04, 8.12, 8.21, 8.59, & 
         & 8.67, 8.78, 8.85, 8.96, 9.06, 9.13, 9.25, 9.33, 9.46, 9.57, & 
         & 9.66, 9.77, 9.87, 9.98,10.12,10.22,10.34,10.47,10.55,10.71, & ! A_res=152
         & 10.83,10.96,11.08,11.22,11.31,11.47,11.63,11.75,11.90,12.06, &
         & 12.16,12.35,12.48,12.64,12.80,12.94,13.11,13.25,13.42,13.60, &
         & 13.71,13.91,14.06,14.23,14.44,14.60,14.78,14.96,15.14,15.31, &
         & 15.52,15.70,15.88,16.10,16.30,16.47,16.66,16.86,15.42,15.45, &
         & 15.47,15.53,19.97,19.98,20.02,25.02,24.94,25.06,25.05,25.02, & ! A_res=202
         & 25.10,25.05,25.14,25.10,25.08,25.16,25.14,25.20,25.18,25.15, &
         & 25.22,25.21,25.28,25.27,25.22,25.29,25.28,25.33,25.31,25.29, &
         & 25.38,25.34,25.39,25.38,25.36,25.44,25.40,25.47,25.45,25.43, &
         & 25.48,25.45,25.53,25.50,25.49,25.55,25.52,25.58,25.56,25.55, &
         & 25.58,25.58,25.66,25.63,25.60,25.66,25.63,25.68,25.67,25.65, & ! A_res=252...
         & 25.68,25.67,25.77,25.71,25.71,25.77,25.77,25.82,25.79,25.79, &
         & 25.85,25.82,25.85,25.85,25.84,25.88,25.88,25.88,25.88,25.88, &
         & 25.96,25.93,25.96,25.96,25.94,25.99,25.97,26.03,26.01,25.99, &
         & 26.08,26.04,26.10,26.08,26.06,26.12,26.08,26.14,26.13,26.08   ]

    real(real64), parameter, dimension(300) :: scalf = [   &
         & 0.000,0.6757,0.5922,1.111,1.057,1.099,1.144,1.122,1.130,1.089, & 
         & 0.9642,1.062,1.045,1.055,1.040,1.050,1.059,1.047,1.031,1.040, &
         & 1.050,1.038,1.026,1.031,1.037,1.028,1.016,1.017,1.019,1.008, &
         & 0.9962,0.9965,0.9967,0.9851,0.9856,0.9841,0.9812,0.9680,0.9675,1.028, &
         & 1.031,1.034,1.037,1.040,1.043,1.046,1.048,1.050,1.052,1.054, &
         & 1.049,1.051,1.052,1.054,1.055,1.056,1.057,1.058,1.051,1.052, & ! A_res = 52...
         & 1.053,1.054,1.054,1.055,1.055,1.058,1.059,1.061,1.061,1.062, &
         & 1.063,1.064,1.065,1.066,1.066,1.067,1.068,1.068,1.069,1.070, &
         & 1.070,1.071,1.071,1.072,1.073,1.073,1.074,1.074,1.074,1.074, &
         & 1.074,1.074,1.074,1.074,1.074,1.074,1.074,1.073,1.073,1.073, & 
         & 1.073,1.073,1.073,1.072,1.072,1.072,1.072,1.071,1.071,1.070, & ! A_res = 102...
         & 1.070,1.069,1.069,1.069,1.068,1.067,1.067,1.066,1.066,1.065, &
         & 1.064,1.064,1.063,1.063,1.062,1.061,1.060,1.060,1.059,1.058, &
         & 1.057,1.056,1.055,1.054,1.053,1.052,1.051,1.050,1.048,1.047, &
         & 1.046,1.045,1.044,1.042,1.041,1.033,1.039,1.037,1.036,1.035, &
         & 1.033,1.031,1.030,1.029,1.028,1.026,1.025,1.027,1.022,1.021, & ! A_res = 152...
         & 1.020,1.018,1.017,1.017,1.015,1.013,1.012,1.012,1.009,1.008, &
         & 1.007,1.005,1.004,1.002,1.002,0.9997,0.9992,0.9980,0.9968,0.9956, &
         & 0.9944,0.9930,0.9920,0.9911,0.9899,0.9888,0.9868,0.9856,0.9845,0.9838, &
         & 0.9827,0.9815,0.9799,0.9788,0.9776,0.9764,0.9752,0.9741,0.8902,0.8902, &
         & 0.8903,0.8906,0.8992,0.8993,0.8990,0.9041,0.9042,0.9748,0.9740,0.9576, & ! A_res=202
         & 0.9585,0.9591,0.9600,0.9591,0.9583,0.9577,0.9571,0.9565,0.9704,0.9708, &
         & 0.9701,0.9695,0.9689,0.9683,0.9544,0.9538,0.9532,0.9680,0.9672,0.9526, &
         & 0.9660,0.9527,0.9661,0.9656,0.9647,0.9516,0.9638,0.9659,0.9652,0.9646, &
         & 0.9522,0.9516,0.9511,0.9506,0.9500,0.9505,0.9510,0.9517,0.9522,0.9516, & 
         & 0.9511,0.9513,0.9520,0.9513,0.9507,0.9514,0.9518,0.9513,0.9508,0.9512, & ! A_res=252..
         & 0.9518,0.9513,0.9509,0.9513,0.9518,0.9524,0.9519,0.9514,0.9508,0.9502, &
         & 0.9508,0.9513,0.9518,0.9512,0.9507,0.9502,0.9507,0.9512,0.9507,0.9502, &
         & 0.9497,0.9501,0.9495,0.9490,0.9485,0.9491,0.9486,0.9481,0.9476,0.9480, &
         & 0.9486,0.9490,0.9484,0.9478,0.9472,0.9481,0.9476,0.9471,0.9478,0.9472   ]
    ! end Kalbach

! ==============================================================================

    ! Calculate energy-dependent coulomb barrier
    call cb_nasa(z_res_after, a_res_after, z_frag, a_frag, t_cm, coul_bar)

    ! Assign X_m, T_1, D, and R_c for various reactions
    if (nint(z_frag) == 0) then
! *********** neutrons ************
       if (a_res_after > 2. .and. a_res_after < 300.) then
          if (t_a_frag < enmax(nint(a_res_after) - 1) ) then
             ! Calculate cross section by Kalbach systematics for low-energy neutrons
             x_sec = kalbach(1, 0.0_real64, 1.0_real64, z_res_after, a_res_after, t_a_frag) * &
                  & scalf(nint(a_res_after)-1)
             return
          end if
       end if
       rc = 1.0
       ! Assign X_m
       if (nint(a_res_after) < 200 ) then ! A_t < 200
          if (nint(z_res_after) == 2 .and. nint(a_res_after) == 4) then
             ! He4
             x1 = 5.2
          else
             x1 = 2.83 - 0.031*(a_res_after) + 0.00017*(a_res_after)**2
          end if
          if (nint(a_res_after) <= 4) then
             ! Light systems
             sl = 1.2 + 1.6*(1.0 - exp(-t_a_frag/15.))
          else if (nint(a_res_after) < 12) then
             sl = 0.6
          else if (nint(a_res_after) == 12) then
             sl = 1.6
          else
             sl = 1.0
          end if
          xm = 1. - x1*exp(-t_a_frag/(x1*sl))
       else
          ! A_t >= 200
          xm = (1. - 0.3*exp(-(t_a_frag - 1.)/15.))*(1. - exp(0.9-t_a_frag))
       end if
       ! Assign T_1
       if (nint(a_res_after) <= 4) then
          ! Light systems
          t1 = 18.
       else if (nint(a_res_after) >= 11 .and. nint(a_res_after) <= 40) then
          ! 11 <= A_t <= 40
          t1 = 30.
       else
          t1 = 40.0
       end if
       ! Assign D
       dg = 0.538/( 0.23873*a_frag/( 1.29 * r_rms( nint(z_frag+1.0), &
            & nint(a_frag) ) )**3 + 0.23873 * a_res_after/( 1.29 * &
            & r_rms( nint(z_res_after+1.0), nint(a_res_after) )  )**3 )
       if (nint(a_res_after) <= 4) then
          ! Light systems
          d = 1.85 + 0.16/(1.0 + exp((500. - t_a_frag)/200.) )
       else if (nint(a_res_after) <= 40 ) then
          ! A_t <=40
          d = dg - 1.5*(a_res_after - 2.*z_res_after)/a_res_after + 0.25 / &
               & ( 1. + exp((t_a_frag - 170.)/100.))
       else if (nint(a_res_after) > 40 .and. nint(a_res_after) < 60) then
          ! 40 < A_t < 60
          d = dg - 1.5*(a_res_after - 2.*z_res_after)/a_res_after
       else if (nint(z_res_after) > 82) then
          ! Z_t > 82
          d = dg - z_res_after/(a_res_after - z_res_after)
       else
          d = dg
       end if
    else 
! *********** charged particles ***********
       xm = 1.0
       rc = 1.0
       t1 = 40.
       d = 1.75 * ( 0.23873 * a_frag / ( 1.29 * &
            & r_rms( nint(z_frag+1.0), nint(a_frag) )  )**3 + &
            & 0.23873 * a_res_after / ( 1.29 * &
            & r_rms( nint(z_res_after+1.0), nint(a_res_after) ) )**3 ) / &
            & ( 2.*0.23873*12.0/(1.29*r_rms(7,12))**3 )
       if (nint(a_frag) == 1 .and. nint(z_frag) == 1) then
          ! protons
          if (nint(a_res_after) <= 4) then
             ! p + Light system
             t1 = 23.
             d = 1.85 + 0.16/(1. + exp((500. - t_a_frag)/200.) )
             if (nint(a_res_after) == 2 .and. nint(z_res_after) == 1) then
                ! p + d
                rc = 13.5
             else if (nint(a_res_after) == 3 .and. nint(z_res_after) == 2) then
                ! p + 3He
                rc = 21.
             else if (nint(a_res_after) == 4 .and. nint(z_res_after) == 2) then
                ! p + 4He
                rc = 27.
             end if
          else
             d = 2.05
             if (nint(z_res_after) == 3) then
                ! p + Li
                rc = 3.5
             else if (nint(z_res_after) == 6) then
                ! p + C
                rc = 3.5
             end if
          end if
       else if (nint(a_frag) == 2 .and. nint(z_frag) == 1) then
          ! deuterons
          if (nint(a_res_after) <= 4) then
             ! d + Light system
             t1 = 23.
             d = 1.65 + 0.1/(1. + exp((500. - t_a_frag)/200.) )
             if (nint(a_res_after) == 2 .and. nint(z_res_after) == 1) then
                ! d + d
                rc = 13.5
             else if (nint(a_res_after) == 4 .and. nint(z_res_after) == 2) then
                ! d + 4He
                rc = 13.5
             end if
          else
             if (nint(z_res_after) == 6) then
                ! d + C
                rc = 6.0
             end if
          end if
       else if (nint(a_frag) == 3 .and. nint(z_frag) == 2) then
          ! 3He
          if (nint(a_res_after) <= 4) then
             ! 3He + Light system
             d = 1.55
          end if
       else if (nint(a_frag) == 4 .and. nint(z_frag) == 2) then
          ! 4He
          g = 75.
          if (nint(a_res_after) == 4 .and. nint(z_res_after) == 2) then
             ! 4He + 4He
             g = 300.
          else if (nint(z_res_after) == 4) then
             ! 4He + Be
             t1 = 25.
             g = 300.
          else if (nint(z_res_after) == 7) then
             ! 4He + N
             g = 500.
          else if (nint(z_res_after) == 13) then
             ! 4He + Al
             t1 = 25.
             g = 300.
          else if (nint(z_res_after) == 26) then
             ! 4He + Fe
             g = 300.
          else if (nint(z_res_after) == 73) then
             ! 4He + Ta
             rc = 0.6
          else if (nint(z_res_after) == 79) then
             ! 4He + Au
             rc = 0.6
          end if
          d = 2.77 - 0.008*a_res_after + 0.000018*(a_res_after)**2 - 0.8 / &
               & (1. + exp((250. - t_a_frag)/g))
       else if (nint(z_frag) == 3) then
          ! Li
          d = d/3.
       end if
    end if

! Calculate Delta_E
    if (t_cm < 1.0d-15) then
       write(*,1000) "203"
       x_sec = 1.0
       return
    end if

    ce = d*(1.0 - exp(-t_a_frag/t1)) - 0.292*exp(-t_a_frag/792.) * &
         & cos( 0.229 * (t_a_frag**0.453) )
    s = a_frag**(1./3.)*a_res_after**(1./3.) / ( a_frag**(1./3.) + &
         & a_res_after**(1./3.) )
    delta = 1.85*s + 0.16*s/t_cm**(1./3.) - ce + 0.91 * &
         ( a_res_after - 2.*z_res_after ) * z_frag/(a_res_after * a_frag)

    ! Calculate inverse cross section
    x_sec = 3.14159*(1.21)*(10.0) * ( a_frag**(1./3.) + a_res_after**(1./3.) &
         & + delta)**2 * (1. - rc*coul_bar/t_cm) * xm   ! (mb)
    x_sec = max(0.0_real64, x_sec)

    return
! ==============================================================================
1000 format(3x, "Warning: Divide by zero error detected in ", &
          & "'inverse_x_section.f90', subroutine 'nasa', line ", A, ".")
! ==============================================================================
  end subroutine nasa


  function kalbach(kp, zp, ap, zt, at, elab)

! ==============================================================================
!
!     written in 1982; revised 1990 by K.K. Gudima
!
!     calculate optical model reaction cross sections
!     using empirical parameterization
!     of Narasimha Murthy, Chaterjee, and Gupta
!     going over to the geometrical limit at high energy
!
!           proton cross sections scaled down with signor for a<100
!           (appropriate for Becchetti-Greenlees potential)
!           neutron cross sections scaled down sith signor for a<40
!           (appropriate for Mani et al potential)
!
!     parameter values set in subroutine sigpar
! 
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    integer(int32), intent(in)   :: kp
    real(real64),   intent(in)   :: zp
    real(real64),   intent(in)   :: ap
    real(real64),   intent(in)   :: zt
    real(real64),   intent(in)   :: at
    real(real64),   intent(in)   :: elab
    real (real64)                :: kalbach

    integer(int32) :: jout
    real(real64)   :: a, ares, athrd, b, c, cut, ec, ecsq, ecut2, etest, flow, &
         & geom, p, ra, sig, signor, signor2, spill, w, xlamb, xmu, xnu, xnulam
    real(real64)   :: ecut

! ==============================================================================

    real(real64), parameter, dimension(6) :: xl0 = &
         & [ 12.10, 0.00437, 0.00619, 0.0186, 0.00459, 0.0643 ]
    real(real64), parameter, dimension(6) :: xl1 = &
         & [ -11.27,  -16.58,   -7.54,  -8.90,   -8.93, -13.96 ]
    real(real64), parameter, dimension(6) :: xm0 = &
         & [ 234.1,   244.7,   583.5,  686.3,   611.2,  781.2 ]
    real(real64), parameter, dimension(6) :: xm1 = &
         & [ 38.26,   0.503,   0.337,  0.325,    0.35,   0.29 ]
    real(real64), parameter, dimension(6) :: xn0 = &
         & [ 1.55,   273.1,   421.8,  368.9,   473.8, -304.7 ]
    real(real64), parameter, dimension(6) :: xn1 = &
         & [ -106.1,  -182.4,  -474.5, -522.2,  -468.2,  -470. ]
    real(real64), parameter, dimension(6) :: xn2 = &
         & [ 1280.8,  -1.872,  -3.592, -4.998,  -2.225, -8.580 ]
    real(real64), parameter, dimension(6) :: xp0 = &
         & [ -312.,   15.72,   0.798, -21.45,   -2.88,  10.95 ]
    real(real64), parameter, dimension(6) :: xp1 = &
         & [ 0.,    9.65,   420.3,  484.7,   205.6, -85.21 ]
    real(real64), parameter, dimension(6) :: xp2 = &
         & [ 0.,   -300.,  -1651., -1608.,  -1487.,  1146. ]

! ==============================================================================
    
    flow = 1.e-18
    spill = 1.e+18
    jout = nint(ap)
    ares = at+ap
    athrd =ares**0.3333
    signor = 1.
! signor reduces p and n result for light targs as per expt.
    if (kp == 1) then
       if (ares < 40.) signor=0.7+ares*0.0075
       xlamb = xl0(1)/athrd + xl1(1)
       xmu = xm0(1)*athrd + xm1(1)*athrd*athrd
       xnu = xn0(1)*athrd*ares + xn1(1)*athrd*athrd + xn2(1)
!       ec = 2.4
!       ecsq = 5.76
       ec = 0.5
       ecsq = 0.25
!       ec = 1.
!       ecsq = 1.
       p = xp0(1)
       xnulam = 1.
       etest = 32.
! etest is the energy above which the rxn cross section is
!    compared with the geometrical limit and the max taken.
!    xnulam here is a dummy value to be used later.
       ra = 0.
    else
       ra = 1.20
       if (kp == 2) then
          ra = 0.
          if (ares < 60.) then
             signor = 0.92
          else if (ares < 100.) then
             signor = 0.8 + ares*0.002
          end if
       end if
       ec = 1.44*zp*zt/(1.5*athrd+ra)
       ecsq = ec * ec
       p = xp0(kp) + xp1(kp)/ec + xp2(kp)/ecsq
       xlamb = xl0(kp)*ares + xl1(kp)
       a = ares**xm1(kp)
       xmu = xm0(kp) * a
       xnu = a* (xn0(kp)+xn1(kp)*ec+xn2(kp)*ecsq)
       if (jout == 2) ra=0.8
       if (jout == 3) ra=0.8
! new values of ra are for calculating the geometrical limit
!    to the cross section.
       if (kp == 2) then
          c = min(3.15_real64,ec*0.5)
          w = 0.7 * c / 3.15
! c and w are for the global corr'n factor for elab<ec
!    for light targs they are scaled down from global values
       end if
       xnulam = xnu / xlamb
       if (xnulam > spill) xnulam=0.
       if (xnulam < flow) go to 20
       if (kp == 2) then
          etest = sqrt(xnulam) + 7
       else
          etest = 1.2 * sqrt(xnulam)
       end if
       ! for xnulam > 0, sig reaches a maximum at sqrt(xnulam).
    end if
20  a = -2.*p*ec + xlamb - xnu/ecsq
    b = p*ecsq + xmu + 2.*xnu/ec
    ecut = 0.
    cut = a*a - 4.*p*b
    if (cut > 0.) ecut = sqrt(cut)
    ecut = (ecut-a) / (p+p)
    ecut2 = ecut

!    if (ecut < -0.05) then
!       c = -ecut * 0.5
!       w = -ecut * 0.1
!    else if (cut < 0) then
!       ecut2 = ecut * 0.25
!    end if

    if (cut < 0) ecut2 = ecut - 2.

!    sigmin = b - 0.25*a*a/p
! ecut is the energy where sigma=0 (if cut>0).  below ecut2
!    sigma is set identically to zero to avoid unphysical values.
    sig = 0.
    if (elab <= ec) then
       if (elab > ecut2) then
          sig = (p*elab*elab+a*elab+b) * signor
          if (kp == 2) then
             signor2 = (ec - elab - c) / w
             signor2 = 1 + exp(signor2)
             sig = sig / signor2

!             if (ecut < -0.05) then
!                if (elab < -ecut) then
!                   signor2 = (c - elab) / w
!                   signor2 = 1 + exp(signor2)
!                   sig = sig / signor2
!                end if
!             end if
          end if

! first signor gives empirical global corr'ns at low elab
!    second signor corrects values near elab=0; light nuclei
       end if
    else
       sig = (xlamb*elab+xmu+xnu/elab) * signor
       geom = 0.0_real64
       if (xnulam < flow) go to 36
       if (elab < etest) go to 36
       geom = sqrt(ap*elab)
       geom = 1.23*athrd + ra + 4.573/geom
       geom = 31.416 * geom * geom
       sig = max(geom,sig)
    end if
36  kalbach=sig
    return

  end function kalbach
! ==============================================================================

end module inverse_x_section

