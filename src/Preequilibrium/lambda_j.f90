
module lambda_j
! Written by L. Kerby, 9/2014

  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private

  ! Parameters used in kin_energy, golden_section
  real(real64), public,  parameter :: error       =   0.2_real64               ! Maximum will be in an interval <= 2*error
  real(real64), public,  parameter :: goldenRatio =   1.6180339887_real64      ! Golden Ratio
  real(real64), public,  parameter :: invGolden   =   1.0_real64 / goldenRatio ! Inverse of Golden Ratio
  real(real64), private, parameter :: div0Lim     =   1.0d-15                  ! #/0 protection limit

  abstract interface
     function RANDOM() result(rang)
       use, intrinsic:: iso_fortran_env, only: real64
       implicit none
       real(real64) :: rang
     end function RANDOM
  end interface


  public  :: gamagu3
  public  :: kin_energy
  private :: g_laguerre_8
  private :: gu8
  private :: wbek
  private :: golden_section
  private :: rejection_gamma

contains

! ======================================================================

  subroutine gamagu3 (gamagu, j, int_energy, bind_energy, coul_bar, &
       & num_exciton, a_frag, z_frag, gam_beta, r0, ac, z_res_after, &
       & a_res_after, level_dens)

! ======================================================================
!
!   Finds the emission width (Gamma_j) for a preequilbrium particle of
!   type "j" by integrating lambda_j over excitation energy. (Eq. 3 of 
!   LA-UR-14-26657, by L. Kerby)
!
!   Called by: PRECOF
!
!   Calls: GU8, g_laguerre_8
!
!   CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk   LANL  T-2,  February, 1996.
!    Edited by AJS, August, 1997.
!    Modified by LMK, July 2014, to allow for any cross section.
!
! ======================================================================
 
    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    real(real64),   intent(  out) :: gamagu
    integer(int32), intent(in   ) :: j
    real(real64),   intent(in   ) :: int_energy
    real(real64),   intent(in   ) :: bind_energy
    real(real64),   intent(inout) :: coul_bar
    real(real64),   intent(in   ) :: num_exciton
    real(real64),   intent(in   ) :: a_frag
    real(real64),   intent(in   ) :: z_frag
    real(real64),   intent(in   ) :: gam_beta
    real(real64),   intent(in   ) :: r0
    real(real64),   intent(in   ) :: ac
    real(real64),   intent(in   ) :: z_res_after
    real(real64),   intent(in   ) :: a_res_after
    real(real64),   intent(in   ) :: level_dens

    real (real64)  :: a1, b1, y

! ======================================================================

    a1 = coul_bar
    b1 = int_energy - bind_energy
    if (b1 < div0Lim) then
       write(*,1000)
       gamagu = 0.
       return
    end if
    if ( num_exciton <= 15.) then
       ! Use an 8-point Gaussian Quadrature
       call gu8 (a1, b1, j, y, int_energy, bind_energy, coul_bar, &
            & num_exciton, a_frag, z_frag, gam_beta, r0, ac, z_res_after, &
            & a_res_after, level_dens)
    else
       ! Use an 8-point Gauss-Laguerre Quadrature
       call g_laguerre_8(j, int_energy, bind_energy, coul_bar, &
            & num_exciton, a_frag, z_frag, gam_beta, r0, ac, z_res_after, &
            & a_res_after, level_dens, y)
    end if
    gamagu = y
    return 

! ======================================================================
1000 format(3x, "Warning: b1 is zero in gamagu3 subroutine.")
! ======================================================================
  end subroutine gamagu3


  subroutine g_laguerre_8(part_type, int_energy, bind_energy, &
       & coul_bar, num_exciton, a_frag, z_frag, gam_beta, r0, ac, &
       & z_res_after, a_res_after, level_dens, gamma_j)

! ======================================================================
! 8-point Gauss-Laguerre integration of wbek (lambda_j)
! Written by L. Kerby, 9/2014
!    
!  Called by: gamagu3
!
!  Calls: wbek
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    integer(int32), intent(in   ) :: part_type
    real(real64),   intent(in   ) :: int_energy
    real(real64),   intent(in   ) :: bind_energy
    real(real64),   intent(inout) :: coul_bar
    real(real64),   intent(in   ) :: num_exciton
    real(real64),   intent(in   ) :: a_frag
    real(real64),   intent(in   ) :: z_frag
    real(real64),   intent(in   ) :: gam_beta
    real(real64),   intent(in   ) :: r0
    real(real64),   intent(in   ) :: ac
    real(real64),   intent(in   ) :: z_res_after
    real(real64),   intent(in   ) :: a_res_after
    real(real64),   intent(in   ) :: level_dens
    real(real64),   intent(  out) :: gamma_j     ! Probability of emitting fragment type j, integral of lambda_j over all emission energies, and "emission width"

    integer(int32) :: i
    real (real64)  :: increm
    real (real64)  :: e        ! epsilon, or T of emitted fragment, sampled/integrated over
    real (real64)  :: lambda_j ! Probability of emitting fragment type j with kinetic energy e, integrand of Gamma_j, and "partial transmission probability" 

! ======================================================================

    ! Abscissas for Gauss-Laguerre quadrature
    real(real64), parameter, dimension(8) :: x_i = &
         & [   0.1702796323051009997889, 0.903701776799379912186, &
         &     2.251086629866130689307,  4.266700170287658793649, &
         &     7.045905402393465697279, 10.75851601018099522406,  &
         &    15.74067864127800457803,  22.8631317368892641057   ]
    ! Weights for Gauss-Laguerre quadrature [Real weights * exp(x_i)]
    real(real64), parameter, dimension(8) :: w_i = &
         & [   0.4377234104929113732326, 1.0338693476655976425, &
         &     1.66970976565877574915, 2.376924701758599480959, &
         &     3.2085409133479262842, 4.2685755108251321986,    &
         &     5.8180833686719219283, 8.9062262152922114065   ]

! ======================================================================

    ! Setup variables
    gamma_j = 0.
    increm = int_energy / max( 1.5*(num_exciton-a_frag-1.) , &
         & (x_i(8)+1.)*int_energy/(int_energy-bind_energy-coul_bar) )

    ! Integrate via Gauss-Laguerre Quadrature
    do i=1,8
       e = coul_bar + increm*x_i(i)
       call wbek(part_type,e,lambda_j, int_energy, bind_energy, &
            & coul_bar, num_exciton, a_frag, z_frag, gam_beta, r0, &
            & ac, z_res_after, a_res_after, level_dens)

       ! Add portion of quadrature
       gamma_j = gamma_j + lambda_j*w_i(i)

    end do

    gamma_j = gamma_j * increm

    return
! ======================================================================
  end subroutine g_laguerre_8


  subroutine gu8 (a, b, j, y, int_energy, bind_energy, coul_bar, &
       & num_exciton, a_frag, z_frag, gam_beta, r0, ac, z_res_after, &
       & a_res_after, level_dens)

! ======================================================================
!
!    8 point Gaussian quadrature of function wbe.
!    y = Integral from a to b of wbe(j,e) de
!
!    Called by: gamagu3
!
!    Calls: wbek
!
!    CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk   LANL  T-2,  February, 1996.
!    Modified by L. Kerby, 8/2014.
!
! ======================================================================
 
    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    real(real64),   intent(in   ) :: a
    real(real64),   intent(in   ) :: b
    integer(int32), intent(in   ) :: j
    real(real64),   intent(  out) :: y
    real(real64),   intent(in   ) :: int_energy
    real(real64),   intent(in   ) :: bind_energy
    real(real64),   intent(inout) :: coul_bar
    real(real64),   intent(in   ) :: num_exciton
    real(real64),   intent(in   ) :: a_frag
    real(real64),   intent(in   ) :: z_frag
    real(real64),   intent(in   ) :: gam_beta
    real(real64),   intent(in   ) :: r0
    real(real64),   intent(in   ) :: ac
    real(real64),   intent(in   ) :: z_res_after
    real(real64),   intent(in   ) :: a_res_after
    real(real64),   intent(in   ) :: level_dens

    integer(int32) :: k
    real (real64)  :: d1, d2, e, wbe

! ======================================================================

    ! Weights for Gaussian quadrature [Real weights * exp(x_i)]
    real(real64), parameter, dimension(8) :: w = &
         & [   0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834, &
         &     0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363   ]
    ! Abscissas for Gaussian quadrature
    real(real64), parameter, dimension(8) :: fiks = &
         & [   0.9602898565, 0.7966664774, 0.5255324099, 0.1834346425, &
         &    -0.1834346425,-0.5255324099,-0.7966664774,-0.9602898565   ]

! ======================================================================

    ! Setup variables
    y = 0.0
    d1 = 0.5*(b - a)
    d2 = 0.5*(b + a)

    ! Integrate via Gauss-Laguerre Quadrature
    do k = 1,8
       e = d1*fiks(k) + d2
       call wbek(j, e, wbe,int_energy, bind_energy, coul_bar, &
            & num_exciton, a_frag, z_frag, gam_beta, r0, ac, &
            z_res_after, a_res_after, level_dens)

       ! Add portion of quadrature
       y = y + w(k)*wbe

    end do

    y = y*d1

    return
! ======================================================================
  end subroutine gu8


  subroutine wbek (j, e, wbe, int_energy, bind_energy, coul_bar, &
       & num_exciton, a_frag, z_frag, gam_beta, r0, ac, z_res_after, &
       & a_res_after, level_dens)

! ======================================================================
!
!   The integral of lambda_j (over KE of emitted fragment j) gives the 
!   preequilibrium particle emission width (Gamma_j) for particle of type "j"
!   lambda_j = Eq. 3 of LA-UR-14-26657, by L. Kerby
!
!   Called from: gu8, g_laguerre_8, kin_energy
!   Calls: nasa
!
!   CEM95 written by S. G. Mashnik
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Edited by AJS, August, 1997.
!   Modified by LMK, July 2014, to allow for any cross section.
!
! ======================================================================
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use preequilibriumParams, only: proton_mass
    use inverse_x_section, only: nasa

    implicit none
    integer(int32), intent(in   ) :: j
    real(real64),   intent(in   ) :: e   ! lambda_j, fragment energy T
    real(real64),   intent(  out) :: wbe ! lambda_j, fragment energy T
    real(real64),   intent(in   ) :: int_energy
    real(real64),   intent(in   ) :: bind_energy
    real(real64),   intent(inout) :: coul_bar
    real(real64),   intent(in   ) :: num_exciton
    real(real64),   intent(in   ) :: a_frag
    real(real64),   intent(in   ) :: z_frag
    real(real64),   intent(in   ) :: gam_beta
    real(real64),   intent(in   ) :: r0
    real(real64),   intent(in   ) :: ac
    real(real64),   intent(in   ) :: z_res_after
    real(real64),   intent(in   ) :: a_res_after
    real(real64),   intent(in   ) :: level_dens

    real(real64) :: x_sec    ! Inverse cross section (fm^2)
    real(real64) :: t_a_frag ! T of emitted fragment, per nucleon (MeV/A)
    real(real64) :: w        ! Total energy in center-of-momentum frame (MeV); LMK 07/2014
    real(real64) :: t_cm     ! Kinetic energy in center-of-momentum frame (MeV); LMK 07/2014
    real(real64) :: temp

! ======================================================================

    ! Ensure valid arguments were passed in
    if (int_energy < div0Lim .or. a_frag < div0Lim .or. r0 < div0Lim &
         & .or. a_res_after < div0Lim .or. ac < div0Lim) then
       ! ISSUE FOUND!
       ! Check where issue is
       if ( int_energy < div0Lim ) then
          ! Energy issue
          write(*,1000) "Energy issue in wbek, energy = ", int_energy
       else if ( a_frag < div0Lim ) then
          ! A issue
          write(*,1000) "A_frag issue in wbek, a_frag = ", a_frag
       else if ( r0 < div0Lim ) then
          ! r0 issue
          write(*,1000) "Radius issue in wbek, r0 = ", r0
       else if ( a_res_after < div0Lim ) then
          ! A issue 
          write(*,1000) "Residue (A) issue in wbek, a_res_after = ", a_res_after
       else if ( ac < div0Lim ) then
          ! ac issue
          write(*,1000) "ac issue in wbek, ac = ", ac
       end if

       wbe = 0.0
       write(*,2000) "324"
       return
    end if

! Used for testing Dostrovsky cross section
!  z = z_frag + z_res_after
!  a = a_frag + a_res_after
!      call orig_dost(j, z, a, e, z_frag, a_frag, coul_bar, x_sec)        ! Dostrovsky/Preeq Cross Section
! Transform from lab frame to center-of-momentum frame, LMK 07/2014
!********************************************************************************
! Relativistic Kinematics
!
! four-vector, momentum basis: (E, px*c, py*c, pz*c)
! W^2 = defined as = E^2 - p(dot)p*c^2    (this is invariant under transformations)
! W = total energy in center-of-momentum frame
!   = Tcm + ma*c^2 + mb*c^2
! 
! For a system of particles: (lab frame)
! W^2 = [sum(Ei)]^2 - [sum(Pi*c)]^2    (i for different particles in system)
!     = (Ea + Eb)^2 - (Pa*c)^2        (lab frame, Pb = 0)
!     = (Ea + Mb*c^2)^2 - (Pa*c)^2    (rest mass energy of target=Mb*c^2)
!     = (Ma*c^2)^2 + (Mb*c^2)^2 + 2*Ea*Mb*c^2
!    where Ea = sqrt((Pa*c)^2 + (Ma*c^2)^2) = Ta + ma*c^2
!
! mass of proton = 938 MeV/c^2
!********************************************************************************
    w = sqrt( (a_frag*proton_mass)**2 + (a_res_after*proton_mass)**2 + &
         & 2.*a_res_after*proton_mass*(e + proton_mass*a_frag) )
    t_cm = w - (proton_mass*a_frag + proton_mass*a_res_after)
    t_a_frag = e/a_frag

    call nasa(z_res_after, a_res_after, z_frag, a_frag, t_a_frag, t_cm, &
         & coul_bar, x_sec)   ! NASA Cross Section 

    x_sec = 0.1 * x_sec   ! convert from mb to fm^2

    if (j > 2) then
       if (e < -bind_energy) then
          wbe = 0.0
       else
          temp = sqrt(int_energy*a_frag) * int_energy * (r0**3) * &
               & a_res_after  !r_rms(nint(z+1),nint(a))**3) &
          wbe = 0.0330 * gam_beta*e*x_sec / temp &
               & * ( (e + bind_energy)/int_energy )**(a_frag-1.5) * &
               & (1 - (e+bind_energy)/int_energy)**( nint(num_exciton) &
               & - nint(a_frag)-1 ) * level_dens
       endif
    else
!  Protons and neutrons:
       wbe = 0.0000745 * level_dens/(ac*int_energy*a_res_after) * &
            & (1. - (bind_energy + e)/int_energy)**(nint(num_exciton)-2) &
            & * e * x_sec
    endif

    ! Ensure portion of integral is positive
    wbe = max(wbe, 0.0_real64)

    return
! ======================================================================
1000 format(3x, "Warning: ", A, f8.3)
2000 format(3x, "Warning: Divide by zero error found in 'lambda_j.f90', ", &
          & "line ", A)
! ======================================================================
  end subroutine wbek


  subroutine kin_energy(j, t_frag, int_energy, bind_energy, &
       & coul_bar, num_exciton, a_frag, z_frag, gam_beta, r0, ac, &
       & z_res_after, a_res_after, level_dens, rngPointer)

! ======================================================================
!
! Find the kinetic energy of emitted fragment j
! Written by L. Kerby, 09/2014

! Called by: precof
! Calls: cb_nasa, wbek, rngPointer
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use coulomb_barrier, only: cb_nasa

    implicit none
    integer(int32), intent(in   ) :: j
    real(real64),   intent(  out) :: t_frag      ! Kinetic energy of fragment (lab frame)
    real(real64),   intent(in   ) :: int_energy
    real(real64),   intent(in   ) :: bind_energy
    real(real64),   intent(inout) :: coul_bar
    real(real64),   intent(in   ) :: num_exciton
    real(real64),   intent(in   ) :: a_frag
    real(real64),   intent(in   ) :: z_frag
    real(real64),   intent(in   ) :: gam_beta
    real(real64),   intent(in   ) :: r0
    real(real64),   intent(in   ) :: ac
    real(real64),   intent(in   ) :: z_res_after
    real(real64),   intent(in   ) :: a_res_after
    real(real64),   intent(in   ) :: level_dens
    procedure(RANDOM), intent(in   ), pointer :: rngPointer

    integer(int32) :: iflag
    real(real64)   :: t_a, t_b            ! Minimum and maximum possible kinetic energy of fragment
    real(real64)   :: t_cm                ! Kinetic energy of fragment (cm frame)
    real(real64)   :: lambda_max, t_max   ! Maximum of lambda_j distribution and it's corresponding T
    real(real64)   :: x1, x2, x3, x4, y2, y3
    real(real64)   :: x, z, z_max, y_gam_max, y_gam, y_lamb

! ======================================================================

! Calculate coulomb barrier (NASA) without (emitted fragment) energy-dependence to find minimum t_a
    t_cm = 1.
    call cb_nasa(z_res_after, a_res_after, z_frag, a_frag, t_cm, coul_bar)

    t_a = coul_bar
    t_b = int_energy - bind_energy

! Finds the maximum (using the Golden Section method) of lambda_j on the interval (t_a, t_b)
    x1 = t_a
    x4 = t_b
    do while ((x4-x1) > error*2.)
       x2 = x4 - (x4-x1)*invGolden        
       x3 = x1 + (x4-x1)*invGolden
       call wbek(j, x2, y2, int_energy, bind_energy, coul_bar, &
            & num_exciton, a_frag, z_frag, gam_beta, r0, ac, &
            & z_res_after, a_res_after, level_dens)
       call wbek(j, x3, y3, int_energy, bind_energy, coul_bar, &
            & num_exciton, a_frag, z_frag, gam_beta, r0, ac, z_res_after, &
            & a_res_after, level_dens)  
       if (y2 > y3) then
          x4 = x3        ! Choose points 1, 2, 3
       else
          x1 = x2        ! Choose points 2, 3, 4
       end if
    end do
    t_max = (x1+x4)/2.
    call wbek(j, t_max, lambda_max, int_energy, bind_energy, &
         & coul_bar, num_exciton, a_frag, z_frag, gam_beta, r0, ac, &
         & z_res_after, a_res_after, level_dens)

! Samples the kinetic energy of the emitted fragment from the lambda_j distribution, using the Rejection Method with a Gamma Distribution (a=2) as the comparison function. 
! The Gamma Distribution has the form Gamma = lambda_max*e*x*exp(-x), where x = z/zmax and z = (T-Vj)/(uej-Bj-Vj)
! The maximum occurs at x=1. 
! Vj=t_a and (uej-Bj-Vj)=t_b-t_a, the minimum and range of the possible kinetic energy of the fragment.
! See NUMERICAL RECIPES, Section 7.3 Rejection Method: Gamma, Poisson, Binomial Deviates, p. 281.
    if (t_b-t_a < div0Lim) then
       iflag = 1
       t_frag = 0.01
       write(*,1000)
    else
       z_max = (t_max-t_a)/(t_b-t_a)
       iflag = 0
    end if
    do while (iflag == 0)
       ! Sample uniform point within Gamma distribution (t_frag, y_gam)
       x = -log(rngPointer()*rngPointer())
       y_gam_max = 1.01*lambda_max*2.7182818*x*exp(-x)   ! 1% cushion to ensure gamma > lambda
       y_gam = y_gam_max*rngPointer()
       z = x*z_max
       t_frag = z*(t_b-t_a) + t_a
       if (t_frag >= t_a .and. t_frag <= t_b) then
          ! Compare to lambda_j distribution
          call wbek(j, t_frag, y_lamb, int_energy, bind_energy, &
               & coul_bar, num_exciton, a_frag, z_frag, gam_beta, r0, &
               & ac, z_res_after, a_res_after, level_dens)
          if (y_gam <= y_lamb) then
             iflag = 1
          end if
       end if
    end do

    return
! ======================================================================
1000 format (3x, "Warning: Divide by zero error in 'lambda_j.f90', ", &
          & "subroutine 'kin_energy'.")
! ======================================================================
  end subroutine kin_energy



!Not used (in kin_energy)
  subroutine golden_section(j, x_a, x_b, ymax, xmax)

! ======================================================================
! Finds the maximum (using the Golden Section method) of lambda_j on the interval (x_a, x_b), LMK 08/2014
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    integer(int32), intent(in)  :: j
    real(real64),   intent(in)  :: x_a
    real(real64),   intent(in)  :: x_b
    real(real64),   intent(out) :: ymax ! (X,Y) point where maximum occurs
    real(real64),   intent(out) :: xmax ! (X,Y) point where maximum occurs

    real(real64) :: x1, x2, x3, x4, y2, y3

    ! Variables NOT passed in (i.e. with NO value assigned)
    real(real64) :: int_energy, bind_energy, coul_bar, num_exciton, &
         & a_frag, z_frag, gam_beta, r0, ac, z_res_after, a_res_after, &
         level_dens

! ======================================================================
    x1 = x_a
    x4 = x_b

    do while ((x4-x1) > error*2.)
       x2 = x4 - (x4-x1)*invGolden
       x3 = x1 + (x4-x1)*invGolden
       call wbek(j, x2, y2, int_energy, bind_energy, coul_bar, &
            & num_exciton, a_frag, z_frag, gam_beta, r0, ac, &
            & z_res_after, a_res_after, level_dens)
       call wbek(j, x3, y3, int_energy, bind_energy, coul_bar, &
            & num_exciton, a_frag, z_frag, gam_beta, r0, ac, &
            & z_res_after, a_res_after, level_dens)  
       if (y2 > y3) then
          x4 = x3   ! Choose points 1, 2, 3
       else
          x1 = x2   ! Choose points 2, 3, 4
       end if
    end do

    xmax = (x1+x4)/2.
    call wbek(j, xmax, ymax, int_energy, bind_energy, coul_bar, &
         & num_exciton, a_frag, z_frag, gam_beta, r0, ac, &
         & z_res_after, a_res_after, level_dens)


    return
! ======================================================================
  end subroutine golden_section


!Not used (in kin_energy)
  subroutine rejection_gamma(j, t_max, lambda_max, t_frag, t_a, t_b, rngPointer)

! ======================================================================
! Samples the kinetic energy of the emitted fragment from the lambda_j distribution, using the Rejection Method with a Gamma Distribution (a=2) as the comparison function. 
! The Gamma Distribution has the form Gamma = lambda_max*e*x*exp(-x), where x = z/zmax and z = (T-Vj)/(uej-Bj-Vj)
! The maximum occurs at x=1. 
! Vj=t_a and (uej-Bj-Vj)=t_b-t_a, the minimum and range of the possible kinetic energy of the fragment.
! See NUMERICAL RECIPES, Section 7.3 Rejection Method: Gamma, Poisson, Binomial Deviates, p. 281.
! LMK 09/2014
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    integer(int32), intent(in   ) :: j
    real(real64),   intent(in   ) :: t_max
    real(real64),   intent(in   ) :: lambda_max
    real(real64),   intent(  out) :: t_frag     ! Kinetic Energy of the emitted fragment
    real(real64),   intent(in   ) :: t_a
    real(real64),   intent(in   ) :: t_b
    procedure(RANDOM), intent(in   ), pointer :: rngPointer

    integer(int32) :: iflag
    real(real64)   :: x, z, z_max, y_gam_max, y_gam, y_lamb

    ! Variables NOT passed in (i.e. with NO value assigned)
    real(real64)   :: int_energy, bind_energy, coul_bar, num_exciton, &
         & a_frag, z_frag, gam_beta, r0, ac, z_res_after, a_res_after, &
         & level_dens

! ======================================================================
    if (t_b-t_a < div0Lim) then
       iflag = 1
       t_frag = 0.0
       write(*,1000)
    else
       z_max = (t_max-t_a)/(t_b-t_a)
       iflag = 0
    end if
    do while (iflag == 0)
       ! Sample uniform point within Gamma distribution (t_frag, y_gam)
       x = -log(rngPointer()*rngPointer())
       y_gam_max = 1.01*lambda_max*2.7182818*x*exp(-x)   ! 1% cushion to ensure gamma > lambda
       y_gam = y_gam_max*rngPointer()
       z = x*z_max
       t_frag = z*(t_b-t_a) + t_a
       if (t_frag >= t_a .and. t_frag <= t_b) then
          ! Compare to lambda_j distribution
          call wbek(j, t_frag, y_lamb, int_energy, bind_energy, &
               & coul_bar, num_exciton, a_frag, z_frag, gam_beta, r0, &
               & ac, z_res_after, a_res_after, level_dens)
          if (y_gam <= y_lamb) then
             iflag = 1
          end if
       end if
    end do

    return
! ======================================================================
1000 format (3x, "Warning: Divide by zero error in 'lambda_j.f90', ", &
          & "subroutine 'rejection_gamma'.")
! ======================================================================
  end subroutine rejection_gamma



! ======================================================================

end module lambda_j
