! Leslie M. Kerby, 08/2014, LANL

module coulomb_barrier
! Calculates Coulomb Barrier (in MeV) the following ways:
!  1) NASA [1, 2, 3]
!  1) Current/Old preequilibrium (Original Dostrovsky [4] with Vc = vk(l)*(1.439976/radncl)*zj(l)*zfj(l)/(ajthr(l) + afjthr(l))*(1 – u/(81*a*am)) )
!  2) Current GEM2 (Modified Dostrovsky[4] with Matsuse parameters)
!
! DEFAULT is ___NASA__ 
!
! [1] Tripathi et al., "Accurate universal parameterization of absorption cross sections," NIM B 117 (1996) 347.
! [2] Tripathi et al., "^... - neutron absorption cross sections," NIM B 129 (1997) 11.
! [3] Tripathi et al., "Universal Parameterization of Absorption Cross Sections: Light Systems," NASA/TP-1999-209726.
! [4] Dostrovsky et al., "Monte Carlo Calculations of Nuclear Evaporation Processes. III. Applications to Low-Energy Reactions",
!                        Physical Review, 116 (1959) 683.

  use fund_data, only: r_rms

  implicit none
  private 
! SUBROUTINE declarations
  public :: cb_nasa         ! Calculates Coulomb barrier by NASA model
  public :: cb_old_preeq    ! Calculates Coulomb barrier as in CEM03.03 preequilibrium 
  public :: cb_gem2         ! Calculates Coulomb barrier as in cCEM03.03 GEM2 evaporation (mod. Dostrovsky)

contains

! ==============================================================================

  subroutine cb_nasa(z_res_after, a_res_after, z_frag, a_frag, t_cm, coul_bar)

! ==============================================================================
! Calculates Coulomb Barrier for NASA parametrization
! Is energy dependent, calculates only one specfic fragment type at a time (not all 66) for a specific energy
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use numbers, only: thrd, ato3rd

    implicit none
    real (real64), intent(in) :: z_res_after   ! Z of residual nucleus AFTER emission of fragment
    real (real64), intent(in) :: a_res_after   ! A of residual nucleus AFTER emission of fragment
    real (real64), intent(in) :: z_frag        ! Z of emitted fragment
    real (real64), intent(in) :: a_frag        ! A of emitted fragment
    real (real64), intent(in) :: t_cm          ! Kinetic Energy in center of momentum frame (MeV)
    real (real64), intent(out) :: coul_bar     ! Coulomb barrier (MeV)

    integer(int32) :: ia_frag, ia_resAfter, iz_frag, iz_resAfter
    real (real64)  :: radius  ! Radius for coulomb barrier calculation

! ==============================================================================

    ia_frag     = nint(a_frag)
    iz_frag     = nint(z_frag + 1)
    ia_resAfter = nint(a_res_after)
    iz_resAfter = nint(z_res_after + 1)

    if ( iz_frag == 1 .or. iz_resAfter == 1) then  
! *********** neutrons ************
       coul_bar = 0.0
    else 
! *********** charged particles ***********
       ! For low-energy systems, Coulomb barrier is near zero
       if (t_cm < 1.0d-15) then
          ! write(*,1000) "52"
          coul_bar = 0.0_real64
          return
       end if

       ! Determine radius
       radius = 1.29*r_rms( iz_frag, ia_frag ) + &
            & 1.29*r_rms( iz_resAfter, ia_resAfter ) + &
            & 1.2*( ato3rd(ia_frag) + ato3rd(ia_resAfter) )/(t_cm)**(thrd)
       ! Calculate coulomb barrier
       coul_bar = 1.44 * z_frag * z_res_after/radius    ! (MeV)
    end if

    return
! ==============================================================================
1000 format(3x, "Warning: Divide by zero error in 'coulomb_barrier.f90', ", &
          & "line ", A)
! ==============================================================================
  end subroutine cb_nasa


  subroutine cb_old_preeq(z_res_before, a_res_before, t_frag, num_preeq_part, &
       z_type, a_type, coul_bar)

! ==============================================================================
! The Coulomb Barrier calculation as found in the old/current preequilibrium with 
! Coulomb Barrier Vc = vk(l)*(1.439976/radncl)*zj(l)*zfj(l)/(ajthr(l) + afjthr(l))*(1 – u/(81*a*am))   (MeV).
! Note that we simplify "am" by setting it equal to 0.1 instead of calculating it, 
! due to the complicated nature of the calculation of "am" (calls to fam, molnix, bf, and the calculation of the array egs).
! Values were also ballpark checked with the Coulomb Barrier calculated with reaction/fusion at nrv.jinr.ru/nrv/webnrv/qcalc/
!
! Linear interpolation and flat extrapolation used in the old preequilibrium code. Lagrange interpolation/extrapolation used here.
! Original Dostrovsky parameters:
!  z  kp  kalpha  
!  10  0.42  0.68
!  20  0.58  0.82
!  30  0.68  0.91
!  50  0.77  0.97
!  70  0.80  0.98
!  (100    1.00)  (used to have better extrapolation)
! Lagrange polynomial: 
! k(2)=  kp = -1.25e-8z^4 + 4.29e-6z^3 -5.2625e-4z^2 +0.02897z +0.17875
! k(6)=  kalpha = -5.291e-10z^5 +1.202e-7z^4 -7.754e-6z^3 -3.837e-5z^2 +1.8932e-2z +0.5011
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    real (real64),   intent(in)  :: z_res_before   ! Z of residual nucleus BEFORE emission of fragment
    real (real64),   intent(in)  :: a_res_before   ! A of residual nucleus BEFORE emission of fragment
    real (real64),   intent(in)  :: t_frag         ! Kinetic Energy of emitted fragment
    integer (int32), intent(in)  :: num_preeq_part ! Number of preequilibrium particles (1-66)
    real (real64),   intent(in)  :: a_type(num_preeq_part)     ! A for the particle types 1-66, used for emitted particle
    real (real64),   intent(in)  :: z_type(num_preeq_part)     ! Z for the particle types 1-66, used for emitted particle
    real (real64),   intent(out) :: coul_bar(num_preeq_part)   ! = kj*vj = Coulomb barrier (including parameter) [MeV]

    integer (int32) :: i
    real (real64)   :: z_res_after, a_res_after    ! Z and A of residual nucleus AFTER emission of fragment
    real (real64), dimension(66) :: kj = 1.0       ! Dostrovsky's original parametrization kj, 1=n, 2=p, 3=d...
    real (real64), dimension(66) :: vj = 1.0       ! Vj=(1.439976/radncl)*zj(l)*zfj(l)/(ajthr(l) + afjthr(l)) 

! ==============================================================================

! Neutron "coulomb barrier" = -beta in Original Dostrovsky
    ! (Setup vj, kj, coul_bar)
    vj(1) = 1.0
    kj(1) = -(2.12*(a_res_before**(-2./3.)) - 0.050)/(  0.76 + 2.2*( &
         & a_res_before**(-1./3.) )  ) !-beta
    coul_bar(1) = kj(1)*vj(1)

! Charged particle channel coulomb barrier
    kj(2) = -1.25e-8*(z_res_before**4) + 4.29e-6*(z_res_before**3) - &
         & 5.2625e-4*(z_res_before**2) + 0.02897*(z_res_before) + 0.17875
    kj(6) = -5.291e-10*(z_res_before**5) +1.202e-7*(z_res_before**4) - &
         & 7.754e-6*(z_res_before**3) - 3.837e-5*(z_res_before**2) + &
         & 1.8932e-2*(z_res_before) +0.5011
    kj(3) = kj(2) + 0.06
    kj(4) = kj(2) + 0.12
    kj(5) = kj(6) - 0.06
    kj(7:66) = 1.0
    do i=2,num_preeq_part  
       if (z_res_before < z_type(i) .or. a_res_before < a_type(i)) then 
          exit  
       end if
       z_res_after = z_res_before - z_type(i)
       a_res_after = a_res_before - a_type(i)
       vj(i) = (1.439976/1.5)*z_res_after*z_type(i)/(a_res_after**(1./3.) + &
            & a_type(i)**(1./3.))*(1. - t_frag/(81.*a_res_before*0.1))
       coul_bar(i) = kj(i)*vj(i)      ! MeV
    end do

    return
! ==============================================================================
  end subroutine cb_old_preeq


  subroutine cb_gem2(z_res_before, a_res_before, num_preeq_part, z_type, a_type, coul_bar)

! ==============================================================================
! The Coulomb Barrier calculation as found in GEM2 -> Vc = 1.43997*Zd*Zj/radius  (MeV)
! (1.43997 = h_bar*c/e*e, converts to MeV)
! 
! Linear interpolation and flat extrapolation used in the GEM2 code. Exponential or linear regression used here.
! Modified Dostrovsky parameters:
!  z  kp  kalpha  
!  20  0.51  0.81
!  30  0.60  0.85
!  40  0.66  0.89
!  50  0.68  0.93
!
! Exponential/linear regression: 
! k(2)=  k_p = 1.0 - 0.64215746*exp(-0.01477696*Z)
! k(6)=  k_alpha = 0.004*Z + 0.73
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    real (real64),   intent(in)  :: z_res_before   ! Z and A of residual nucleus BEFORE emission of fragment
    real (real64),   intent(in)  :: a_res_before   ! Z and A of residual nucleus BEFORE emission of fragment
    integer (int32), intent(in)  :: num_preeq_part ! Number of preequilibrium particles (1-66)
    real (real64),   intent(in)  :: z_type(num_preeq_part)     ! Z and A for the particle types 1-66, used for emitted particle, 1=n, 2=p, 3=d...
    real (real64),   intent(in)  :: a_type(num_preeq_part)     ! Z and A for the particle types 1-66, used for emitted particle, 1=n, 2=p, 3=d...
    real (real64),   intent(out) :: coul_bar(num_preeq_part)   != kj*vj = Coulomb barrier (including parameter)   (MeV)

    integer (int32)              :: i
    real (real64), dimension(num_preeq_part) :: z_res_after  ! Z of residual nucleus AFTER emission of fragment
    real (real64), dimension(num_preeq_part) :: a_res_after  ! A of residual nucleus AFTER emission of fragment
    real (real64), dimension(66) :: kj = 1.0           ! Dostrovsky's modified parametrization kj
    real (real64), dimension(66) :: vj = 1.0           ! Vj=1.43997*Zd*Zj/radius  (MeV) 
    real (real64), dimension(66) :: radius = 1.0       ! Radius => (r0*A^(1/3) + pj) for <= He4, (1.12*(A_d**(1./3.)-0.86*(A_d)**(-1./3.)) &
                                                       !    + (1.12*(A_f)**(1./3.)-0.86*(A_f)**(-1./3.)) + 3.75 for > He4 (fm)
    real (real64), dimension( 6) :: rj = 1.0           ! pj (or rj) in the radius calculation => 0 for n,p; 1.2 for d-He4

! ==============================================================================
    do i=1,num_preeq_part  
       if (z_res_before < z_type(i) .or. a_res_before < a_type(i)) then 
          exit  
       end if
       z_res_after(i) = z_res_before - z_type(i)
       a_res_after(i) = a_res_before - a_type(i)
    end do

    ! Neutron "coulomb barrier" = -beta in Modified Dostrovsky/GEM2
    vj(1) = 1.0
    kj(1) = -(1.66*(a_res_after(1)**(-2./3.)) - 0.050)/(  0.76 + 1.93 * &
         & ( a_res_after(1)**(-1./3.) )  ) !-beta
    coul_bar(1) = kj(1)*vj(1)

    ! Charged particle channel coulomb barrier
    kj(2) = 1.0 - 0.64215746*exp(-0.01477696*z_res_before)
    kj(6) = 0.004*z_res_before + 0.73
    kj(3) = kj(2) + 0.06
    kj(4) = kj(2) + 0.12
    kj(5) = kj(6) - 0.06
    kj(7:66) = 1.0
    rj(1:2) = 0.0
    rj(3:6) = 1.2
    radius(1:6) = 1.7*(a_res_after(1:6))**(1./3.) + rj(1:6) !fm
    radius(7:num_preeq_part) = (1.12*(a_res_after(7:num_preeq_part))**(1./3.)-0.86 * &
         & (a_res_after(7:num_preeq_part))**(-1./3.)) + (1.12*(a_type(7:num_preeq_part))**(1./3.) - &
         & 0.86 * ( a_type(7:num_preeq_part) )**(-1./3.) ) + 3.75    !fm
    vj(2:num_preeq_part) = (1.439976)*z_res_after(2:num_preeq_part)*z_type(2:num_preeq_part)/(radius(2:num_preeq_part))
    coul_bar(2:num_preeq_part) = kj(2:num_preeq_part)*vj(2:num_preeq_part) ! MeV
    return

  end subroutine cb_gem2
! ==============================================================================

end module coulomb_barrier
