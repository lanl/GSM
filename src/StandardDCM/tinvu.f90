
  subroutine tinvu (sDCM, tin, cmin, tn, cmn, sintin, costin, sinfin, &
       & cosfin, sintn, costn, sinfn, cosfn, v, u, tin1)

! ======================================================================
!
!     Calculation of tin1, v, & u
!     tin1 is kinetic energy of the cascade particle in a frame in
!         which the nuclear partner particle is at rest.
!     u is the total energy available in the center-of-momentum frame.
!     v is the velocity of the center-of-momentum frame of the
!         interacting particles with respect to the nucleus.
!
!     The four variables ....in correspond to the direction of the
!     momentum of the cascade particle; those ....n to the direction
!     of the target particle (or two-particle state, for absorption)
!     in the lab.
!
!     Called by: ABSORP POINTE
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February-March, 1996.
!   Edited by AJS  July, 1997.
!   Edited by A. J. Sierk, LANL T-16  October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: one, two

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    real(real64), intent(in   ) :: tin      ! Kinetic and rest energies of particle 1
    real(real64), intent(in   ) :: cmin
    real(real64), intent(in   ) :: tn       ! Kinetic and rest energies of particle 2
    real(real64), intent(in   ) :: cmn
    real(real64), intent(in   ) :: sintin   ! Angles of one particle 1
    real(real64), intent(in   ) :: costin
    real(real64), intent(in   ) :: sinfin
    real(real64), intent(in   ) :: cosfin
    real(real64), intent(in   ) :: sintn    ! Angles of one particle 2
    real(real64), intent(in   ) :: costn
    real(real64), intent(in   ) :: sinfn
    real(real64), intent(in   ) :: cosfn
    real(real64), intent(  out) :: v(3)     ! Velocity in CP system of interacting particles
    real(real64), intent(  out) :: u        ! Energy available in CP system
    real(real64), intent(  out) :: tin1     ! Kinetic energy of cascade partner

    real(real64) :: ein, en, fac, pin, pinx, piny, pinz, pn, pnx, pny, &
         & pnz, temp, v2

! ======================================================================

!   Momentum and total energy of cascade particle in lab frame:
    pin = sqrt(abs(tin*(tin + two*cmin)))
    pinx = pin*sintin*cosfin
    piny = pin*sintin*sinfin
    pinz = pin*costin
    ein = tin + cmin
!   Momentum and total energy of Fermi sea particle (or of pair of
!   nucleons for pion absorption) in lab rframe:
    pn = sqrt(abs(tn*(tn + two*cmn)))
    pnx = pn*sintn*cosfn
    pny = pn*sintn*sinfn
    pnz = pn*costn
    en = tn + cmn
!   Total energy in lab frame.
    fac = ein + en
!   Velocity of the center-of-momentum system of the 2 (or 3)
!   particles in the lab frame.
    if (fac < div0Lim .and. fac > -div0Lim) then
       fac = div0Lim
       write(sDCM%io%message,1000) "65-67"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    v(1) = (pinx + pnx)/fac
    v(2) = (piny + pny)/fac
    v(3) = (pinz + pnz)/fac
    v2 = v(1)**2 + v(2)**2 + v(3)**2
!   Lorentz transformation of total energy into C.M. frame.
    u = fac*sqrt(abs(one - v2))
!   KINETIC energy of projectile in target rest frame
    temp = two*cmn
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message, 1000) "77"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    tin1 = (u**2 - (cmin + cmn)**2)/(temp)
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'tinvu.f90', line(s) ", A)
! ======================================================================
  end subroutine tinvu
