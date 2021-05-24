
  subroutine abel (sDCM, partin, v, u, pist, pnst, cti, fii, cmi, cmn)

! ======================================================================
!
!     Calculates momenta of secondary particles
!     in center-of-mass system for absorption (called from ABSORP) and
!     elastic or pion charge-exchange scattering (called from ELEX).
!
!   Called by: ABSORP ELEX
!
!   Calls: ROTOR
!
!   CEM95 written by S. G. Mashnik
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Edited by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================
!
!  v is the velocity of the center-of-momentum frame in the lab frame.
!     pi + 2N system for absorption; 2N system for elastic scattering;
!     pi + N system for SCX scattering.
!  u is the TOTAL energy available in the center-of-momentum frame.
!  cti is random value of cos(theta); fii is random value of phi.
!
!  partin refers to the incident pion for pion absorption.
!
!  Definition of partin:
!                       partin(1); x coordinate of particle
!                       partin(2); y coordinate of particle
!                       partin(3); z coordinate of particle
!                       partin(4); sin(theta), direction of momentum
!                       partin(5); cos(theta), direction of momentum
!                       partin(6); sin(phi), direction of momentum
!                       partin(7); cos(phi), direction of momentum
!                       partin(8); (Lab) kinetic energy of particle
!                       partin(9); rest mass of particle
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: real64
    use standardDCMParams, only: zro, one, two

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    real(real64), intent(in   ) :: partin(9)
    real(real64), intent(in   ) :: v(3)
    real(real64), intent(in   ) :: u
    real(real64), intent(  out) :: pist(3)
    real(real64), intent(  out) :: pnst(3)
    real(real64), intent(in   ) :: cti
    real(real64), intent(in   ) :: fii
    real(real64), intent(in   ) :: cmi
    real(real64), intent(in   ) :: cmn

    real(real64) :: cfi=zro, ei=zro, fac=zro, fac2=zro, pim=zro, &
         & sfi=zro, spv=zro, sti=zro, temp=zro, temp1=zro, temp2=zro, &
         & temp3=zro, temp4=zro, v2=zro
    real(real64), dimension(3) :: pi = zro, pins = zro

! ======================================================================

! u = 0.0000000000001
    temp = two*u
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "56"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    ei = (u**2 + cmi**2 - cmn**2)/(temp)
!   |momentum| of each of the final state particles in the CM system:1
    temp = ei**2 - cmi**2
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "67"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    pim = sqrt(temp)
    v2 = v(1)**2 + v(2)**2 + v(3)**2
!   temp3 = cascade particle initial momentum in the lab frame.
    temp = partin(8)*(partin(8) + two*partin(9))
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "75"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    temp3 = sqrt(temp)
!   temp4 = transverse component of the cascade particle momentum in
!           the lab frame.
    temp4 = temp3*partin(4)

!   Dot product of the lab velocity of CM with the lab momentum of
!   cascade particle.
    spv = temp4*partin(7)*v(1) + temp3*partin(5)*v(3) + &
         & temp4*partin(6)*v(2)
!   fac = gamma (of CM frame with respect to lab!)
    temp = one - v2
    if (temp < div0Lim) then
       temp = 0.01d0
       write(sDCM%io%message,1200) "90"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    fac = one/sqrt(temp)
!   Temp1 = p (dot) v * (gamma - 1)/v**2
    temp = v2
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "97"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    temp1 = spv*(fac - one)/temp
!   Temp2 = Ein*gamma
    temp2 = (partin(8) + partin(9))*fac
    fac2 = temp1 - temp2
!   pins is the momentum vector of the particle transformed into the
!   CM system. (See Jackson, 2nd ed., p541 for Lor. trans. matrix)
    pins(1) = temp4*partin(7) + v(1)*fac2
    pins(2) = temp4*partin(6) + v(2)*fac2
    pins(3) = temp3*partin(5) + v(3)*fac2
    temp = one - cti**2
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "111"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    sti = sqrt(temp)
    cfi = cos(fii)
    sfi = sin(fii)
!   pins and v define a coordinate system in the CM frame,
!   which is rotated w.r.t. the lab coordinate system.
!   pi is (one) final state particle's momentum (with a random direction)
!   in this rotated CM frame; CM' (see ROTOR for definition of x' axes).
    pi(1) = pim*sti*cfi
    pi(2) = pim*sti*sfi
    pi(3) = pim*cti
!   Rotor finds the components of the final state particle's momentum in
!   the CM frame rotated to have its axes parallel to the corresponding
!   axes in the lab frame.
    call sDCM%rotor (pins, v, pi, pist)
    pnst(1) = -pist(1)
    pnst(2) = -pist(2)
    pnst(3) = -pist(3)
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'abel.f90', line(s) ", A)
1100 format("Square root error prevented in 'abel.f90', line(s) ", A)
1200 format("Divide by zero and square root error prevented in ", &
          & "'abel.f90', line(s) ", A)
! ======================================================================
  end subroutine abel
