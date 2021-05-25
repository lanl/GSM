
  function wim (sDCM, clientTarg, iq, partin, ipatin)

! ======================================================================
!
!    Imaginary part of optical potential; protons for iq = 1, neutrons
!    for iq = 0.
!
!   Called by CASCAD
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February-March, 1996.
!   Edited by AJS, August, 1997.
!   Modified by AJS, March, 1999.
!
!   "Last" change: 13-AUG-2003 by NVM
!   Modified by A. J. Sierk, October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================
!
!  Definition of partin:
!                       partin(1); x coordinate of particle
!                       partin(2); y coordinate of particle
!                       partin(3); z coordinate of particle
!                       partin(8); kinetic energy of particle
!                       partin(9); rest mass of particle
!
!  Definition of ipatin:
!                       ipatin(1); charge of particle
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two, twthrd, twpi, emneut, &
         & emprot
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    class(StandardDCMData), intent(inout) :: clientTarg
    integer(int32),         intent(in   ) :: iq
    real(real64),           intent(in   ) :: partin(9)
    integer(int32),         intent(in   ) :: ipatin(5)
    real(real64)                          :: wim

    integer(int32) :: i1, i2
    real(real64)   :: a1a, am, am2, b, bs, ct, dz, et, fac, fctt, fi, &
         & p, p0, pp, ppx, ppy, ppz, r, r1, rho, s, sqf, st, t, temp, &
         & t0gev, tf, tp, tt1, u, vx, vy, vz, y, z

! ======================================================================

    real(real64), parameter :: a1 = 0.58_real64
    real(real64), parameter :: a2 = 0.51_real64

! ======================================================================

    i1 = ipatin(5)
    i2 = ipatin(1)
    am = partin(9)
    r = clientTarg%zoneBoundR( clientTarg%numZones() ) * &
         & sqrt(partin(1)**2 + partin(2)**2 + partin(3)**2)
    if (iq <= 0) then
!   Becchetti & Greenlees neutron imaginary potential
       r1 = 1.26d0*clientTarg%aTargThrd()
       temp = one + exp((r - r1)/a1)
       fac = one/(temp)
       tf = clientTarg%neutFermiMom(i1)
       rho = clientTarg%neutronDensity(1)*fac
       am2 = emneut
    else
!   Becchetti & Greenlees proton imaginary potential
       r1 = 1.32d0*clientTarg%aTargThrd()
       a1a = a2 + 0.7d0*(one - two * &
            & clientTarg%numProtons() / clientTarg%numBaryons() )
       if (a1a < div0Lim .and. a1a > -div0Lim) then
          a1a = div0Lim
          write(sDCM%io%message,1000) "93"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       temp = one + exp((r - r1)/a1a)
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message,1000) "98"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       fac = one/(temp)
       tf = clientTarg%protFermiMom(i1)
       rho = clientTarg%protonDensity(1)*fac
       am2 = emprot
    endif
    tp = tf*sDCM%rang()**twthrd
    t0gev = partin(8)
!   eps is the average particle separation energy (distance of Fermi
!   surface below 0 energy).
!   z is the asymptotic kinetic energy; y the momentum; pp the momentum
!   of the Fermi sea particle.
    z = abs(t0gev - tf - clientTarg%getSepEnergy() )
    y = sqrt(abs(z*(z + two*am)))
    pp = sqrt(abs(tp*(tp + two*am2)))
    ct = one - two*sDCM%rang()
    st = sqrt(abs(one - ct**2))
    fi = twpi*sDCM%rang()
!   Randomly selected Fermi motion of particle in nucleus:
    ppx = pp*st*cos(fi)
    ppy = pp*st*sin(fi)
    ppz = pp*ct
    p0 = sqrt(abs(t0gev*(t0gev + two*am)))
    tt1 = am + am2
    et = t0gev + tt1 + tp
    if (et < div0Lim .and. et > -div0Lim) then
       et = div0Lim
       write(sDCM%io%message,1000) "126-128"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    vx = ppx/et
    vy = ppy/et
    vz = ppz/et + p0/et
    u = et*sqrt(abs(one - vx**2 - vy**2 - vz**2))
    if (am < div0Lim .and. am > -div0Lim) then
       am = div0Lim
       write(sDCM%io%message,1000) "134"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    t = abs((u**2 - tt1**2)/(two*am))
    p = sqrt(abs(t*(t + two*am)))
!  b = velocity of projectile in target nucleon rest frame
    temp = t + am
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "142"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    b = p/(temp)
    if (b.ne.zro) then
       if (i2 == iq) then
          s = 10.63d0/b**2 - 29.92d0/b + 42.9d0
       else
          s = 34.10d0/b**2 - 82.20d0/b + 82.2d0
       endif
    else
       s = 1000.d0
    endif
    if (t0gev < div0Lim .and. t0gev > -div0Lim) then
       t0gev = div0Lim
       write(sDCM%io%message,1000) "156"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    bs = tf/t0gev
    dz = one - 1.4d0*bs
    if (bs > 0.5d0) then
       fctt = two - one/bs
       sqf = sqrt(fctt)
       dz = dz + 0.4d0*bs*sqf*fctt**2
    endif
    temp = z + am
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "168"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    wim = 10.d0*(y/(temp))*dz*rho*s

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'wim.f90' line(s) ", A)
! ======================================================================
  end function wim
