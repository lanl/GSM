
  subroutine statSDCM (sDCM, u, v, partin, ipatne, mv, np, results)

! ======================================================================
!
!     Determining secondary particles characteristics for
!     gamma + N --> 2 pi interaction with statistical model.
!
!     See Barashenkov, et al. Nucl. Phys. A231 462 (1974).
!
!  Definition of partin:
!                       partin(4); sin(theta), direction of momentum
!                       partin(5); cos(theta), direction of momentum
!                       partin(6); sin(phi), direction of momentum
!                       partin(7); cos(phi), direction of momentum
!                       partin(8); kineti! energy of particle
!                       partin(9); rest mass of particle
!
!  Definition of ipatne:
!                       ipatne(1); charge of target nucleon
!
!    CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Edited by AJS,  August, 1997.
!    Modification: 22-oct-1998 by NVMokhov
!   "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two, thrd, twthrd, pi, twpi, &
         & emneut, emprot, massPiPM, massPi0
    use standardDCMData,   only: photonEG

    implicit none
    class(StandardDCM), intent(inout) :: sDCm
    real(real64),   intent(in   ) :: u
    real(real64),   intent(in   ) :: v(3)
    real(real64),   intent(in   ) :: partin(9)
    integer(int32), intent(in   ) :: ipatne(5)
    integer(int32), intent(in   ) :: mv
    integer(int32), intent(  out) :: np
    class(StandardDCMResults), intent(inout) :: results

    real(real64) :: aux, b1, ct1, ct2, ct3, e1, e2, emnu, emnupi, empi, &
         & empi2, epim, f1, fi1, fi2, fi3, p1, p2, p3, st1, st2, t1, &
         & t2, t3, temp, temp1, temp2, temp3, temp4, tpim
    real(real64), dimension(3) :: pv1=zro, pv2=zro, pv3=zro, ps1=zro, &
         & ps2=zro, ps3=zro, pin=zro, pinst=zro

! ======================================================================

    results%imemo(2,mv+1) = 0
    results%imemo(3,mv+1) = 0
    results%imemo(4,mv+1) = 0
    results%imemo(2,mv+2) = 0
    results%imemo(3,mv+2) = 0
    results%imemo(4,mv+2) = 0
    results%imemo(2,mv+3) = 0
    results%imemo(3,mv+3) = 0
    results%imemo(4,mv+3) = 1
!  Determine randomly the resulting particle types
    temp4 = sDCM%rang()
    if (ipatne(1) <= 0) then
!  I.  gamma + n case:
       results%imemo(1,mv+3) = 0
       results%pmemo(9,mv+3) = emneut
       if (temp4 <= thrd) then
!   n + 2 pi0
          results%imemo(1,mv+1) = 0
          results%pmemo(9,mv+1) = massPi0
          results%imemo(1,mv+2) = 0
          results%pmemo(9,mv+2) = massPi0
       elseif (temp4 < twthrd) then
!   p + pi0 + pi-
          results%imemo(1,mv+1) = 0
          results%pmemo(9,mv+1) = massPi0
          results%imemo(1,mv+2) = -1
          results%pmemo(9,mv+2) = massPiPM
          results%imemo(1,mv+3) = 1
          results%pmemo(9,mv+3) = emprot
       else
!   n + pi- + pi+
          results%imemo(1,mv+1) = -1
          results%pmemo(9,mv+1) = massPiPM
          results%imemo(1,mv+2) = 1
          results%pmemo(9,mv+2) = massPiPM
       endif
    else
!  II. gamma + p case:
       results%imemo(1,mv+3) = 1
       results%pmemo(9,mv+3) = emprot
       if (temp4 <= thrd) then
!   p + pi0 + pi0
          results%imemo(1,mv+1) = 0
          results%pmemo(9,mv+1) = massPi0
          results%imemo(1,mv+2) = 0
          results%pmemo(9,mv+2) = massPi0
       elseif (temp4 <= twthrd) then
!    n + pi+ + pi0
          results%imemo(1,mv+1) = 0
          results%pmemo(9,mv+1) = massPi0
          results%imemo(1,mv+2) = 1
          results%pmemo(9,mv+2) = massPiPM
          results%imemo(1,mv+3) = 0
          results%pmemo(9,mv+3) = emneut
       else
!    p + pi+ + pi-
          results%imemo(1,mv+1) = 1
          results%pmemo(9,mv+1) = massPiPM
          results%imemo(1,mv+2) = -1
          results%pmemo(9,mv+2) = massPiPM
       endif
    endif
    np = 3

!  Other combinations might be possible!?
!  HOW, if we restrict to N + 2 pi?  AJS  10/2/03
    empi = max(results%pmemo(9,mv+1), results%pmemo(9,mv+2))
    empi2 = empi**2
    emnu = results%pmemo(9,mv+3)
    emnupi = emnu + min(results%pmemo(9,mv+1), results%pmemo(9,mv+2))
    temp = u
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "130 and 137"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    epim = (temp**2 + empi2 - emnupi**2)/(two*temp)
    tpim = epim - empi
10  t1 = sDCM%rang()*tpim
    t2 = sDCM%rang()*tpim
    e1 = t1 + results%pmemo(9,mv+1)
    e2 = t2 + results%pmemo(9,mv+2)
    aux = u - e1 - e2
    f1 = 27.d0*e1*e2*aux/(temp**3)
    b1 = sDCM%rang()
    if (b1 < f1)then
       t3 = aux - results%pmemo(9,mv+3)
    else
       go to 10
    endif
    if (t3 <= zro) go to 10
    p1 = sqrt(abs(t1*(t1 + two*results%pmemo(9,mv+1))))
    p2 = sqrt(abs(t2*(t2 + two*results%pmemo(9,mv+2))))
    p3 = sqrt(abs(t3*(t3 + two*results%pmemo(9,mv+3))))
    temp1 = (p1 + p2 - p3)*(p1 - p2 + p3)*(p2 + p3 - p1)
    if (temp1 <= zro) go to 10
    ct3 = one - two*sDCM%rang()
    fi3 = twpi*sDCM%rang()
    temp2 = sqrt(abs(one - ct3**2))
    pv3(1) = p3*temp2*cos(fi3)
    pv3(2) = p3*temp2*sin(fi3)
    pv3(3) = p3*ct3
    temp3 = sqrt(abs(partin(8)*(partin(8) + two*partin(9))))
    pin(1) = temp3*partin(4)*partin(7)
    pin(2) = temp3*partin(4)*partin(6)
    pin(3) = temp3*partin(5)
    call photonEG%lorentzTransform (pin, v, pinst, partin(8), partin(9))
    call sDCM%rotor (pinst, v, pv3, ps3)
    temp = two*p3*p1
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "167"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    ct1 = (p2**2 - p1**2 - p3**2)/(temp)
    temp = two*p3*p2
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "173"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    ct2 = (p1**2 - p2**2 - p3**2)/(temp)
    fi1 = twpi*sDCM%rang()
    fi2 = pi + fi1
    st1 = sqrt(abs(one - ct1**2))
    st2 = sqrt(abs(one - ct2**2))
    pv1(1) = p1*st1*cos(fi1)
    pv1(2) = p1*st1*sin(fi1)
    pv1(3) = p1*ct1
    call sDCM%rotor (ps3, v, pv1, ps1)
    pv2(1) = p2*st2*cos(fi2)
    pv2(2) = p2*st2*sin(fi2)
    pv2(3) = p2*ct2
    call sDCM%rotor (ps3, v, pv2, ps2)
    results%pmemo(1,mv+1) = ps1(1)
    results%pmemo(2,mv+1) = ps1(2)
    results%pmemo(3,mv+1) = ps1(3)
    results%pmemo(1,mv+2) = ps2(1)
    results%pmemo(2,mv+2) = ps2(2)
    results%pmemo(3,mv+2) = ps2(3)
    results%pmemo(1,mv+3) = ps3(1)
    results%pmemo(2,mv+3) = ps3(2)
    results%pmemo(3,mv+3) = ps3(3)
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'stat.f90' line(s) ", A)
! ======================================================================
  end subroutine statSDCM
