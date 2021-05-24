
  subroutine isobar (sDCM, u, v, tin1, partin, ipatne, mv, np, results)

! ======================================================================
!
!     Determining secondary particle characteristics for gamma-n
!     interaction with (3/2,3/2) isobar production.
!
!     See Barashenkov, et al. Nucl. Phys. A231 462 (1974).
!
!    CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk   LANL  T-2,  February, 1996.
!    Modified by SGM ????
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================
!
!  Definition of partin:
!                       partin(1); x coordinate of particle
!                       partin(2); y coordinate of particle
!                       partin(3); z coordinate of particle
!                       partin(4); sin(theta), direction of momentum
!                       partin(5); cos(theta), direction of momentum
!                       partin(6); sin(phi), direction of momentum
!                       partin(7); cos(phi), direction of momentum
!                       partin(8); kinetic energy of particle
!                       partin(9); rest mass of particle
!
!  Definition of ipatin:
!                       ipatin(1); charge of particle
!                       ipatin(2); (Always 0 in CEM95!)
!                       ipatin(3); strangeness of particle
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, hlf, one, two, massPiPM, massPiPM2, &
         & twpi, emneut, emprot
    use standardDCMData,   only: photonEG

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    real(real64),   intent(in   ) :: u
    real(real64),   intent(in   ) :: v(3)
    real(real64),   intent(in   ) :: tin1
    real(real64),   intent(in   ) :: partin(9)
    integer(int32), intent(in   ) :: ipatne(5)
    integer(int32), intent(in   ) :: mv
    integer(int32), intent(  out) :: np
    class(StandardDCMResults), intent(inout) :: results

    real(real64) :: a1, a2, a3, alpha, b1, bms, ctilpi, ctpi, edn, emnu, &
         & emnu2pi, emnupi, ent, epim, epit, f, f1, fipi, ftilpi, &
         & p, pim, pmt, r1, temp, temp1, temp2, tn, tpi, ts
    real(real64), dimension(3) :: vt=zro, ppim=zro, ppt=zro, ppit=zro, &
         & ppi=zro, pp=zro, pin=zro, pinst=zro, ppimst=zro, ppist=zro, &
         & ppst=zro

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
    results%pmemo(9,mv+1) = massPiPM
    results%pmemo(9,mv+2) = massPiPM
    if (ipatne(1) > 0) then
       results%imemo(1,mv+1) = -1
       results%imemo(1,mv+2) = 1
       results%imemo(1,mv+3) = 1
       results%pmemo(9,mv+3) = emprot
    else
       results%imemo(1,mv+1) = 1
       results%imemo(1,mv+2) = -1
       results%imemo(1,mv+3) = 0
       results%pmemo(9,mv+3) = emneut
    endif
    np = 3

    emnu = results%pmemo(9,mv+3)
    emnupi = emnu + massPiPM
    emnu2pi = emnupi + massPiPM
    temp = two*u
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "105"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    a1 = (u**2 + massPiPM2 - emnupi**2)/(temp)
    a2 = sqrt(abs(a1**2 - massPiPM2))
    a3 = u - a1
    temp = u
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "113"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    f1 = a1*a2*a3/temp
    alpha = 200.d0*f1
10  bms = sDCM%rang()*(u - emnu2pi) + emnupi
    temp = two*u
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "121"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    epim = (u**2 + massPiPM2 - bms**2)/(temp)
    edn = u - epim
    pim = sqrt(abs(epim**2 - massPiPM2))
    temp = u
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "129"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    f = (pim*epim*edn)/temp
    temp = emnu
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "135"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    ts = hlf*(bms**2 - emnupi**2)/temp
    temp = alpha
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "141"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    p = f*sDCM%qints(ts,9)/temp
    b1 = sDCM%rang()
    if (p > b1) then
       tpi = epim - massPiPM
    else
       go  to 10
    endif
    r1 = sDCM%rang()
    if (tin1 < one) then
       ctpi = sDCM%costa (28, tin1, r1)
    else
       ctpi = sDCM%costa (29, tin1, r1)
    endif
    fipi = twpi*sDCM%rang()
    temp = two*bms
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "160"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    epit = (bms**2 + massPiPM2 - emnu**2)/(temp)
    ent = bms - epit
    temp1 = sqrt(abs(one - ctpi**2))
    ppim(1) = pim*temp1*cos(fipi)
    ppim(2) = pim*temp1*sin(fipi)
    ppim(3) = pim*ctpi
    temp2 = tpi + massPiPM - u
    if (temp2 < div0Lim .and. temp2 > -div0Lim) then
       temp2 = div0Lim
       write(sDCM%io%message,1000) "171-173"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    vt(1) = ppim(1)/temp2
    vt(2) = ppim(2)/temp2
    vt(3) = ppim(3)/temp2
    ctilpi = one - two*sDCM%rang()
    ftilpi = twpi*sDCM%rang()
    temp1 = sqrt(abs(one - ctilpi**2))
    pmt = sqrt(abs(epit**2 - massPiPM2))
    ppit(1) = pmt*temp1*cos(ftilpi)
    ppit(2) = pmt*temp1*sin(ftilpi)
    ppit(3) = pmt*ctilpi
    tpi = epit - massPiPM
    vt(1) = -vt(1)
    vt(2) = -vt(2)
    vt(3) = -vt(3)
    call photonEG%lorentzTransform (ppit, vt, ppi, tpi, massPiPM)
    ppt(1) = -ppit(1)
    ppt(2) = -ppit(2)
    ppt(3) = -ppit(3)
    tn = ent - emnu
    call photonEG%lorentzTransform (ppt, vt, pp, tn, emnu)
    temp2  =  sqrt(abs(partin(8)*(partin(8) + two*partin(9))))
    pin(1) = temp2*partin(4)*partin(7)
    pin(2) = temp2*partin(4)*partin(6)
    pin(3) = temp2*partin(5)
    call photonEG%lorentzTransform (pin, v, pinst, partin(8), partin(9))
    call sDCM%rotor (pinst, v, ppim, ppimst)
    results%pmemo(1,mv+1) = ppimst(1)
    results%pmemo(2,mv+1) = ppimst(2)
    results%pmemo(3,mv+1) = ppimst(3)
    call sDCM%rotor (pinst, v, ppi, ppist)
    results%pmemo(1,mv+2) = ppist(1)
    results%pmemo(2,mv+2) = ppist(2)
    results%pmemo(3,mv+2) = ppist(3)
    call sDCM%rotor (pinst, v, pp, ppst)
    results%pmemo(1,mv+3) = ppst(1)
    results%pmemo(2,mv+3) = ppst(2)
    results%pmemo(3,mv+3) = ppst(3)

    return

! ======================================================================
1000 format("Divide by zero error in 'isobar.f90' line(s) ", A)
! ======================================================================
  end subroutine isobar
