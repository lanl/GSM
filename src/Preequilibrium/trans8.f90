
  subroutine trans8 (preeqObj, p, h, am, c1, c2, c3, results)

! ======================================================================
!
!   Calculation of exciton transition rates; c1 for delta-n = +2,
!   c2 for delta-n = -2, and c3 for delta-n = 0.
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Modified by AJS, March, 1999.
!   Modified by SGM, 2000-2001, to get CEM2k
!
!   "Last" change: 13-AUG-2003 by NVM
!   Edited by A. J. Sierk, LANL T-16  October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use preequilibriumParams, only: hlf, zro, one, two, thr, four, fortyfive, &
         & emnucm

    implicit none
    class(Preequilibrium), intent(inout) :: preeqObj
    real(real64),          intent(in   ) :: p
    real(real64),          intent(in   ) :: h
    real(real64),          intent(in   ) :: am
    real(real64),          intent(  out) :: c1
    real(real64),          intent(  out) :: c2
    real(real64),          intent(  out) :: c3
    type(preequilibriumResults), intent(inout) :: results

    integer(int32) :: n
    real(real64)   :: aph, aph1, b1, b2, b3, en, est, f0, fct, fct0, &
         & fct2, ge, hm1, hm2, hm3, hp1, pn1, pp1, pp2, sf, spn, spp, &
         & svv, t, t1, temp

! ======================================================================

    t1 = p + h
    if (t1 < div0Lim .and. t1 > -div0Lim) then
       t1 = div0Lim
       write(preeqObj%io%message,1000) "42"
       call preeqObj%io%print(4, 3, preeqObj%io%message)
    end if
    est = 1.6d0*fortyfive + results%residual%kinEnergy/t1
    b2 = two*est/emnucm
    b1 = sqrt(abs(b2))
    if (b1 < div0Lim .and. b1 > -div0Lim) then
       b1 = div0Lim
       write(preeqObj%io%message,1000) "54, 74"
       call preeqObj%io%print(4, 3, preeqObj%io%message)
    end if
    if (b2 < div0Lim .and. b2 > -div0Lim) then
       b2 = div0Lim
       write(preeqObj%io%message,1000) "58"
       call preeqObj%io%print(4, 3, preeqObj%io%message)
    end if
    spp = 10.63d0/b2 - 29.92d0/b1 + 42.9d0
    spn = 34.10d0/b2 - 82.20d0/b1 + 82.2d0
    sf = (spp + spn)/two
    if (est < div0Lim .and. est > -div0Lim) then
       est = div0Lim
       write(preeqObj%io%message,1000) "64"
       call preeqObj%io%print(4, 3, preeqObj%io%message)
    end if
    b3 = fortyfive/est
    t = one - 1.4d0*b3
    if (b3 > hlf) then
       fct0 = two - one/b3      !No check needed because of if statement
       fct =  sqrt(fct0)*fct0**2
       t = t + 0.4d0*b3*fct
    endif
!   Average value of  (sigma * Vrel)  over excited states, and using
!   exclusion principle; divided by interaction volume. 
!   (lambda +) from Eq. 18 of Appendix 4 of CEM95 manual.
    temp = 1.2d0 + one/(4.7d0*b1)
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(preeqObj%io%message,1000) "79"
       call preeqObj%io%print(4, 3, preeqObj%io%message)
    end if
    svv = 0.00332d0*sf*t*sqrt(abs(est))/(temp)**3
    c1 = svv
    en = t1 + one
    n = nint(en)
    ge = am*results%residual%numBaryons*results%residual%kinEnergy
    pp1 = p + one
    hp1 = h + one
    pp2 = p + two
    hm1 = h - one
    hm2 = h - two
    hm3 = h - thr
    aph = (p*pp1 + h*hm3)/four
    aph1 = (pp1*pp2 + hp1*hm2)/four
    f0 = ge - aph
    temp = ge - aph1
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(preeqObj%io%message,1000) "98"
       call preeqObj%io%print(4, 3, preeqObj%io%message)
    end if
    pn1 = (f0/(temp))**n
    pn1 = min(1.d20, pn1)
    fct2 = c1*pn1*en

!   (lambda -) from Eq. 18 of Appendix 4 of CEM95 manual.
    if (f0 < div0Lim .and. f0 > -div0Lim) then
       f0 = div0Lim
       write(preeqObj%io%message,1000) "107"
       call preeqObj%io%print(4, 3, preeqObj%io%message)
    end if
    c2 = fct2*p*h*(t1 - two)/f0**2
    c2 = max (zro, c2)
    
!   (lambda 0) from Eq. 18 of Appendix 4 of CEM95 manual.
    temp = t1*f0
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(preeqObj%io%message,1000) "116"
       call preeqObj%io%print(4, 3, preeqObj%io%message)
    end if
    c3 = fct2*(p*(p - one) + four*p*h + h*hm1)/(temp)
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'trans.f90', line ", &
          & A, ".")
! ======================================================================
  end subroutine trans8
