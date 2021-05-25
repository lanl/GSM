
  function costa (sDCM, j, t, r1)

! ======================================================================
!
!     Cosine calculation for elastic and charge-exchange reactions.
!     For gamma + N, finds ang. distr. for N + 2 pi final state.
!
!   Called by: ABSORP COSEL COSEX DIRECT8 ISOBAR
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk  LANL  T-2  February, 1996.
!   Removed call to RNDM to calling program, AJS, March, 2004.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two
    use standardDCMData,   only: ankj

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    integer(int32), intent(in   ) :: j
    real(real64),   intent(in   ) :: t
    real(real64),   intent(in   ) :: r1
    real(real64)                  :: costa

    integer(int32) :: k, n
    real(real64)   :: cta, s1, s2, temp, temp1, term
    real(real64)   :: tk(4), r1n(5)

! ======================================================================

    s1 = zro
    r1n(1) = one
    r1n(2) = r1
    r1n(3) = r1*r1
    r1n(4) = r1*r1n(3)
    r1n(5) = r1n(3)*r1n(3)
    tk(1) = one
    tk(2) = t
    tk(3) = t*t
    tk(4) = t*t*t
    s2 = zro
    do n = 1,4
       do k = 1,4
          term = ankj(n,k,j)*tk(k)
          s1 = s1 + term*r1n(n)
          s2 = s2 + term
       end do
    end do
    temp = r1
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "54"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    cta = two*sqrt(temp)*(s1 + (one - s2)*r1n(5)) - one
    temp1 = abs(cta)
    if (temp1 <= one) then
       costa = cta
    else
       costa = sign(one, cta)
    endif
    return

! ======================================================================
1100 format("Square root error prevented in 'costa.f90', line(s) ", A)
! ======================================================================
  end function costa


  function costan (sDCM, j0, t0, r1)

! ======================================================================
!
!     Cosine calculation for elastic and charge-exchange reactions.
!     costan "n" for NEW; extrapolates old approximations from finite
!     angles different from 0 and 180 deg.  The old approximations
!     near 0 and 180 exhibit unphysical "pole-like" behavior.
!
!   Called by: COSEL COSEX
!
!   CEM95 written by S. G. Mashnik
!   COSTA Edited by A. J. Sierk  LANL  T-2  February, 1996.
!   Modified from COSTA by KK Gudima, Fall, 2003.
!   Edited by A. J. Sierk, LANL T-16, March, 2004.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two
    use standardDCMData, only: ankj

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    integer(int32), intent(in   ) :: j0
    real(real64),   intent(in   ) :: t0
    real(real64),   intent(in   ) :: r1
    real(real64)                  :: costan

    integer(int32) :: j, k, n
    real(real64)   :: a, b, c, cta, d, d12, d23, d31, da, db, dc, &
         & s1, s2, t, temp, temp1, term, x1, x2, x3, y1, y2, y3
    real(real64)   :: tk(4), r1n(5)

! ======================================================================

    real(real64), parameter :: &
         & rl1 = 0.20_real64, &
         & rl2 = 0.40_real64, &
         & rr1 = 0.90_real64, &
         & rr2 = 0.95_real64

 ! ======================================================================

    t = t0
    j = j0
    if (sDCM%options%newAngularDistributions > 0) then
       ! Use new version for angular distributions, if applicable


!    j = 12 to 16 never in calling list!!  AJS  (5/13/04)
!    Covered by the subroutine cosgamn.
       if (r1 < rl1 .or. r1 > rr2) then

          ! New angular distributions can be used
          if (r1 < rl1) then
             x1 = zro
             y1 = -one
             x2 = rl1
             y2 = sDCM%costa (j, t, x2)
             x3 = rl2
             y3 = sDCM%costa (j, t, x3)
          else
             x1 = one
             y1 = one
             x2 = rr2
             y2 = sDCM%costa (j, t, x2)
             x3 = rr1
             y3 = sDCM%costa (j, t, x3)
          endif
          d23 = x2 - x3
          d31 = x3 - x1
          d12 = x1 - x2
          d  = d23*x1*x1 + d31*x2*x2 + d12*x3*x3
          da = d23*y1    + d31*y2    + d12*y3
          db = (y2 - y3)*x1*x1 + (y3 - y1)*x2*x2 + (y1 - y2)*x3*x3
          dc = (x2*y3 - x3*y2)*x1*x1 + (x3*y1 - x1*y3)*x2*x2 + &
               & (x1*y2 - x2*y1)*x3*x3
          temp = d
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "137-139"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          a = da/temp
          b = db/temp
          c = dc/temp
          cta = a*r1*r1 + b*r1 + c
          go to 20
       endif

    end if

    !  Old version of angular distributions;
    !  Used for r1 >= rl1 and r1 <= rr2; or if sDCM%options%newAngularDistributions <= 0:
    r1n(1) = one
    r1n(2) = r1
    r1n(3) = r1*r1
    r1n(4) = r1*r1*r1
    r1n(5) = r1n(3)*r1n(3)
    tk(1) = one
    tk(2) = t
    tk(3) = t*t
    tk(4) = t*t*t
    s1 = zro
    s2 = zro
    do n = 1,4
       do k = 1,4
          term = ankj(n,k,j)*tk(k)
          s1 = s1 + term*r1n(n)
          s2 = s2 + term
       end do
    end do
    temp = r1
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "169"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    cta = two*sqrt(temp)*(s1 + (one - s2)*r1n(5)) - one


20  temp1 = abs(cta)
    if (temp1 <= one) then
       costan = cta
    else
       costan = sign(one, cta)
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'costa.f90', line(s) ", A)
1100 format("Square root error prevented in 'costa.f90', line(s) ", A)
! ======================================================================
  end function costan


  function cosex (sDCM, l, t, cm, r1, photoData)

! ======================================================================
!
!     Random cosine of scattering angle, weighted by the angular
!     distribution for charge-exchange scattering.
!     l is hard wired to be 0!!
!     l /=/ 0 statements from CEM95 reinserted by AJS (10/10/03)
!     Modified to use an extrapolated version of the angular
!     distributions for angles near 0 and pi.
!     Use improved angular distributions for gammas; KKG 2/04.
!
!   Called by: ELEX
!
!   Calls: COSGAMN COSTA COSTAN
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk  LANL  T-2  February, 1996.
!    Edited by AJS, July, 1997.
!    "Last" change: 12-AUG-2003 by N. V. Moknov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by K. K. Gudima (Ang. dist. approx.), Feb., 2004.
!    Added r1 to calling list to reduce number of routines using RNDM.
!      AJS, March, 2004.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: hlf, one, two, emnucg

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    integer(int32), intent(in   ) :: l
    real(real64),   intent(in   ) :: t
    real(real64),   intent(in   ) :: cm
    real(real64),   intent(in   ) :: r1
    type(sDCMPhotonCrossSections), intent(in   ) :: photoData
    real(real64)                  :: cosex

    real(real64) :: tb1, tem, temp, tmax

! ======================================================================

    if (l.ne.0) then
       if (t <= hlf) then
!  KKG 02/12/04:
          cosex = sDCM%cosgamn (14, r1, photoData)
       elseif (t <= one) then
          cosex = sDCM%cosgamn (15, r1, photoData)
       else
          cosex = sDCM%cosgamn (16, r1, photoData)
       endif
    else
       if (t <= 0.08d0) then
          cosex = sDCM%costan (17, t, r1)
       elseif (t <= 0.3d0) then
          cosex = sDCM%costan (18, t, r1)
       elseif (t <= one) then
          cosex = sDCM%costan (10, t, r1)
       elseif (t <= 2.4d0) then
          cosex = sDCM%costan (11, t, r1)
       else
!  KKG 02/16/04 (exact tmax):
          tem = two*emnucg
          temp = tem*t + (emnucg + cm)**2
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "254"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          tmax = tem*t*tem*(t + two*cm)/(temp)
          tb1 = 7.5d0*tmax
          temp = tb1
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "262"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          cosex = one + (two*log(one + r1 * &
               & (exp(-tb1) - one)))/temp
       endif
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'costa.f90', line(s) ", A)
! ======================================================================
  end function cosex


  function cosgamn (sDCM, j0, r1, photoData)

! ======================================================================
!
!     Cosine calculation for elastic and
!     charge-exchange gamma + N reactions using linear interpolation.
!     Energy of gamma is fixed in INIGAM.
!
!   Called by: COSEL COSEX
!
!   Written by K. K. Gudima, Nov. 2003
!   Edited by A. J. Sierk, LANL T-16, March, 2004.
!   Removed call to RNDM to calling subprogram, AJS, March, 2004.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, radianToDegree
    use standardDCMData,   only: thetai, cthetai

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    integer(int32), intent(in   ) :: j0
    real(real64),   intent(in   ) :: r1
    type(sDCMPhotonCrossSections), intent(in   ) :: photoData
    real(real64)                  :: cosgamn

    integer(int32) :: ir, ir1, ir2, jg
    real(real64)   :: cth, delx, dely, rrri, temp, th, x1, y1

! ======================================================================

    cosgamn = 0
    if (j0 == 12 .or. j0 == 13)  then
       jg = 2
    elseif (j0 == 14 .or. j0 == 15 .or. j0 == 16)  then
       jg = 1
    else
       ! Invalid value - warn client/user
       write(sDCM%io%message, 2000) j0
       call sDCM%io%print(2, 3, sDCM%io%message)
       write(sDCM%io%message, 2010) j0
       call sDCM%io%print(2, 3, sDCM%io%message)

       ! Approximate cross section:
       jg = 1
       if ( jg == 1 ) then
          write(sDCM%io%message, 2010)
       else if ( jg == 2 ) then
          write(sDCM%io%message, 2020)
       end if
       call sDCM%io%print(2, 3, sDCM%io%message)

    endif
    do ir = 1,181
       rrri = r1 - photoData%ri(jg,ir)
       if (abs(rrri) < 1.0d-5)  then
          cth = cthetai(ir)
          go to 20
       endif
    enddo
    do ir = 2,181
       rrri = r1 - photoData%ri(jg,ir)
       if (rrri < zro) then
          ir1 = ir - 1
          ir2 = ir
          go to 10
       endif
    enddo
    ir1 = 180
    ir2 = 181
10  continue
    x1 = photoData%ri(jg,ir1)
    delx = photoData%ri(jg,ir2) - x1
    y1 = thetai(ir1)
    dely = thetai(ir2) - y1
    temp = delx
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim 
       write(sDCM%io%message,1000) "337"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    th = y1 + dely*((r1 - x1)/temp)
    cth = cos(th*radianToDegree)
20  if (abs(cth) <= one) then
       cosgamn = cth
    else
       cosgamn = sign(one, cth)
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'costa.f90', line(s) ", A)
2000 format("Invalid 'j0' flag (", i5, ") in 'cosgamn' procedure.")
2010 format("   Approximating the cross section as: ", &
          & "gamma + p -> n + pi+")
2020 format("   Approximating the cross section as: ", &
          & "gamma + p -> p + pi0")
! ======================================================================
  end function cosgamn
