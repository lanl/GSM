
  function qints (sDCM, x, lq)

! ======================================================================
!
!      Parabolic interpolation of tabulated cross sections.
!      For lq = 1-4 and 15-19, use a parabola in log(sigma) vs. log(E)
!      for energies (x) less than the 2nd point.
!
!    Called by: BINEL ISOBAR SIGMAT8
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Modified for new cross section tables, June, 1997 (AJS)
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Modified to define once an array of a, b, c coefficients
!      to speed up code.  AJS  February, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one
    use standardDCMData,   only: sigma, argus, aa, bb, cc, iiQints

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    real(real64),   intent(in   ) :: x
    integer(int32), intent(in   ) :: lq
    real(real64)                  :: qints

    integer(int32) :: i, lpha, lpmax, lqEff
    real(real64)   :: y
    logical        :: ilog

! ======================================================================

    qints = zro
    lpmax = 31
    ilog = .false.
    lqEff = lq
    if (lqEff <= 4) then
       ilog = .true.
       lpmax = 30
    elseif (lqEff > 14 .and. lqEff <= 19) then
       ilog = .true.
    elseif (lqEff == 20 .or. lqEff == 21 .or. lqEff == 28) then
       ilog = .true.
    elseif (lqEff >= 23 .and. lqEff <= 27) then
       ilog = .true.
    elseif (lqEff < 1 .or. lqEff > 28) then
       write(sDCM%io%message, 1000) lqEff
       call sDCM%io%print(2, 3, sDCM%io%message)

       ! Approximate to nearest valid point
       if ( lqEff > 28 ) lqEff = 28
       if ( lqEff <  1 ) lqEff =  1

       write(sDCM%io%message, 1010) lqEff
       call sDCM%io%print(2, 3, sDCM%io%message)

    endif
    i = iiQints(lqEff)
    lpha = 1

10  if (abs(x - argus(lpha,i)) <= 1.d-4) then
!  If x =~ table value, set cross section to tabulated one.
       qints = sigma(lpha,lqEff)
       return
    elseif (x > argus(lpha,i)) then
!  Increment lpha until x is <= tabulated x.
       if (lpha >= lpmax-1) then
          lpha = lpmax - 1
       else
          lpha = lpha + 1
          go to 10
       endif
    elseif (x < argus(lpha,i)) then
       if (lpha <= 1 .and. .not.ilog .and. i == 3) then
!  If x < first table value and not using parabolic log interpolation;
!  below threshold; set cross section to 0.
          qints = zro
          return
       elseif (lpha <= 1 .and. i == 2) then
!  For pi+N reactions with energy < 0; use 0-energy cross section:
          qints = sigma(1,lqEff)
          return
       else
          lpha = max(2,lpha)
          lpha = min(lpmax-1,lpha)
       endif
    endif
    y = x
    if (ilog .and. lpha == 2) y = log(x)
    qints = aa(lpha,lqEff)*y**2 + bb(lpha,lqEff)*y + cc(lpha,lqEff)
    if (ilog .and. lpha == 2) then
       qints = exp(qints)
    else
       qints = max(qints, zro)
    endif
    return

! ======================================================================
1000 format("An invalid 'lq' index (", i3, ") was detected in 'qints.f90'.")
1010 format("   Approximating lq to nearest valid point (", i3, ").")
! ======================================================================
  end function qints
