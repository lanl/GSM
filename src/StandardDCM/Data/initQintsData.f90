
  function initializeQintsData () result(errorFlag)

! ======================================================================
!
! Sets up data used by the "qints" function (described below):
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

    implicit none
    integer(int32) :: errorFlag

    integer(int32) :: ij, lpha0, lpmax0, lq0
    real(real64) :: a, b, c, delta, deltaa, deltab, deltac, fact, phi1, &
         & phi2, phi3, psi1, psi2, psi3
    logical ilog1

! ======================================================================

    errorFlag = sDCMDataSetup

    do lq0 = 1,28

       ! Determine proper value for "ilog1"
       ilog1 = .false.
       lpmax0 = 31
       if (lq0 <= 4) then
          ilog1 = .true.
          lpmax0 = 30
       elseif (lq0 > 14 .and. lq0 <= 19) then
          ilog1 = .true.
       elseif (lq0 == 20 .or. lq0 == 21 .or. lq0 == 28) then
          ilog1 = .true.
       elseif (lq0 >= 23 .and. lq0 <= 27) then
          ilog1 = .true.
       endif


       ij = iiQints(lq0)
       do lpha0 = 2,lpmax0-1
          phi1 = sigma(lpha0-1,lq0)
          psi1 = argus(lpha0-1,ij)
          phi2 = sigma(lpha0,lq0)
          psi2 = argus(lpha0,ij)
          phi3 = sigma(lpha0+1,lq0)
          psi3 = argus(lpha0+1,ij)
          if (ilog1 .and. lpha0 == 2) then
             phi1 = log(sigma(lpha0-1,lq0))
             psi1 = log(max(argus(lpha0-1,ij), 1.d-6))
             phi2 = log(sigma(lpha0,lq0))
             psi2 = log(argus(lpha0,ij))
             phi3 = log(sigma(lpha0+1,lq0))
             psi3 = log(argus(lpha0+1,ij))
          endif
          a = psi2 - psi3
          b = psi3 - psi1
          c = psi1 - psi2
          delta = a*psi1**2 + b*psi2**2 + c*psi3**2
          deltaa = phi1*a + phi2*b + phi3*c
          deltab = (phi2 - phi3)*psi1**2 + (phi3 - phi1)*psi2**2 + &
               & (phi1 - phi2)*psi3**2
          deltac = (psi2*phi3 - psi3*phi2)*psi1**2 + &
               & (psi3*phi1 - psi1*phi3)*psi2**2 + &
               & (psi1*phi2 - psi2*phi1)*psi3**2
          fact = zro
          if (delta.ne.zro) fact = one/delta

          ! Initialize "qints" data (aa, bb, cc)
          aa(lpha0,lq0) = deltaa*fact
          bb(lpha0,lq0) = deltab*fact
          cc(lpha0,lq0) = deltac*fact
       end do
    end do

    return

! ======================================================================
  end function initializeQintsData
