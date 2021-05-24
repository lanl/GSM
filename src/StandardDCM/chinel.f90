
  subroutine chinel (sDCM, ipatin, l, ms, mb, ksi, np, mv, tin1, me, &
       & ipatne, z, a, results)

! ======================================================================
!
!     Determinine secondary particles' charges in inelastic (pion prod.)
!     reaction.
!
!   Called by: BINEL
!
!   Calls: SIGMAT8
!
!   CEM95 written by S. G. Mashnik
!   Edited by A. J. Sierk  LANL T-2  February, 1996.
!   Edited by AJS, July-August, 1997.
!   "Last" change: 12-AUG-2003 by NVMokhov
!   Edited by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================
!
!  Definition of ipatin (ipatne similar for 2nd particle):
!                       ipatin(1); charge of particle
!                       ipatin(2); non-zero for photon interactions
!                       ipatin(3); strangeness of particle
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: thrd, twthrd, massPiPM, massPi0, &
         & emneut, emprot

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    integer(int32), intent(in   ) :: ipatin(5)
    integer(int32), intent(in   ) :: l
    integer(int32), intent(in   ) :: ms
    integer(int32), intent(in   ) :: mb
    integer(int32), intent(in   ) :: ksi
    integer(int32), intent(in   ) :: np
    integer(int32), intent(in   ) :: mv
    real(real64),   intent(in   ) :: tin1
    integer(int32), intent(in   ) :: me
    integer(int32), intent(in   ) :: ipatne(5)
    real(real64),   intent(in   ) :: z
    real(real64),   intent(in   ) :: a
    type(StandardDCMResults), intent(inout) :: results

    integer(int32) :: i, ic, itemp, lambda, &
         & migq, mtemp
    real(real64)   :: b1, b2, bpi0, bpiex, spi0, sth, temp, temp1, temp2

! ======================================================================

    if (np == 3) then
!   Three-particle final state:
!   spi0 = pi0 production cross section with no change of incident
!          particle identity;
!   pp --> p p pi0; pn --> p n pi0;nn --> n n pi0;pi+ p --> p pi+ pi0
!   pi- p --> p pi- pi0; pi0 p --> p 2pi0; pi0 n --> n 2pi0
       spi0 = sDCM%sigmat8 (l, ms, mb, ksi, 4, tin1)
!   Total one-pion production cross section:
       sth = sDCM%sigmat8 (l, ms, mb, ksi, 7, tin1)
       bpi0 = spi0/sth
!   Argument 5 --> pp --> p n pi+; pn --> p p pi-; pi- p --> n pi- pi+;
!   pi+ p --> n 2pi+; pi0 p --> p pi+ pi-; pi0 n --> n pi+ pi-
!   (Charged pion production cross section):
       bpiex = (spi0 + sDCM%sigmat8(l, ms, mb, ksi, 5, tin1))/sth
       temp1 = sDCM%rang()
       if (temp1 < bpi0) then
!   Incident cascade particle goes into mv + 3; nucleon from Fermi sea
!   goes into mv + 1; pi0 goes into mv + 2:
          results%imemo(1,mv+1) = ipatne(1)
          results%imemo(1,mv+2) = 0
          results%pmemo(9,mv+2) = massPi0
          results%imemo(1,mv+3) = ipatin(1)
       else
!   bpiex = 1.0 for mb,ksi = 1,1  or  2,1 pi+ + p; pi- + n; pp; or nn
!   bpiex < 1.0 for 1,3  or  2,2  pi0 + N or p + n
          if (temp1 >= bpiex) then
             results%imemo(1,mv+1) = ipatne(1) - (ipatin(4) - 1)*ipatin(1)
             results%imemo(1,mv+3) = (ipatin(4) - 1)*ipatin(1)**2 - &
                  & ipatin(4)*ipatin(1) + 1
             if (results%imemo(1,mv+3).ne.0 .and. results%pmemo(9,mv+3) < 0.15d0) &
                  & results%pmemo(9,mv+3) = massPiPM
          else
             results%imemo(1,mv+1) = 1 - ipatne(1)
             results%imemo(1,mv+3) = ipatin(1)
             if (results%imemo(1,mv+3).ne.0 .and. results%pmemo(9,mv+3) < 0.15d0) &
                  & results%pmemo(9,mv+3) = massPiPM
          endif
          results%imemo(1,mv+2) = me - results%imemo(1,mv+1) - results%imemo(1,mv+3)
          if (results%imemo(1,mv+2) == 0) results%pmemo(9,mv+2) = massPi0
       endif
       ic = results%imemo(1,mv+1)
       if (ic == 0) results%pmemo(9,mv+1) = emneut
       if (ic == 1) results%pmemo(9,mv+1) = emprot
       return
    endif
!   More than three particles in the final state:
10  continue
    b1 = sDCM%rang()
    temp = a
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message, 1000) "105, 115"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    if (b1 <= z/temp) then
       results%imemo(1,mv+1) = 1
    else
       results%imemo(1,mv+1) = 0
    endif
    ic = results%imemo(1,mv+1)
    if (ic == 0) results%pmemo(9,mv+1) = emneut
    if (ic == 1) results%pmemo(9,mv+1) = emprot
    if (mb >= 1) then
       b2 = sDCM%rang()
       if (b2 <= z/temp) then
          results%imemo(1,mv+3) = 1
       else
          results%imemo(1,mv+3) = 0
       endif
    endif
    lambda = 2
20  continue
    if (mb <= 1 .or. lambda.ne.3) then
       temp2 = sDCM%rang()
       mtemp = mv + lambda
       if (temp2 < thrd) then
          results%imemo(1,mtemp) = 1
          results%pmemo(9,mtemp) = massPiPM
       elseif (temp2 < twthrd) then
          results%imemo(1,mtemp) = 0
          results%pmemo(9,mtemp) = massPi0
       else
          results%imemo(1,mtemp) = -1
          results%pmemo(9,mtemp) = massPiPM
       endif
    endif
    if (lambda < np) then
       lambda = lambda + 1
       go to 20
    endif
    migq = 0
    do i = 1,np
       itemp = mv + i
       migq = migq + results%imemo(1,itemp)
    end do
    if (me.ne.migq) go to 10
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'chinel.f90', line(s) ", A)
! ======================================================================
  end subroutine chinel
