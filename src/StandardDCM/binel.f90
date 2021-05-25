
  subroutine binel (sDCM, partin, ipatin, ipatne, l, ms, mb, ksi, me, v, &
       & u, tin1, mv, np, nin, z, a, sigdif, sig1pi, r1, results)

! ======================================================================
!
!     Inelastic interaction (pion production!) subroutine.
!     For photon reactions, gam + N --> N + 2 pi.
!     (3 or more)-body final state.
!
!   Called by: TYPINT
!
!   Calls: CHINEL DIRECT8 ISOBAR STAT VMNSP
!
!     mv = last index filled in pmemo, imemo arrays.
!     np = number of particles in final state,
!     nin = flag indicating successful completion (0) or failure (2),
!     Photon interactions reinserted by AJS 10/21/03.
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL T-2, February, 1996.
!   Edited by AJS, August-September, 1997.
!   Modified by AJS, December, 1998.
!   Modified by A. J. Sierk, LANL T-16, October, 2003.
!   Moved RNDM call to calling routine, March, 2004 (AJS).
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================
!
!  Definition of partin (partne similar for 2nd particle):
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
!  Definition of ipatin (ipatne similar for 2nd particle):
!                       ipatin(1); charge of particle
!                       ipatin(2); Photon index
!                       ipatin(3); strangeness of particle
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two, massPiPM
    use standardDCMData,   only: pionProdThresh

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    real(real64),   intent(in   ) :: partin(9)
    integer(int32), intent(in   ) :: ipatin(5)
    integer(int32), intent(in   ) :: ipatne(5)
    integer(int32), intent(in   ) :: l
    integer(int32), intent(in   ) :: ms
    integer(int32), intent(in   ) :: mb
    integer(int32), intent(in   ) :: ksi
    integer(int32), intent(in   ) :: me
    real(real64),   intent(in   ) :: v(3)
    real(real64),   intent(in   ) :: u
    real(real64),   intent(in   ) :: tin1
    integer(int32), intent(in   ) :: mv
    integer(int32), intent(  out) :: np
    integer(int32), intent(  out) :: nin
    real(real64),   intent(in   ) :: z
    real(real64),   intent(in   ) :: a
    real(real64),   intent(in   ) :: sigdif
    real(real64),   intent(in   ) :: sig1pi
    real(real64),   intent(in   ) :: r1
    class(StandardDCMResults), intent(inout) :: results

    integer(int32) :: id, ik, ith, kp, lp
    real(real64)   :: betais, betath, temp, th

! ======================================================================

    nin = 0
    ik = 0
    if (ipatin(2) > 0) then
!  photon:
       betais = sDCM%qints (tin1, 26)/sDCM%qints (tin1, 27)
       if (r1 > betais) then
          call sDCM%statSDCM (u, v, partin, ipatne, mv, np, results)
       else
          call sDCM%isobar (u, v, tin1, partin, ipatne, mv, np, results)
       endif
       return
    endif
    id = ipatin(4) + 1
    if (tin1 <= pionProdThresh(id)) then
       betath = one
    else
       temp = sigdif
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message,1000) "116"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       betath = sig1pi/temp
    endif
    if (r1 < betath) then
!   One pion produced.
       ith = 1
       th = one
    else
!   N pions produced, where N is determined via Monte-Carlo in VMNSP.
       ith = 0
       th = two
    endif
    if ((u - partin(9) - massPiPM*th)  >  0.96d0) then
10     call sDCM%vmnsp (partin, ipatin, u, mv, np, ith, mb, tin1, lp, results)
       if (np <= 0) then
!   np = 0 means memory block has overflowed.
          nin = 2
          return
       endif
       if (lp == 0) then
          call sDCM%direct8 (v, tin1, mb, mv, np, partin, kp, ith, results)
          if (kp.ne.0) then
!   kp .ne. 0; abnormal return from DIRECT8.
             ik = ik + 1
             if (ik >= 50) then
                nin = 2
             else
                go to 10
             endif
          else
!   Find charges of outgoing particles:
             call sDCM%chinel (ipatin, l, ms, mb, ksi, np, mv, tin1, me, &
                  & ipatne, z, a, results)
          endif
       else
!  lp = 2; error from VMNSP:
          nin = 2
       endif
    else
       nin = 2
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'binel.f90', line(s) ", A)
! ======================================================================
  end subroutine binel
