
  function jtypb (ith, mb, lamb)

! ======================================================================
!
!     Determining type of coefficients b(n,k);
!     for angular distribution calculation.
!
!     ith = 0 or 1;
!     mb = baryon number of channel;
!     lamb = 1 for "recoiling" Fermi-sea nucleon,
!            2 for first pion produced,
!            3 for "recoiling" Projectile,
!          > 3 for subsequent pions for multiple pion production.
!
!    Called by: VMNSP
!
!    CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Edited by AJS, July, 1997.
!    Edited by AJS, December, 1998.
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32

    implicit none
    integer(int32), intent(in   ) :: ith
    integer(int32), intent(in   ) :: mb
    integer(int32), intent(in   ) :: lamb
    integer(int32)                :: jtypb

! ======================================================================

    if (ith.ne.0) then
!   One pion produced.
       if (mb <= 1) then
!   pi + N incident channel:
          if (lamb > 1) then
             jtypb = 6
          else
             jtypb = 5
          endif
       else
!   N + N incident channel:
          if (lamb <= 1 .or. lamb == 3) then
!  nucleon:
             jtypb = 1
          else
!  pion:
             jtypb = 2
          endif
       endif
    else
!   All other pion production cross sections:
       if (mb <= 1) then
!   pi + N incident channel:
          if (lamb <= 1) then
!  nucleon:
             jtypb = 7
          else
!  pion:
             jtypb = 8
          endif
       else
!   N + N incident channel:
          if (lamb <= 1 .or. lamb == 3) then
             jtypb = 3
          else
             jtypb = 4
          endif
       endif
    endif
    return

! ======================================================================
  end function jtypb
