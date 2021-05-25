
  function jtypa (ith, mb, lamb)

! ======================================================================
!
!     Determining type of angular distribution coefficients a(n,k).
!     mb = baryon number.
!
!    Called by: DIRECT8
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32

    implicit none
    integer(int32), intent(in   ) ::  ith
    integer(int32), intent(in   ) ::  mb
    integer(int32), intent(in   ) ::  lamb
    integer(int32)                ::  jtypa

! ======================================================================

    if (ith.ne.0) then
!   iks = 7 cross section:
       if (mb <= 1) then
!   pi + N:
          if (lamb > 1) then
             jtypa = 25
          else
             jtypa = 24
          endif
       else
!   N + N
          if (lamb <= 1 .or. lamb == 3) then
             jtypa = 20
          else
             jtypa = 21
          endif
       endif
    else
!   iks = 4, 5, or 6 cross section:
       if (mb <= 1) then
!   pi + N:
          if (lamb <= 1) then
             jtypa = 26
          else
             jtypa = 27
          endif
       else
!   N + N
          if (lamb <= 1 .or. lamb == 3) then
             jtypa = 22
          else
             jtypa = 23
          endif
       endif
    endif
    return

! ======================================================================
  end function jtypa
