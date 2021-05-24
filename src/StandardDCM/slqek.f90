
  subroutine slqek(l, ms, mb, ksi, me, lin, msin, mbin, mein, ln, &
       & msn, mbn, men)

! ======================================================================
!
!     Determine cross section type.
!
!     Called by: POINTE TYPINT
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February-March, 1996.
!   Edited by AJS, July, 1997.
!   Edited by A. J. Sierk, LANL T-16  October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    integer(int32), intent(  out) :: l
    integer(int32), intent(  out) :: ms
    integer(int32), intent(  out) :: mb
    integer(int32), intent(  out) :: ksi
    integer(int32), intent(  out) :: me
    integer(int32), intent(in   ) :: lin
    integer(int32), intent(in   ) :: msin
    integer(int32), intent(in   ) :: mbin
    integer(int32), intent(in   ) :: mein
    integer(int32), intent(in   ) :: ln
    integer(int32), intent(in   ) :: msn
    integer(int32), intent(in   ) :: mbn
    integer(int32), intent(in   ) :: men

! ======================================================================

!  ms, msin, & msn are always zero (strangeness; hardwired in cem95).
    ms = msin + msn
!  Total baryon number in system:
    mb = mbin + mbn
!  Total charge in system:
    me = mein + men
!  l, lin, & ln used for gamma reactions
    l = lin + ln
    if (l > 0) then
       ksi = 1
    elseif (mb > 1) then
!   baryon on baryon
       if (me.ne.1) then
!   2 neutrons or 2 protons
          ksi = 1
       else
!   1 neutron and 1 proton
          ksi = 2
       endif
    else
!   meson interacting with baryon
       if (me == 2 .or. me == -1) then
!   pi+ on proton or pi- on neutron
          ksi = 1
       elseif (me.ne.0) then
!   total charge = 1
          if (mein == 1) then
!   pi+ incident on neutron
             ksi = 2
          else
!   pi0 incident on proton
             ksi = 3
          endif
       else
!   total charge = 0
          if (mein == -1) then
!   pi- incident on proton
             ksi = 2
          else
!   pi0 incident on neutron
             ksi = 3
          endif
       endif
    endif
    return

! ======================================================================
! 1000 format ("Strange particle detected in 'slqek.f90'.")
! ======================================================================
  end subroutine slqek
