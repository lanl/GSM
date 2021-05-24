
  function shellEnergy (molObj, a, z)

! ======================================================================
!
!    Ground state shell and pairing corrections vs. Z and A.
!
!    Called by: DELTAM PRECOF
!
!    Calls: MOLNIX
!
!    CEM95 written by S. G. Mashnik
!
!    Added ish to agrument list, so shell can be used either for ground
!    state shell and pairing correction to barrier heights, or for
!    shell correction to level densities.
!    Removed ish as an argument  AJS  (10/10/03)
!
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Modified by A. J. Sierk,  LANL  T-2  April-May, 1996.
!    "Last" change: 13-AUG-2003 by NVM
!    Modified by A. J. Sierk, LANL T-16  October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use molnixParams, only: zro

    implicit none
    class(Molnix), intent(inout) :: molObj
    real(real64),  intent(in   ) :: a
    real(real64),  intent(in   ) :: z
    real(real64)                 :: shellEnergy

    integer(int32) :: in, ina, iz
    real(real64)   :: un

! ======================================================================

    iz = nint(z)
    un = a - z
    in = nint(un)

!   Microscopic correction to ground state mass from the Finite Range
!   Liquid Drop Model of Moller, Nix, Myers, & Swiatecki, [Atomic Data
!   Nucl. Data Tables, 59, 185 (1995)].
    if (iz < 8)  then
       ina = in - molObj%nmina(iz) + 1_int32
       if (ina <= 0)  then 
          shellEnergy = zro
       else 
          shellEnergy = molObj%defineEnergy (iz, in, 1)
       endif
    else 
       shellEnergy = molObj%defineEnergy (iz, in, 1)
    endif
    return

! ======================================================================
  end function shellEnergy

