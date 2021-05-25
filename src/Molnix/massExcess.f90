
  function massExcess (molObj, a, z) result(deltaM)

! ======================================================================
!
!  Calculation of mass excess.                                         *
!
! ======================================================================
!
!   Called by: BINDNUC PRECOF PREQAUX RENORM
!
!   Calls: MOLNIX
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk  LANL  T-2  February, 1996.
!    Modified by AJS, May, 1996 (adding exptl. or Moller-Nix masses).
!    Modified by AJS, March, 1999
!   "Last" change: 12-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(Molnix), intent(inout) :: molObj
    real(real64),  intent(in   ) :: a
    real(real64),  intent(in   ) :: z
    real(real64)                 :: deltaM

    integer(int32) :: iz, n

! ======================================================================

    iz = nint(z)
    n  = nint(a - z)
    deltaM = molObj%defineEnergy(iz, n, 2)
    return

! ======================================================================
  end function massExcess
