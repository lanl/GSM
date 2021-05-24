
  function pairingGap (molObj, a, z)

! ======================================================================
!
!    Pairing gap energy shift
!
!    Calls: MOLNIX
!
!    CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Modified by A. J. Sierk,  LANL  T-2  April-May, 1996.
!    "Last" change: 13-AUG-2003 by NVM
!    Modified by A. J. Sierk, LANL T-16  October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    class(Molnix), intent(inout) :: molObj
    real(real64),  intent(in   ) :: a
    real(real64),  intent(in   ) :: z
    real(real64)                 :: pairingGap

    integer(int32) :: in, iz

! ======================================================================

    iz = nint(z)
    in = nint(a - z)

    pairingGap = molObj%defineEnergy (iz, in, 3)
    return

! ======================================================================
  end function pairingGap

