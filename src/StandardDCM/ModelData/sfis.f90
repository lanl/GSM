
  subroutine sfis (sxfi, bs, rs, fis, dataObj)

! ======================================================================
!
!   Calculates r**2 * (the Coulomb potential) as a function of r.
!
!   CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk  LANL  T-2  February, 1996.
!    Modified by A. J. Sierk  LANL  T-2  March, 1996.
!    Modified by AJS, March, 1999.
!   "Last" change: 14-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMDataParams, only: zro

    implicit none
    real(real64), intent(in   ) :: sxfi
    real(real64), intent(in   ) :: bs
    real(real64), intent(in   ) :: rs
    real(real64), intent(  out) :: fis
    type(StandardDCMData), intent(inout) :: dataObj

    real(real64)   :: ab1, ab2, rsm
    procedure(IntegralInterface), pointer :: ss1Ptr => ss1
    procedure(IntegralInterface), pointer :: ss2Ptr => ss2

! ======================================================================

    rsm = dataObj%target%zoneBoundR( dataObj%options%numZones )

    call sfint1 (ss2Ptr, sxfi, zro, bs, rs, ab1, dataObj)
    call sfint1 (ss1Ptr, rsm, sxfi, bs, rs, ab2, dataObj)
    fis = (ab1 + ab2*sxfi)*sxfi

    return

! ======================================================================
  end subroutine sfis

