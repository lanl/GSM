
  subroutine sfint (sfname, flim1, flim2, bs, rs, fint, dataObj)

! ======================================================================
!
!    Gaussian quadrature; integral from flim2 to flim1 of fname.
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk  LANL  T-2  February, 1996.
!    "Last" change: 14-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMDataParams, only: zro, hlf, gaussWgt, gaussAbsc, numQuads

    implicit none
    procedure(IntegralInterface), intent(in   ), pointer :: sfname
    real(real64), intent(in   ) :: flim1
    real(real64), intent(in   ) :: flim2
    real(real64), intent(in   ) :: bs
    real(real64), intent(in   ) :: rs
    real(real64), intent(  out) :: fint
    class(StandardDCMData), intent(inout) :: dataObj

    integer(int32) :: it
    real(real64)   :: fname, fy, t1, t2, temp

! ======================================================================

    temp = zro
    t1 = hlf*(flim1 - flim2)
    t2 = hlf*(flim1 + flim2)
    do it = 1, numQuads
       fy = t1*gaussAbsc(it) + t2
       call sfname (fy, bs, rs, fname, dataObj)
       temp = temp + gaussWgt(it)*fname
    end do
    fint = temp*t1
    return

! ======================================================================
  end subroutine sfint


  subroutine sfint1 (sfname, flim1, flim2, bs, rs, fint1, dataObj)

! ======================================================================
!
!   Gaussian quadrature; integral from flim2 to flim1 of fname.
!   Identical to FINT; Needed to prevent a recursion when FINT (FIS...)
!   is called, since FIS uses FINT1.
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk  LANL  T-2  February-March, 1996.
!   "Last" change: 14-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMDataParams, only: zro, hlf, gaussWgt, gaussAbsc, numQuads

    implicit none
    procedure(IntegralInterface), intent(in   ), pointer :: sfname
    real(real64), intent(in   ) :: flim1
    real(real64), intent(in   ) :: flim2
    real(real64), intent(in   ) :: bs
    real(real64), intent(in   ) :: rs
    real(real64), intent(  out) :: fint1
    class(StandardDCMData), intent(inout) :: dataObj

    integer(int32) ::  it
    real(real64)   :: fname, fy, t1, t2, temp

! ======================================================================

    temp = zro
    t1 = hlf*(flim1 - flim2)
    t2 = hlf*(flim1 + flim2)
    do it = 1, numQuads
       fy = t1*gaussAbsc(it) + t2
       call sfname (fy, bs, rs, fname, dataObj)
       temp = temp + gaussWgt(it)*fname
    end do
    fint1 = temp*t1
    return

! ======================================================================
  end subroutine sfint1
