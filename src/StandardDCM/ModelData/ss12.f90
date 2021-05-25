
  subroutine ss1 (sx, bs, rs, s1, dataObj)

! ======================================================================
!
!   r * Woods-Saxon function
!
!   CEM95 written by S. G. Mashnik
!   Edited by A. J. Sierk  LANL  T-2  February, 1996.
!   "Last" change: 14-AUG-2003 by NVMokhov
!   Edited by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    real(real64), intent(in   ) :: sx
    real(real64), intent(in   ) :: bs
    real(real64), intent(in   ) :: rs
    real(real64), intent(  out) :: s1
    type(StandardDCMData), intent(inout) :: dataObj

    real(real64) :: temp, temp1

! ======================================================================

    temp = bs
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(dataObj%io%message,1000) "26"
       call dataObj%io%print(4, 3, dataObj%io%message)
    end if
    temp1 = 1.d0 + exp((sx - rs)/temp)
    s1 = sx/(temp1)
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'ss12.f90' line(s) ", A)
! ======================================================================
  end subroutine ss1


  subroutine ss2 (sx, bs, rs, s2, dataObj)

! ======================================================================
!
!   r**2 * Woods-Saxon function.
!   (For integrating density radially.)
!
!   CEM95 written by S. G. Mashnik
!   Edited by A. J. Sierk  LANL  T-2  February, 1996.
!   Last change: 14-AUG-2003 by NVMokhov
!   Edited by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    real(real64), intent(in   ) :: sx
    real(real64), intent(in   ) :: bs
    real(real64), intent(in   ) :: rs
    real(real64), intent(  out) :: s2
    type(StandardDCMData), intent(inout) :: dataObj

    real(real64) :: temp, temp1

! ======================================================================

    temp = bs
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(dataObj%io%message,1000) "59"
       call dataObj%io%print(4, 3, dataObj%io%message)
    end if
    temp1 = 1.d0 + exp((sx - rs)/temp)
    s2 = sx**2/(temp1)
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'ss12.f90' line(s) ", A)
! ======================================================================
  end subroutine ss2
