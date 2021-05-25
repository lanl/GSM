
  subroutine lorentzTransform (photonEG, p, v, pstar, t, cm)

! ======================================================================
!
!     Momentum calculation in system which has a relative velocity
!     v to given one.
!     Lorentz transformation; see Jackson, 2nd ed., p541.
!
!    Called by: ISOBAR STAT
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk  LANL  T-2  February, 1996.
!    Edited by AJS, August, 1997.
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================


    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(PhotonEventGenerator), intent(inout) :: photonEG
    real(real64),   intent(in   ) :: p(3)
    real(real64),   intent(in   ) :: v(3)
    real(real64),   intent(  out) :: pstar(3)
    real(real64),   intent(in   ) :: t
    real(real64),   intent(in   ) :: cm

    real(real64) :: spv, sv, temp1, temp2, v2

! ======================================================================

    v2 = v(1)**2 + v(2)**2 + v(3)**2
    temp1 = sqrt(abs(one - v2))
    spv = p(1)*v(1) + p(2)*v(2) + p(3)*v(3)
    if (temp1 < divZerLim ) then
       temp1 = divZerLim
       write(photonEG%io%message,1000) "45 and 46"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    temp2 = spv*(one/temp1 - one)/v2
    sv = (t + cm)/temp1

    pstar(1) = p(1) + v(1)*temp2 - v(1)*sv
    pstar(2) = p(2) + v(2)*temp2 - v(2)*sv
    pstar(3) = p(3) + v(3)*temp2 - v(3)*sv
    return

! ======================================================================
1000 format("Divide by zero error prevented in ", &
          & "'lorentzTransform.f90', line(s) ", A)
! ======================================================================
  end subroutine lorentzTransform
