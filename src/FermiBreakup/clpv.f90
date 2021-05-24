
  subroutine clpv  (fbuObj, p1, v, p2, e1)

! ======================================================================
!
!    Called by: RASTAR
!
!   Transforms momentum p1 in frame moving with velocity v with respect
!   to the lab frame into the lab frame momentum p2.
!
!    Last change: 13-Aug-2003 BY NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fermiBreakupParams, only : zro, one

    implicit none
    class(FermiBreakup),        intent(inout) :: fbuObj
    real(real64), dimension(3), intent(in   ) :: p1
    real(real64), dimension(3), intent(in   ) :: v
    real(real64), dimension(3), intent(  out) :: p2
    real(real64),               intent(in   ) :: e1

    integer(int32) :: i, k
    real(real64)   :: gam, spv, temp, temp2, v2

! ======================================================================

    spv = zro
    v2 = zro
    do i = 1,3
       v2 = v2 + v(i)**2
       spv = spv + p1(i)*v(i)
    end do
    temp = one - v2
    if (temp < div0Lim) then
       temp = 0.01d0
       write(fbuObj%io%message,1050) "265"
       call fbuObj%io%print(4, 3, fbuObj%io%message)
    end if
    gam = one/sqrt(temp)
    temp = gam + one
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(fbuObj%io%message,1000) "271"
       call fbuObj%io%print(4, 3, fbuObj%io%message)
    end if
    temp2 = gam*(spv*gam/(temp) + e1)
    do k = 1,3
       p2(k) = p1(k) + v(k)*temp2
    end do
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'clpv.f90', line ", A)
1050 format("Square root/divide by zero error prevented in ", &
          & "'clpv.f90', line ", A)
! ======================================================================
  end subroutine clpv
