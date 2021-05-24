
  subroutine kinema (sDCM, pstar, v, p, ct, st, cfi, sfi, t, cm)

! ======================================================================
!
!     Kinematics routine.
!     Lorentz transformation of the 3-momentum pstar in the center
!     of mass system to p in the laboratory system; t is the lab
!     kinetic energy.  v is the velocity of the CM in the lab system.
!
!    Called by: PAULIP PRECOF COALESL
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL T-2, February, 1996.
!    Edited by AJS, October, 1997.
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: real64
    use standardDCMParams, only : zro, one

    implicit none
    class(StandardDCM),         intent(inout) :: sDCM
    real(real64), dimension(3), intent(in   ) :: pstar
    real(real64), dimension(3), intent(in   ) :: v
    real(real64), dimension(3), intent(  out) :: p
    real(real64),               intent(  out) :: ct
    real(real64),               intent(  out) :: st
    real(real64),               intent(  out) :: cfi
    real(real64),               intent(  out) :: sfi
    real(real64),               intent(  out) :: t
    real(real64),               intent(in   ) :: cm

    real(real64) :: pm, pms, pmstars, spv, temp, temp1, temp2, temp3, &
         & temp4, v2

! ======================================================================

    spv = pstar(1)*v(1) + pstar(2)*v(2) + pstar(3)*v(3)
    v2 = v(1)**2 + v(2)**2 + v(3)**2
    temp1 = sqrt(abs(one - v2))
    if (temp1 < div0Lim .and. temp1 > -div0Lim) then
       temp1 = div0Lim
       write(sDCM%io%message,1000) "46, 49"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    temp = v2
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "46, 49"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    temp2 = spv*(one/temp1 - one)/temp
    pmstars = pstar(1)**2 + pstar(2)**2 + pstar(3)**2
    temp3 = sqrt(abs(pmstars + cm**2))
    temp1 = temp3/temp1
    temp3 = temp1 + temp2
    p(1) = pstar(1) + v(1)*temp3
    p(2) = pstar(2) + v(2)*temp3
    p(3) = pstar(3) + v(3)*temp3
    pms = p(1)**2 + p(2)**2 + p(3)**2
    pm = sqrt(abs(pms))
    temp = pm
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "61"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    ct = p(3)/temp
    temp4 = one - ct**2
    if (temp4 <= zro) then
       ct = one
       cfi = one
       sfi = zro
       st = zro
    else
       st = sqrt(abs(temp4))
       temp3 = pm*st
       if (temp3 < div0Lim .and. temp3 > -div0Lim) then
          temp3 = div0Lim
          write(sDCM%io%message,1000) "75, 76"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       cfi = p(1)/temp3
       sfi = p(2)/temp3
    endif
    t = sqrt(abs(pms + cm**2)) - cm

    return

! ======================================================================
1000 format("Divide by zero error prevented in ", &
          & "'kinema.f90', line(s) ", A)
! ======================================================================
  end subroutine kinema
