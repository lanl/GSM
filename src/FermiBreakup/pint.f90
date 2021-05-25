
  subroutine pint (fbuObj, p, ct, st, cf, sf, t, cm)

! ======================================================================
!
!    Called by: RASTAR
!
!   Defining angles from the 3-momentum components.
!
!    Last change: 13-Aug-2003 BY NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: real64
    use fermiBreakupParams, only : zro, one

    implicit none
    class(FermiBreakup), intent(inout) :: fbuObj
    real(real64),        intent(in   ), dimension(3) :: p
    real(real64),        intent(  out) :: ct
    real(real64),        intent(  out) :: st
    real(real64),        intent(  out) :: cf
    real(real64),        intent(  out) :: sf
    real(real64),        intent(  out) :: t
    real(real64),        intent(in   ) :: cm

    real(real64) :: ctq, pmh, pq, pz, temp

! ======================================================================

    pz = p(3)**2
    pq = p(1)**2 + p(2)**2 + pz
    temp = pq
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(fbuObj%io%message,1000) "175"
       call fbuObj%io%print(4, 3, fbuObj%io%message)
    end if
    ctq = pz/temp
    temp = pq + cm**2
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(fbuObj%io%message,1000) "181"
       call fbuObj%io%print(4, 3, fbuObj%io%message)
    end if
    t = sqrt(temp) - cm
    if (ctq >= one) then
       st = zro
       ct = one
       sf = zro
       cf = one
    else
       temp = ctq
       if (temp < 0.0d0) then
          temp = 0.01d0
          write(fbuObj%io%message,1050) "193"
          call fbuObj%io%print(4, 3, fbuObj%io%message)
       end if
       ct = sqrt(temp)
       if (p(3) <= zro) ct = -ct
       st = sqrt(one - ctq)
       temp = pq
       if (temp < 0.0d0) then
          temp = 0.01d0
          write(fbuObj%io%message,1100) "201"
          call fbuObj%io%print(4, 3, fbuObj%io%message)
       end if
       pmh = st*sqrt(temp)
       temp = pmh
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(fbuObj%io%message,1000) "207, 208"
          call fbuObj%io%print(4, 3, fbuObj%io%message)
       end if
       cf = p(1)/temp
       sf = p(2)/temp
    endif

    return
! ======================================================================
1000 format("Divide by zero error prevented in 'pint.f90', line ", A)
1050 format("Square root/divide by zero error prevented in ", &
          & "'pint.f90', line ", A)
1100 format("Square root error prevented in 'pint.f90', line ", A)
! ======================================================================
  end subroutine pint
