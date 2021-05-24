
  subroutine drein1 (evapObj, j, s, a, eye1, eye0)

! ======================================================================
!
! ********This routine is the same as the drein1 in the LAHET code******
!
!     Compute statistical theory emission integrals.
!
!     for s<0.5 use a series expansion;
!     for s>0.5 the explicit relationship;
!
!     returns eye0 for neutrons only;
!     returns eye1 for all.
!
!    Correction of c1 by R. E. Prael
!    Extended precision of constants; AJS
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, hlf, one, thr

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    integer(int32),     intent(in   ) :: j
    real(real64),       intent(in   ) :: s
    real(real64),       intent(in   ) :: a
    real(real64),       intent(  out) :: eye1
    real(real64),       intent(  out) :: eye0

    integer(int32) :: n
    real(real64)   :: b, b0, c, exps, temp

! ======================================================================

! Coeficients for series expansions:
    real(real64), dimension(7), parameter :: c0 = &
         & [ 0.6666666666667d0, 0.2500000000000d0, 0.0666666666667d0, &
         &   0.0138888888889d0, 0.0023809523810d0, 0.0003472222222d0, &
         &   0.0000440917108d0   ]
    real(real64), dimension(7), parameter :: c1 = &
         & [ 0.5333333333333d0, 0.1666666666667d0, 0.0380952380952d0, &
         &   0.0069444444444d0, 0.0010582010582d0, 0.0001388888889d0, &
         &   0.0000160333494d0   ]

! ======================================================================

    exps = zro
    if (s < 1.d+2) exps = exp(-s)
    if (s >= hlf) then
!   Explicit relation:
       temp = a
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '608'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       b = hlf/temp
       eye1 = b*b*(thr + s*(s - thr) + exps*(hlf*s*s - thr))
       if (j == 1) eye0 = b*(s - one + exps)
       return
    else
!   Small s series expansion:
!  eye1=(1.0/(8*a*a))*(s**4/4)*(sum n=0 to 7:8*s**n/(n!*(n+2)*(n+4))
!   (*) exp(-s)!  AJS  10/16/03
       eye1 = one
       b = one
       do n = 1,7
          b = b*s
          c = b*c1(n)
          if (c < div0Lim) go to 10
          eye1 = eye1 + c
       end do
10     temp = a
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '631'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       b0 = hlf*hlf*s*s/temp
       eye1 = hlf*eye1*exps*b0*b0
       if (j == 1) then
! eye0 (neutrons only) = (.5/a)*s**2/2*(sum n=0 to 7:2*s**n/(n!*(n+2))
!   (*) exp(-s)!  AJS  10/16/03
          eye0 = one
          b = one
          do n = 1,7
             b = b*s
             c = b*c0(n)
             if (c < div0Lim) go to 20
             eye0 = eye0 + c
          end do
20        eye0 = exps*eye0*b0
       endif
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'drein.f90', line ", A)
! ======================================================================
  end subroutine drein1

  subroutine drein2 (evapObj, s, a, eye2)

! ======================================================================
!
! ********This routine is the same as the drein2 in LAHET code**********
!
!     Compute statistical theory emission integrals.
!
!     Compute third integral.
!
!     for s<0.5 use a series expansion.
!     for s>0.5 the explicit relationship
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Extended precision of constants; AJS
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, hlf, one, thr, thrd

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64),       intent(in   ) :: s
    real(real64),       intent(in   ) :: a
    real(real64),       intent(  out) :: eye2

    integer(int32) :: n
    real(real64)   :: b, c, exps, temp

! ======================================================================

!   Coeficients for series expansion:
    real(real64), dimension(7), parameter :: c2 = &
         & [ 0.4571428571428d0, 0.1250000000000d0, 0.0253968253978d0, &
         &   0.0041666666667d0, 0.0005772005772d0, 0.0000694444444d0, &
         &   0.0000074000074d0 ]

! ======================================================================

    exps = zro
    if (s < 1.d+02) exps = exp(-s)
    if (s >= hlf) then
!   Explicit relation:
       b = s*s
       temp = a*a*a
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '703'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       eye2 = 0.25d0*(s*(15.d0 - s*(6.d0 - s)) - 15.d0 + (15.d0 + &
            & 0.125d0*b*(b - 12.d0))*exps)/(temp)
    else
!   Series expansion:
!    eye2 = (1/(32a**3))*s**6/6*(sum n=0 to 7:48*s**n/(n!(n+2)(n+4)(n+6))
!   (*) exp(-s)!  AJS  10/16/03
       eye2 = one
       b = one
       do n = 1,7
          b = b*s
          c = b*c2(n)
          if (c < div0Lim) go to 10
          eye2 = eye2 + c
       end do
10     temp = a
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '723'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       b = hlf*hlf*s*s/temp
       eye2 = thrd*eye2*b*b*b*exps
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'drein.f90', line ", A)
! ======================================================================
  end subroutine drein2
