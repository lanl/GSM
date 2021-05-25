
  function subev (fbObj, u, e, f, n)

! ======================================================================
!
!    Quadratic interpolation of tabular values.
!
!    CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk   LANL  T-2,  February, 1996.
!    Edited by A. J. Sierk, LANL T-16  October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    class(FissionBarrier),        intent(inout) :: fbObj
    integer(int32),               intent(in   ) :: n
    real(real64),                 intent(in   ) :: u
    real(real64),   dimension(n), intent(in   ) :: e
    real(real64),   dimension(n), intent(in   ) :: f
    real(real64)                                :: subev

    integer(int32) :: j
    real(real64)   :: temp1, temp2, temp3, x1, x2, x3, y1, y2, y3

! ======================================================================

    if (u <= e(1)) then
       x1 = e(1)
       x2 = e(2)
       x3 = e(3)
       y1 = f(1)
       y2 = f(2)
       y3 = f(3)
    else
       if (u >= e(n-1)) then
          if (u > e(n)) then
             subev = f(n)
             return
          else
             x1 = e(n-2)
             x2 = e(n-1)
             x3 = e(n)
             y1 = f(n-2)
             y2 = f(n-1)
             y3 = f(n)
          endif
       else
          do j = 2,n
             if (u < e(j)) then
                x1 = e(j-1)
                x2 = e(j)
                x3 = e(j+1)
                y1 = f(j-1)
                y2 = f(j)
                y3 = f(j+1)
                go to 20
             endif
          end do
       endif
    endif
20  temp1 = (x1 - x2)*(x1 - x3)
    if (temp1 < div0Lim .and. temp1 > -div0Lim) then
       temp1 = div0Lim
       write(fbObj%io%message,1000) "75"
       call fbObj%io%print(4, 3, fbOBj%io%message)
    end if
    temp2 = (x2 - x1)*(x2 - x3)
    if (temp2 < div0Lim .and. temp2 > -div0Lim) then
       temp2 = div0Lim
       write(fbObj%io%message,1000) "76"
       call fbObj%io%print(4, 3, fbOBj%io%message)
    end if
    temp3 = (x3 - x1)*(x3 - x2)
    if (temp3 < div0Lim .and. temp3 > -div0Lim) then
       temp3 = div0Lim
       write(fbObj%io%message,1000) "77"
       call fbObj%io%print(4, 3, fbOBj%io%message)
    end if

    subev = y1*(((u - x2)*(u - x3))/(temp1)) + &
         y2*(((u - x1)*(u - x3))/(temp2)) + &
         y3*(((u - x1)*(u - x2))/(temp3))

    return
! ======================================================================
1000 format("Divide by zero error prevented in 'subev.f90', line(s) ", A)
! ======================================================================
  end function subev
