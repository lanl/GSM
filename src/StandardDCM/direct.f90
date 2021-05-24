
  subroutine direct8 (sDCM, v, tin1, mb, mv, np, partin, kp, ith, results)

! ======================================================================
!
!     Determines direction of secondary particle's motion.
!
!     kp = 0 denotes a  satisfactory termination.
!     kp = 2 denotes an unsatisfactory termination.
!
!   Called by: BINEL
!
!   Calls: COSTA JTYPEA ROTOR
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk  LANL  T-2  February, 1996.
!    Edited by AJS, July, 1997.
!
!   "Last" change: 12-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================
!
!  Definition of partin (pmemo is similar):
!                       partin(1); x coordinate of particle
!                       partin(2); y coordinate of particle
!                       partin(3); z coordinate of particle
!                       partin(4); sin(theta), direction of momentum
!                       partin(5); cos(theta), direction of momentum
!                       partin(6); sin(phi), direction of momentum
!                       partin(7); cos(phi), direction of momentum
!                       partin(8); kinetic energy of particle
!                       partin(9); rest mass of particle
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, hlf, one, two, pi, twpi

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    real(real64),   intent(in   ) :: v(3)
    real(real64),   intent(in   ) :: tin1
    integer(int32), intent(in   ) :: mb
    integer(int32), intent(in   ) :: mv
    integer(int32), intent(in   ) :: np
    real(real64),   intent(in   ) :: partin(9)
    integer(int32), intent(  out) :: kp
    integer(int32), intent(in   ) :: ith
    class(StandardDCMResults), intent(inout) :: results

    integer(int32) :: ja, lambda, ltemp, m1, m2, m1temp, m2temp, nd
    real(real64)   :: cfm1, cfm2, ctl, ctm1, ctm2, fl, fm1, fm2, pakvm, &
         & r1, s2, sfm1, sfm2, spv, st, stl, stm1, stm2, temp, temp1, temp2, &
         & temp3, temp4, temp5, temp6, v2
    real(real64), dimension(3) :: pakv=zro, plst=zro, pl=zro, pakst=zro, &
         & pin=zro

! ======================================================================

    nd = 0
    kp = 0
    st = sDCM%rang()
    if (mb > 1) then
!   Two baryons:
       if (st < hlf) then
          m1 = 1
          m2 = 3
       else
          s2 = sDCM%rang()
          if (s2 < hlf) then
             m1 = 2
             m2 = 3
          else
             m1 = 1
             m2 = 2
          endif
       endif
    else
!   One baryon:
       if (st < hlf) then
          m1 = 2
          m2 = 3
       else
          s2 = sDCM%rang()
          if (s2 < hlf) then
             m1 = 1
             m2 = 3
          else
             m1 = 1
             m2 = 2
          endif
       endif
    endif
20  lambda = 1
    pakv(1) = zro
    pakv(2) = zro
    pakv(3) = zro
    m1temp = mv + m1
    m2temp = mv + m2
30  ltemp = mv + lambda
    if (lambda.ne.m1 .and. lambda.ne.m2) then
       ja = jtypa (ith, mb, lambda)
       r1 = sDCM%rang()
       ctl = sDCM%costa (ja, tin1, r1)
       fl = twpi*sDCM%rang()
       temp = one - ctl**2
       if (temp < 0.0d0) then
          temp = 0.01d0
          write(sDCM%io%message,1100) "110"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       stl = sqrt(temp)
       temp1 = cos(fl)
       temp2 = sin(fl)
! pmemo(4,  ) temporarily contains momentum magnitude:
       temp3 = results%pmemo(4, ltemp)
       results%pmemo(1,ltemp) = temp3*stl*temp1
       results%pmemo(2,ltemp) = temp3*stl*temp2
       results%pmemo(3,ltemp) = temp3*ctl
       pakv(1) = pakv(1) + results%pmemo(1,ltemp)
       pakv(2) = pakv(2) + results%pmemo(2,ltemp)
       pakv(3) = pakv(3) + results%pmemo(3,ltemp)
    endif
    if (lambda < np) then
       lambda = lambda + 1
       go to 30
    endif
    pakvm = sqrt(pakv(1)**2 + pakv(2)**2 + pakv(3)**2)
    if (np.ne.3) then
       if (results%pmemo(4,m1temp) >= (pakvm + results%pmemo(4,m2temp))) go to 50
       if (results%pmemo(4,m1temp) <= abs(pakvm - results%pmemo(4,m2temp))) go to 50
    endif
    v2 = v(1)**2 + v(2)**2 + v(3)**2
    temp = partin(8)*(partin(8) + two*partin(9))
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "137"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    temp4 = sqrt(temp)
    spv = temp4*partin(4)*partin(7)*v(1) + temp4*partin(5)*v(3) + &
         & temp4*partin(4)*partin(6)*v(2)
    temp = one - v2
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "145"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    temp6 = sqrt(temp)
    temp = v2
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "150"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    if (temp6 < div0Lim .and. temp6 > -div0Lim) then
       temp6 = div0Lim
       write(sDCM%io%message,1000) "155"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    temp5 = spv*(one/temp6 - one)/temp
    temp = temp6
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "161"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    temp6 = (partin(8) + partin(9))/temp
    temp3 = temp5 - temp6
    pin(1) = temp4*partin(4)*partin(7) + v(1)*temp3
    pin(2) = temp4*partin(4)*partin(6) + v(2)*temp3
    pin(3) = temp4*partin(5) + v(3)*temp3
    lambda = 1
40  ltemp = mv + lambda
    if (lambda.ne.m1 .and. lambda.ne.m2) then
       pl(1) = results%pmemo(1,ltemp)
       pl(2) = results%pmemo(2,ltemp)
       pl(3) = results%pmemo(3,ltemp)
       call sDCM%rotor (pin, v, pl, plst)
       results%pmemo(1,ltemp) = plst(1)
       results%pmemo(2,ltemp) = plst(2)
       results%pmemo(3,ltemp) = plst(3)
    endif
    if (lambda < np) then
       lambda = lambda + 1
       go to 40
    endif
    call sDCM%rotor (pin, v, pakv, pakst)
    temp = two*pakvm*results%pmemo(4,m1temp)
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "187"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    ctm1 = (results%pmemo(4,m2temp)**2 - results%pmemo(4,m1temp)**2 - pakvm**2)/temp
    temp = two*pakvm*results%pmemo(4,m2temp)
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "193"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    ctm2 = (results%pmemo(4,m1temp)**2 - results%pmemo(4,m2temp)**2 - pakvm**2)/temp
    fm1 = twpi*sDCM%rang()
    fm2 = pi + fm1
    temp = one - ctm1**2
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "201"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    stm1 = sqrt(temp)
    temp = one - ctm2**2
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "207"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    stm2 = sqrt(temp)
    cfm1 = cos(fm1)
    sfm1 = sin(fm1)
    cfm2 = cos(fm2)
    sfm2 = sin(fm2)
    pl(1) = results%pmemo(4,m1temp)*stm1*cfm1
    pl(2) = results%pmemo(4,m1temp)*stm1*sfm1
    pl(3) = results%pmemo(4,m1temp)*ctm1
    call sDCM%rotor (pakst, v, pl, plst)
    results%pmemo(1,m1temp) = plst(1)
    results%pmemo(2,m1temp) = plst(2)
    results%pmemo(3,m1temp) = plst(3)
    pl(1) = results%pmemo(4,m2temp)*stm2*cfm2
    pl(2) = results%pmemo(4,m2temp)*stm2*sfm2
    pl(3) = results%pmemo(4,m2temp)*ctm2
    call sDCM%rotor (pakst, v, pl, plst)
    results%pmemo(1,m2temp) = plst(1)
    results%pmemo(2,m2temp) = plst(2)
    results%pmemo(3,m2temp) = plst(3)
    return


50  nd = nd + 1
    if (nd < 100) then
       go to 20
    else
       kp = 2
    endif

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'direct.f90' line(s) ", A)
1100 format("Square root error prevented in 'direct.f90' line(s) ", A)
! ======================================================================
  end subroutine direct8
