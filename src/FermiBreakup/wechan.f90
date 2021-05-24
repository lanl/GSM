
  function wechan (fbuObj, ja, k, ntv, daz, tq)
 
! ======================================================================
!
!    Called by: RAZVAL
!
!    Calls: WEPT
!
!    Last change: 13-Aug-2003 BY NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Modified 03 May 2006 by R E Prael to make the Coulomb barrier slightly
!      penetrable (purely an artificial procedure that reduces the channel
!      probability by 10.d-06 below the Coulomb energy.) Necessary for the 
!      breakup of Be8 at low excitation.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fermiBreakupParams, only : zro, one, eulers_number
    use fermiBreakUpData, only : ms, gaf, gafSize, vak, min_fermi_AData

    implicit none
    class(FermiBreakup), intent(inout) :: fbuObj
    integer(int32),      intent(in   ) :: ja
    integer(int32),      intent(in   ) :: k
    integer(int32),      intent(in   ) :: ntv(fbuObj%options%recNumNucleons)
    real(real64),        intent(in   ) :: daz
    real(real64),        intent(in   ) :: tq

    real(real64)                       :: wechan

    integer(int32) :: i, ia, ik, in, iz, j, ki, mrs, tempk
    real(real64)   :: akm, akv, bnq, gam, rmq, rpm, spm, temp, teq, tn, &
         & vmk, vtk

! ======================================================================
 
    bnq = daz
    rmq = -0.1d0
    spm = one
    akm = one
    akv = dble(ja)   ! A of the current residual? (use no scaling)
    if (fbuObj%options%akpScalingFlag > zro) akv = fbuObj%options%akpScalingFlag*dble(k)   ! Scaling for the considered fragment?
    akv = akv*vak   ! Could this be the system volume? (decaying or ground state)?
    do i = 1, min(k, fbuObj%options%recNumNucleons)
       ia = ntv(i)/100
       iz = ntv(i) - ia*100
       in = ia - iz + 1
       akm = akm*dble(ia)
       if ( in > min_fermi_AData .or. iz > (min_fermi_AData-1) ) then
          ! Array ms(X,Y) will be exceeded; print and approximate instead
          write(fbuObj%io%message, 2100) in, iz
          call fbuObj%io%print(3, 3, fbuObj%io%message)
          write(fbuObj%io%message, 2110)
          call fbuObj%io%print(3, 3, fbuObj%io%message)
          if ( in > min_fermi_AData     ) in = min_fermi_AData     ! Approximate in neutron number
          if ( iz > (min_fermi_AData-2) ) iz = min_fermi_AData-2   ! Approximate in proton  number
       end if
       spm = spm*dble(ms(in,iz+1))
       bnq = bnq - fbuObj%wept (ia, iz)
    end do
    tn = bnq
!  REP, 03 May, 2006:
    if (tn > zro) then
       temp = dble(ja)
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(fbuObj%io%message,1000) "352"
          call fbuObj%io%print(4, 3, fbuObj%io%message)
       end if
       akm = akm/temp
       temp = akm
       if (temp < 0.0d0) then
          temp = 0.01d0
          write(fbuObj%io%message,1100) "358"
          call fbuObj%io%print(4, 3, fbuObj%io%message)
       end if
       akm = akm*sqrt(temp)*spm
       if (k > 2) then
          ki = k - 1
          rpm = one
          vmk = one
          teq = eulers_number*tn/(1.5d0*dble(ki) - one)
          temp = teq
          if (temp < 0.0d0) then
             temp = 0.01d0
             write(fbuObj%io%message,1100) "369"
             call fbuObj%io%print(4, 3, fbuObj%io%message)
          end if
          vtk = akv*teq*sqrt(temp)
          do i = 1,ki
             mrs = 1
             ik = i + 1
             vmk = vmk*vtk
             do j = ik,k
                if (ntv(i) == ntv(j)) mrs = mrs + 1
             end do
             rpm = rpm*dble(mrs)
          end do
          tempk = k
          if ( tempk > gafSize ) then
             ! Array 'gaf' will be exceeded; print warning and approximate
             write(fbuObj%io%message,2000) tempk
             call fbuObj%io%print(3, 3, fbuObj%io%message)
             tempk = min_fermi_AData
          end if
          gam = gaf(tempk)
          temp = teq*rpm
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(fbuObj%io%message,1000) "385"
             call fbuObj%io%print(4, 3, fbuObj%io%message)
          end if
          rmq = vmk*gam*akm/(temp)
       else
          temp = tn
          if (temp < 0.0d0) then
             temp = 0.01d0
             write(fbuObj%io%message,1100) "392"
             call fbuObj%io%print(4, 3, fbuObj%io%message)
          end if
          rmq = 1.1283792d0*akv*akm*sqrt(temp)
          if (ntv(1) == ntv(2)) rmq = 0.5d0*rmq
       endif
!  REP, 03 May, 2006:
       if (tn < tq) rmq = rmq*1.d-06
    endif
    wechan = rmq

    return
! ======================================================================
1000 format("Divide by zero error prevented in 'wechan.f90', line ", A)
1100 format("Square root error prevented in 'wechan.f90', line ", A)
2000 format("Fermi break-up library 'gaf' data array ", &
          & "exceeded (k=", i3, "). Using closest valid value.")
2100 format("The Fermi Break-Up model 'ms' data array ", &
          & "exceeded (N=", i3, ", Z=", i3, ").")
2110 format("   Using closest valid value.")
! ======================================================================
  end function wechan
