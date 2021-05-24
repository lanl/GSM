
  subroutine selecte (evapObj, j, compound, beta, pekin, ekin, calcVars)

! ======================================================================
!
!   Finds the kinetic energy of the outgoing evaporated particle.
!
!   variable   IN/OUT
!     j    :     IN        index of evaporated particle
!     a    :     IN        residual mass before emission            
!     z    :     IN        charge number of res nuclei before emission 
!     u    :     IN        excitation energy of nuclei before emission
!   pekin  :     IN        random number multiplied by decay width r and
!                          excited-state enhancement factor rr.
!    ekin  :     OUT       randomly selected kinetic energy distributed
!                          according to the emission width.
!
!    Called by: STDCAY
!
!    Calls: PE
!
!    NOTE:  This routine is REALLY inefficient.  Should be replaced
!          by a routine with a decent rejection comparison function.
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, two
    use evaporationFissionData, only: ifa, ifz, ux0

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    integer(int32),     intent(in   ) :: j
    type(evapCompound), intent(in   ) :: compound
    real(real64),       intent(in   ) :: beta
    real(real64),       intent(inout) :: pekin
    real(real64),       intent(  out) :: ekin
    type(evaporationCalculation), intent(inout) :: calcVars

    integer(int32) :: i, iaa, ii, imax, izz, nn
    real(real64)   :: aa, alpha, ax, bet, dppt, e0, ex, pp, ppm, ppo, &
         & ppt, ran, sx, tau, temp, ux, uxx, vj, x, xv

! ======================================================================

    aa = compound%numBaryons - dble(ifa(j))
    izz = nint(compound%numProtons) - ifz(j)
    iaa = nint(aa)
    nn = iaa - izz
    if (iaa == 1) then
       ekin = compound%kinEnergy - calcVars%q(j)
       return
    endif

    bet = beta
    if (j.ne.1) bet = -calcVars%v(j)
    ux = ux0(iaa)
    ex = ux + calcVars%delta(j)
    ax = evapObj%data%ax0(izz,nn)
    tau = evapObj%data%taux0(izz,nn)
    temp = tau
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,2000) '1381'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    tau = one/temp
    sx = evapObj%data%sx0(izz,nn)
    e0 = ex - tau*(log(tau) - 0.25d0*log(ax) - 1.25d0*log(ux) + sx)
    ppt = zro
    imax = 1000
    vj = calcVars%v(j)
    if (j > 6 .and. evapObj%options%inverseParameter == -one) vj = vj*0.1d0
    uxx = compound%kinEnergy - calcVars%q(j) - calcVars%v(j)
    do ii = 1,2
       do i = 0,imax
          temp = dble(imax)
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(evapObj%io%message,2000) '1396'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          x = vj + dble(i)*uxx/temp
          pp = evapObj%pe (compound%kinEnergy, calcVars%q(j), calcVars%delta(j), bet, tau, &
               & e0, ex, calcVars%smalla(j), x, calcVars)
          if (j > 6 .and. x < (calcVars%v(j) + one) .and. &
               & evapObj%options%inverseParameter == -one) then
             temp = dble(izz)
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(evapObj%io%message,2000) '1404'
                call evapObj%io%print(4, 3, evapObj%io%message)
             end if
             alpha = 0.869d0 + 9.91d0/temp
             xv = x - calcVars%v(j)
             temp = xv
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(evapObj%io%message,2000) '1411'
                call evapObj%io%print(4, 3, evapObj%io%message)
             end if
             pp = pp/temp
             temp = one + calcVars%v(j)
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(evapObj%io%message,2000) '1417'
                call evapObj%io%print(4, 3, evapObj%io%message)
             end if
             pp = pp/(temp)*exp(alpha*(xv - one))*x
          endif
          
          if (i == 0) then
             ppo = pp
          else
             temp = two*dble(imax)
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(evapObj%io%message,2000) '1428'
                call evapObj%io%print(4, 3, evapObj%io%message)
             end if
             dppt = calcvars%gj(j)*(ppo + pp)*uxx/(temp)
             ppt = ppt + dppt
             ppo = pp
             ppm = ppt
             if (ppm > pekin) go to 30
          endif
       end do
       temp = calcVars%r(j)*calcVars%rr(j)
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,2000) '1440'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       ran = pekin/(temp)
       if (ii == 1) then
          pekin = ppm*ran
       else
          temp = compound%numBaryons
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(evapObj%io%message,2000) '1449'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          write(evapObj%io%message, 1000) 
          call evapObj%io%print(1, 2, evapObj%io%message)
          write(evapObj%io%message, 1100) j, compound%numBaryons, compound%numProtons, &
               & compound%kinEnergy/temp, ran
          call evapObj%io%print(1, 2, evapObj%io%message)
          write(evapObj%io%message, 1200) pekin, ppm, calcVars%rr(j), calcVars%r(j)
          call evapObj%io%print(1, 2, evapObj%io%message)

          ! Assume kinetic energy is 0
          ekin = 0
          write(evapObj%io%message, 1300) ekin
          call evapObj%io%print(1, 2, evapObj%io%message)
          return
       endif
    end do
30  ekin = x
    return

! ======================================================================
1000 format ("Unable to determine kinetic energy of fragment during ", &
          & "evaporation.")
1100 format ('   jemiss = ',i3,', A = ',f4.0,', Z =',f4.0, &
          & ', Ex/A =  ',f8.6,', random number = ',f8.6)
1200 format ('   pekin, ppm, rr, r =', 4f10.6)
1300 format("   Assuming kinetic energy is ", f8.3, ".")
2000 format("Divide by zero error prevented in 'selecte.f90', line ", A)
! ======================================================================
  end subroutine selecte
