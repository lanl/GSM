
  subroutine stdcay (evapObj, compound, iaf, izf, iflag, fisinh, &
       & ekin, bet1, results, calcVars)

! ======================================================================
!
! <STDCAY>
!  Calculate evaporation process and fission process
!    - determine if emission occurs
!    - determine if fission occurs
!        -if fission , call FISGEM
!        -if no fission, determine kinetic energy, etc. for emittor
!                        and caluculate recoil energy, etc.
!
!=====================================================================
!   Called by:  GEMDEC
!
!   Calls:
! <Subroutine>
!     GAMMA  :  Decay width calculation
!    FISGEM  :  Fission calculation
!   SELECTE  :  Select kinetic energy in the CM system
!--------------------------------------------------------------------
! <Function>
!  fprob     : Calculate fission probabirity
!=====================================================================
! <Variables>   (Units of MeV for energy)
!     a    : mass of parent nucleus -> residual nucleus   (IN and OUT)
!     z    : charge of parent nucleus    -> residual      (IN and OUT)
!     u    : excited energy of parent nucleus-> residual  (IN and OUT)
!   iflag  : 1 = no more emission                         (OUT)
!          : 2 = fission
!          : 3 = emission occur
!  fisinh  : true=fission occur, false=no fission         (OUT)
!   erec   : recoil energy in the lab system              (IN and OUT)
!   ekin   : kinetic energy in the lab system             (OUT)
!   bet0   : unit vector of recoil momentum               (IN and OUT)
!   bet1   : unit vector of momentum of emitted particle  (IN and OUT)
!/////////////////////////////////////////////////////////////////////
!
!   "Last" change: 13-AUG-2003 by NVMokhov
!   Modified by A. J. Sierk, LANL T-16, September-October, 2003.
!   Modified by K. K. Gudima, December, 2004.
!   Edited by AJS, January, 2005.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!   Edited by CMJ, XCP-3, July 2018 (Evap. class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, two, amu, twpi
    use evaporationFissionData, only: ifa, ifz

    implicit none
    class(Evaporation),  intent(inout) :: evapObj
    type(evapCompound),  intent(inout) :: compound
    integer(int32),      intent(  out) :: iaf
    integer(int32),      intent(  out) :: izf
    integer(int32),      intent(  out) :: iflag
    real(real64),        intent(inout) :: ekin
    real(real64),        intent(inout) :: bet1(3)
    logical,             intent(inout) :: fisinh
    type(evaporationResults), intent(inout) :: results
    type(evaporationCalculation), intent(inout) :: calcVars
    
    integer(int32) :: ia, iz, j, jemiss, k
    real(real64)   :: am, amf, amr, b1, beta, e1cm, e2cm, gam, gpo, p1x, &
         & p1y, p1z, p2x, p2y, p2z, pcm, pcmx, pcmy, pcmz, pe2, pekin, &
         & ph, pin, pr2, pran, pres, pv, redm, sum, temp, th, tr, &
         & uran, vres, vtr

! ======================================================================

!  Initialization
    iaf = 0
    izf = 0
    ia = nint(compound%numBaryons)
    iz = nint(compound%numProtons)

!   Evaporation calculation starts:
10  continue 
    uran = evapObj%rang()

!  Calculate decay width:
    call evapObj%gamma (compound, uran, beta, calcVars)
    if (calcVars%sigma <= zro) then
!  No more emission, stop evaporation simulation:
       iflag = one
       return
    endif

!    Select ejectile
    sum = zro
    do j = 1, evapObj%options%numEvapType
       k = j
       sum = sum + calcVars%r(j)/calcVars%sigma
       if (sum >= uran) go to 20 ! Found ejectile, evaporate it
    end do
    go to 10 ! Could not find random particle to evaporate, try again

! Start of evaporation process
20  continue
    jemiss = k

!   Fission calculation :
!   fission only occurs once (fission never occurs to post-fission 
!   nucleus):
    if (jemiss == 1 .and. compound%numProtons > 65.d0 .and. .not.fisinh) then
       results%fissionProbability = evapObj%fprob(compound, beta, calcVars)
       b1 = evapObj%rang()

       ! Fission if probability large enough
       if (results%fissionProbability > b1) then

          ! State fission occurred
          iflag = 2

          ! Check if fission already occurred (if so, only flag; if not, undergo fission)
          if ( results%compound%physicsFlag /= fissionFlag ) then

             ! Simulate fission:
             fisinh = .true.
             call evapObj%fisgem (compound, fisinh, results)

             ! If fission was successful, store data:
             if ( fisinh ) then


                ! Store pre-fission compound nucleus information, set physics flag
                results%compound%numBaryons   = compound%numBaryons
                results%compound%numProtons   = compound%numProtons
                results%compound%kinEnergy    = compound%kinEnergy
                results%compound%recEnergy    = compound%recEnergy
                !  Recoil momentum (nonrel):
                am = compound%numBaryons*amu + evapObj%energy (iz, ia)
                results%compound%linearMomTot = sqrt( abs(compound%recEnergy**2 + two*am*compound%recEnergy) )
                ! Flag that compound is a pre-fission compound
                results%compound%physicsFlag  = fissionFlag
             else
                ! Try to fission or evaporate (i.e. obtain new fission probability)
                go to 20
             end if
          end if

          return
       endif
    endif

!  Set A & Z of the ejectile (evaporation):
    iaf = ifa(jemiss)
    izf = ifz(jemiss)
!  Emission occured:
    iflag = 3                             

!  Select kinetic energy in the CM system:
    pran = evapObj%rang()
    pekin = pran*calcVars%r(jemiss)*calcVars%rr(jemiss)
    call evapObj%selecte (jemiss, compound, beta, pekin, ekin, calcVars)
    if (ekin < zro) then
       write (evapObj%io%message, 1100)
       call evapObj%io%print(1, 1, evapObj%io%message)
       ekin = 0
    endif

!   Velocity of Parent nucleus in the LAB system: vres*bet0
    am = compound%numBaryons*amu + evapObj%energy (iz, ia)
    pres = sqrt(abs(compound%recEnergy**2 + two*am*compound%recEnergy))
    temp = compound%recEnergy + am
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
        write(evapObj%io%message,1000) "160"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    vres = pres/(temp)

!   Momentum in the CM system : pcmx, pcmy, pcmz
    amf = dble(iaf)*amu + evapObj%energy (izf, iaf) - dble(izf)*0.511004d0
    amr = dble(ia - iaf)*amu + evapObj%energy (iz-izf, ia-iaf)
    if (am < div0Lim .and. am > -div0Lim) then
       am = div0Lim
        write(evapObj%io%message,1000) "169"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    redm = amf*amr/am
!    redm = amu*emured(ia, iaf)
    pcm = sqrt(abs(ekin**2 + two*redm*ekin))
    th  = acos(two*evapObj%rang() - one)
    ph  = twpi*evapObj%rang()
    pcmx = pcm*sin(th)*cos(ph)
    pcmy = pcm*sin(th)*sin(ph)
    pcmz = pcm*cos(th)

!   Residual nucleus, A, Z, and excitation energy:
    compound%numBaryons = compound%numBaryons - dble(iaf)
    compound%numProtons = compound%numProtons - dble(izf)
    compound%kinEnergy  = compound%kinEnergy  - calcVars%q(jemiss) - ekin
! LMK, 07/2012, Temporary fix for the rare negative-u bug (modified by CMJ)
!    u = max(zro,u)
    if ( compound%kinEnergy < zro ) then
       write(evapObj%io%message,2000)
       call evapObj%io%print( 3, 3, evapObj%io%message)
       write(evapObj%io%message,2010) &
            compound%numBaryons, compound%numProtons, compound%kinEnergy
       call evapObj%io%print( 3, 3, evapObj%io%message)
       compound%kinEnergy = zro
    end if

    gam = sqrt(abs(one - vres**2))
    if (gam < div0Lim .and. gam > -div0Lim) then
       gam = div0Lim
       write(evapObj%io%message,1000) "191"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    gam = one/gam
    pv = pcmx*compound%linearMomFrac(1) + pcmy*compound%linearMomFrac(2) + pcmz*compound%linearMomFrac(3)
    pv = pv*vres
    temp = gam + one
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) "199"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    gpo = gam*pv/(temp)
!   Boost ejectile momentum to Lab system : p1x,p1y,p1z
    e1cm = sqrt(pcm**2 + amf**2)
    tr = gam*(e1cm + gpo)
    vtr = vres*tr
    p1x = pcmx + vtr*compound%linearMomFrac(1)
    p1y = pcmy + vtr*compound%linearMomFrac(2)
    p1z = pcmz + vtr*compound%linearMomFrac(3)
!   Kinetic Energy and velocity of ejectile in Lab system: ekin, bet1
    pe2 = p1x**2 + p1y**2 + p1z**2
    pin = one/sqrt(abs(pe2))
    ekin =  sqrt(abs(amf**2 + pe2)) - amf
    bet1(1) = pin*p1x
    bet1(2) = pin*p1y
    bet1(3) = pin*p1z
!   Boost residual momentum to Lab system: p2x,p2y,p2z
    e2cm = sqrt(pcm**2 + amr**2)
    tr = gam*(e2cm + gpo)
    vtr = vres*tr
    p2x = -pcmx + vtr*compound%linearMomFrac(1)
    p2y = -pcmy + vtr*compound%linearMomFrac(2)
    p2z = -pcmz + vtr*compound%linearMomFrac(3)
!   Kinetic Energy and velocity of Residual in Lab system: erec, bet0
    pr2 = p2x**2 + p2y**2 + p2z**2
    if (pr2 < div0Lim .and. pr2 > -div0Lim) then
       pr2 = div0Lim
       write(evapObj%io%message,1000) "227"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    pin = one/sqrt(abs(pr2))
    compound%recEnergy =  sqrt(abs(amr**2 + pr2)) - amr
    compound%linearMomFrac(1) = pin*p2x
    compound%linearMomFrac(2) = pin*p2y
    compound%linearMomFrac(3) = pin*p2z

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'stdcay.f90', line ", A)
1100 format("Kinetic energy can not be determined in subroutine 'selecte'.")
2000 format("Negative kinetic energy was sampled for an evaporated fragment.")
2010 format("   Fragment has (A = ", f5.2, ", Z = ", f5.2, &
          & ", Ek = ", es15.8, ")")
! ======================================================================
  end subroutine stdcay

