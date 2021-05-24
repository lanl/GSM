
  subroutine evaporationInterface (gsmObj, a, z, ue, trec, pnx, pny, pnz, &
       & ln, bf0, fitaf, fitaf1, gsmRxn)

! ======================================================================
!
!    This routine handles the statistical decay of the compound nucleus
!    by calling the GEMDEC routine to use the GEM2 decay model of
!    Furihata.  It was extracted from the old PRECOF routine from
!    the cem2k+GEM2 code.
!
!    CALLED BY: PRECOF
!    CALLS: GEMDEC MOLNIX
!
!    Written by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by K. K. Gudima, November, 2004.
!    Edited by AJS, January, 2005.
!    Modified by REP, 11/17/05
!    Updated for MARS by K.K. Gudima, March 2007
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!    Edited by CMJ, XCP-3, 2016-2017 (LAQGSM expansion)
! 
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64, error_unit
    use gsm_params, only: zro, one, two, thr, &
         & four, fiv, six, thousand, pi, twpi, emnucb, emnuct
    use evaporationClass,     only: &
         & evaporationResults, newEvaporationResults, &
         & evaporationFragment, &
         & fissionFragment,     &
         & evapCompound,        &
         & evaporationFlag, fermiBreakUpFlag, fissionFlag

    implicit none
    class(GSM),         intent(inout) :: gsmObj
    real(real64),       intent(inout) :: a
    real(real64),       intent(inout) :: z
    real(real64),       intent(inout) :: ue
    real(real64),       intent(inout) :: trec
    real(real64),       intent(in   ) :: pnx
    real(real64),       intent(in   ) :: pny
    real(real64),       intent(in   ) :: pnz
    integer(int32),     intent(in   ) :: ln
    real(real64),       intent(in   ) :: bf0
    real(real64),       intent(in   ) :: fitaf
    real(real64),       intent(in   ) :: fitaf1
    class(GSMReaction), intent(inout) :: gsmRxn

    integer(int32) :: in, inj, iz, izj, nfr, k, kf, prodMech, particleID
    real(real64)   :: cosPhi, cosTheta, e, enfr, eresa, phi, pmod, &
         & remn, sinPhi, sinTheta, theta, temp
    real(real64), dimension(3) :: beta = zro

    ! Used to start the simulation...
    type(evapCompound) :: compound

    ! (for Evaporation results)
    type(evaporationResults)                        :: evapRes
    type(evaporationFragment), dimension( nint(a) ) :: evapBnk
    type(fissionFragment),     dimension(       2 ) :: fissionBnk

    type(GSMResults), pointer:: results => NULL()

! ======================================================================

    iz = max(1, nint(z))
    in = max(1, nint(a - z))
    if (iz > 7 .and. in > 7) then
       remn = a*emnucb + gsmObj%genData%molnixEnergies%defineEnergy(iz, in, 2)/thousand
    else
       remn = a*emnuct + gsmObj%genData%molnixEnergies%defineEnergy(iz, in, 2)/thousand
    endif
!  Total recoil energy and kinetic energy of entire nucleus:
    pmod = sqrt(pnx**2 + pny**2 + pnz**2)
    e = sqrt(pmod**2 + remn**2)
    eresa = (e - remn)*thousand
    temp = a
    if ( abs(temp) < div0Lim ) then
       temp = div0Lim
       write(gsmObj%io%message,2000) "287"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    temp = pmod
    if ( abs(temp) < div0Lim ) then
       temp = div0Lim
       write(gsmObj%io%message,2000) "293-295"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    ! Setup normalized momemtum
    beta(1) = pnx/temp
    beta(2) = pny/temp
    beta(3) = pnz/temp


    ! Construct results data type:
    evapRes = newEvaporationResults( evapBnk, fissionBnk )


    ! Setup compound object
    compound%numBaryons = a
    compound%numProtons = z
    compound%kinEnergy  = ue
    compound%recEnergy  = eresa
    compound%linearMomFrac(:) = beta(:)
    compound%afMultiplier = fitaf
    compound%czMultiplier = fitaf1


    ! Undergo evaporation/fission
    call gsmObj%genModels%evap%simulate ( compound, evapRes )


    ! Store evaporation/fission fragments into spt/parz arrays
    results => gsmRxn%results
    progenyLoop: do nfr = 1, evapRes%numProgeny

       ! Ensure bank isn't full:
       if ( results%numProgeny > results%maxProgenyM1 ) then
          write(gsmObj%io%message, 1000)
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          write(gsmObj%io%message, 1010) evapRes%numProgeny-nfr+1
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          exit progenyLoop
       end if

       ! Obtain interim results:
       izj = nint(evapBnk(nfr)%numProtons)
       inj = nint(evapBnk(nfr)%numBaryons - evapBnk(nfr)%numProtons)
       enfr = dble(inj)
       cosTheta = evapBnk(nfr)%linearMomZ
       if (abs(cosTheta) >= one) then
          cosTheta = sign(one,cosTheta)
          sinTheta = zro
          theta = zro
          if (cosTheta < zro) theta = pi 
       else 
          sinTheta = sqrt(abs(one - cosTheta**2))
          theta = atan2(sinTheta,cosTheta)
       endif
       if (abs(sinTheta) > 1.d-10) then
          temp = sinTheta
          if ( abs(temp) < div0Lim ) then
             temp = div0Lim
             write(gsmObj%io%message,2000) "359, 360"
             call gsmObj%io%print(4, 3, gsmObj%io%message)
          end if
          cosPhi = evapBnk(nfr)%linearMomX/temp
          sinPhi = evapBnk(nfr)%linearMomY/temp
          if (abs(cosPhi) >= one) then
             cosPhi = sign(one,cosPhi)
             if (cosPhi <= zro) phi = pi
             if (cosPhi > zro) phi = zro
             sinPhi = zro
          else
             phi = atan2(sinPhi,cosPhi)
             if (phi < zro) phi = twpi + phi
          endif
       else
          phi = zro
          cosPhi = one
          sinPhi = zro
       endif

       ! Determine the particle 'type'
       particleID = zro
       if (izj == 0 .and. inj == 1) particleID = one
       if (izj == 1 .and. inj == 0) particleID = two
       if (izj == 1 .and. inj == 1) particleID = thr
       if (izj == 1 .and. inj == 2) particleID = four
       if (izj == 2 .and. inj == 1) particleID = fiv
       if (izj == 2 .and. inj == 2) particleID = six
       if (izj > 2 .or.  inj > 2) particleID = thousand*evapBnk(nfr)%numProtons + enfr
       if ( particleID == 0 ) particleID = thousand*evapBnk(nfr)%numProtons + enfr
       ! Production mechanism:
       if ( evapBnk(nfr)%physicsFlag == evaporationFlag ) then
          prodMech = 1000_int32   ! Progeny from compound nucleus
       else if ( evapBnk(nfr)%physicsFlag == fermiBreakUpFlag ) then
          prodMech = 1500_int32   ! Progeny from Fermi Breakup of compound nucleus
       else if ( evapBnk(nfr)%physicsFlag == fissionFlag ) then
          prodMech = 2000_int32   ! Progeny from evaporation of fission fragment
       else
          ! Error; progeny origin wasn't flagged by evaporation class
          write(gsmObj%io%message, 3000)
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          prodMech = 1000_int32
       end if


       ! Tally progeny:
       results%numProgeny = results%numProgeny + 1
       results%progenyBnk(results%numProgeny)%numBaryons = evapBnk(nfr)%numBaryons
       results%progenyBnk(results%numProgeny)%numProtons = evapBnk(nfr)%numProtons
       results%progenyBnk(results%numProgeny)%kinEnergy  = evapBnk(nfr)%kinEnergy / thousand
       results%progenyBnk(results%numProgeny)%restMass   = evapBnk(nfr)%restMass
       results%progenyBnk(results%numProgeny)%phi        = phi
       results%progenyBnk(results%numProgeny)%theta      = theta
       results%progenyBnk(results%numProgeny)%sinTheta   = sinTheta
       results%progenyBnk(results%numProgeny)%cosTheta   = cosTheta
       results%progenyBnk(results%numProgeny)%typeID     = particleID
       results%progenyBnk(results%numProgeny)%prodMech   = prodMech
    end do progenyLoop


    ! Set fission probability
    results%info%wf = evapRes%fissionProbability
    
    ! Progeny stored; now obtain residual information and/or fission information 
    if ( evapRes%compound%physicsFlag == evaporationFlag ) then
       ! Compound nucleus was only evaporated:
       results%info%fusion = zro
       if(results%tallySim) then
          ! Get statistics on evaporated residual
          trec = evapRes%compound%recEnergy
          call gsmObj%ststcs (evapRes%compound%numBaryons, &
               & evapRes%compound%numProtons, &
               & trec, ln, bf0, 4)
       end if

    else if ( evapRes%compound%physicsFlag == fermiBreakUpFlag ) then
       ! Compund underwent Fermi Breakup:
       gsmRxn%outData%ifermi = gsmRxn%outData%ifermi + 1

    else if ( evapRes%compound%physicsFlag == fissionFlag ) then
       ! Compound underwent Fission:
       if(results%tallySim) then
          ! Store pre-fission compound's characteristics in the /residf/ block
          results%info%compound = evapRes%compound

          ! Obtain statistics on pre-fission residual
          results%info%fusion = one
          call gsmObj%ststcs (evapRes%compound%numBaryons, &
               & evapRes%compound%numProtons, &
               & evapRes%compound%kinEnergy, &
               & ln, bf0, 5)
       
          ! Store fission fragment data
          results%info%ifiss = 1_int32
          do k = 1, min(2, evapRes%numFissionFragments)
             if ( evapRes%fissionBnk(k)%physicsFlag == evaporationFlag ) then
                ! Fission fragment underwent evaporation
                results%info%fissFrag(k) = evapRes%fissionBnk(k)
             else if ( evapRes%fissionBnk(k)%physicsFlag == fermiBreakUpFlag ) then
                ! Fermi Break-up of fission fragment occurred
                gsmRxn%outData%ifermi = gsmRxn%outData%ifermi + 1
             else if ( evapRes%fissionBnk(k)%physicsFlag == fissionFlag ) then
                ! Fission fragment was evaporated and flagged to undergo fission again
                write(gsmObj%io%message, *) "A fission fragment underwent fission (not handled)!"
                call gsmObj%io%print(1, 2, gsmObj%io%message)
             else
                ! No fission was undergone
                write(gsmObj%io%message, *) "Unexpected logic found in GEM2 evaporation interface."
                call gsmObj%io%print(1, 2, gsmObj%io%message)
             end if
          end do
       end if

    else
       ! Error: Compound physics processes were not accurately tracked
       write(gsmObj%io%message,*) "The compound was flagged to undergo an unkown process in the evaporation interface."
       call gsmObj%io%print(1, 2, gsmObj%io%message)
       return
    end if

    return

! ======================================================================
1000 format("The GSM progeny bank was filled during the evaporation simulation.")
1010 format("   Unable to bank the remaining ", i3, " progeny.")
2000 format("Divide by zero error prevented in ", &
          & "'evapInterface.f90', line ", A)
3000 format("Unable to determine progeny origin from ", &
          & "the evaporation and fission process. Progeny will be assumed ", &
          & "to originate from evaporation of the compound nucleus.")
! ======================================================================
  end subroutine evaporationInterface
