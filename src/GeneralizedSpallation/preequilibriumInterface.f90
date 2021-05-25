
  subroutine preequilibriumInterface (gsmObj, gsmRxn, gsmNucleus, &
       & exciton)

! ======================================================================
!
!    Interface to create a preequilibrium class and then simulate it.
!
!    Called by: simulateDecay
!
!    Calls: Preequilibrium, simulate, ststcs, restor1
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Edited by A. J. Sierk,  LANL  T-2  November, 1997.
!    Modified by AJS, February-March, 1999.
!    Modified by SGM to include reduced masses, 1999
!    Modified by SGM at 06/16/2000
!    Last modification by SGM of 04/27/01, to get the module for CEM2k
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.  
!    Modified by K. K. Gudima, December, 2004.  
!    Edited by AJS, January, 2005.
!    Edited by SGM 07/09/06 to account for the KKG 06/23/06 changes 
!      to use Fermi break-up model in Preco and Evap when A < 13.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Modified by LMK, 2014-2015, for expanison of preeq to heavy ions
! 
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64, int64
    use gsm_params, only: zro, thousand
    use preequilibriumClass, only : &
         & Preeq => Preequilibrium, &
         & preequilibriumResults, newPreequilibriumResults, &
         & preequilibriumFragment, residualNucleus, preeqExcitonData, &
         ! Flags for progeny origin:
         & preequilibriumProgeny, fermiBreakUpProgenyType, &
         ! Flags for state of the residual nucleus:
         & residualState, compoundState, stableState, fermiBreakupState
    use gsm_derived_types, only: nucleus

    
    ! Excited Residual Nucleus (ERN) Properties:
    implicit none
    class(GSM),         intent(inout) :: gsmObj
    class(GSMReaction), intent(inout) :: gsmRxn
    type(nucleus),      intent(inout) :: gsmNucleus
    type(excitonData),  intent(in   ) :: exciton

    integer(int32) :: kfr, particleID, prodMech

    ! For the simulation:
    type(residualNucleus) :: excitedResidual
    type(preeqExcitonData) :: preeqExc
    type(preequilibriumFragment), dimension( nint(gsmNucleus%ap) ) :: progenyBnk
    type(preequilibriumResults) :: preeqRes

    type(GSMResults), pointer:: results => NULL()

! ======================================================================

    ! Setup residual data type based on the passed in residual
    excitedResidual%numBaryons  = gsmNucleus%ap
    excitedResidual%numProtons  = gsmNucleus%zp
    excitedResidual%kinEnergy   = gsmNucleus%up
    ! (Linear Momentum)
    excitedResidual%linearMomX = gsmNucleus%pnx
    excitedResidual%linearMomY = gsmNucleus%pny
    excitedResidual%linearMomZ = gsmNucleus%pnz
    ! (Angular Momentum)
    excitedResidual%angularMom(1) = gsmNucleus%elx
    excitedResidual%angularMom(2) = gsmNucleus%ely
    excitedResidual%angularMom(3) = gsmNucleus%elz
    excitedResidual%angMomFlag    = gsmNucleus%ln
    ! (Fission barrier height)
    excitedResidual%fissBarr = gsmNucleus%bf0

    ! Exciton data:
    preeqExc%numTotal    = exciton%numTotal
    preeqExc%numNucleons = dble(exciton%numProtons + exciton%numNeutrons)
    preeqExc%numProtons  = dble(exciton%numProtons)
    preeqExc%numHoles    = dble(exciton%numHoles)


    ! Construct results object:
    preeqRes = newPreequilibriumResults( progenyBnk )


    ! Simulate preequilibrium physics
    call gsmRxn%models%preeq%simulate(excitedResidual, preeqExc, preeqRes)


    ! Check for errors during the simulation:
    if ( preeqRes%simState > 0 ) then
       if ( preeqRes%simState == 1 ) then
          ! Preeq. process didn't finish (to few available progeny elements in the bank)
          write(gsmObj%io%message, 3000)
          call gsmObj%io%print(3, 3, gsmObj%io%message)
       else if ( preeqRes%simState == 10 ) then
          ! Preeq. object wasn't constructed
          write(gsmObj%io%message, 3600)
          call gsmObj%io%print(1, 1, gsmObj%io%message)
          write(gsmObj%io%message, 3999)
          call gsmObj%io%print(1, 1, gsmObj%io%message)
          stop
       else if ( preeqRes%simState == 11 ) then
          ! Preeq. data object failed to construct
          write(gsmObj%io%message, 3700)
          call gsmObj%io%print(1, 1, gsmObj%io%message)
          write(gsmObj%io%message, 3999)
          call gsmObj%io%print(1, 1, gsmObj%io%message)
          stop
       else if ( preeqRes%simState == 12 ) then
          ! PreeqRes object failed to construct
          write(gsmObj%io%message, 3800)
          call gsmObj%io%print(1, 1, gsmObj%io%message)
          write(gsmObj%io%message, 3999)
          call gsmObj%io%print(1, 1, gsmObj%io%message)
          stop
       else
          ! Unknown error
          write(gsmObj%io%message, 3900) preeqRes%simState
          call gsmObj%io%print(2, 2, gsmObj%io%message)
          write(gsmObj%io%message, 3999)
          call gsmObj%io%print(2, 2, gsmObj%io%message)
       end if
    end if

    ! Interface simulation results now:
    results => gsmRxn%results

    ! Obtain progeny from class
    progenyLoop: do kfr = 1, preeqRes%numProgeny

       ! Ensure bank can hold fragment:
       if ( results%numProgeny > results%maxProgenyM1 ) then
          write(gsmObj%io%message, 2000)
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          write(gsmObj%io%message, 2010) preeqRes%numProgeny - kfr + 1
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          exit progenyLoop
       end if

       ! Obtain interim results:
       ! (particle ID)
       particleID = 0
       if ( progenyBnk(kfr)%numProtons <= 2 .and. progenyBnk(kfr)%numBaryons <= 5 ) then
          if ( progenyBnk(kfr)%numProtons == 0 ) then
             if ( progenyBnk(kfr)%numBaryons == 1 ) then
                particleID = 1   ! Neutron
             end if
          else if ( progenyBnk(kfr)%numProtons == 1 ) then
             if ( progenyBnk(kfr)%numBaryons == 1 ) then
                particleID = 2   ! Proton
             else if ( progenyBnk(kfr)%numBaryons == 2 ) then
                particleID = 3   ! Deuteron
             else if ( progenyBnk(kfr)%numBaryons == 3 ) then
                particleID = 4   ! Triton
             end if
          else if ( progenyBnk(kfr)%numProtons == 2 ) then
             if ( progenyBnk(kfr)%numBaryons == 3 ) then
                particleID = 5   ! Helion (He-3)
             else if ( progenyBnk(kfr)%numBaryons == 4 ) then
                particleID = 6   ! Helium
             end if
          end if
       else
          ! Use expanded formula
          particleID = nint (   thousand * progenyBnk(kfr)%numProtons + &
               & ( progenyBnk(kfr)%numBaryons - progenyBnk(kfr)%numProtons )   )
       end if
       if( particleID == 0 ) then
          particleID = nint (   thousand * progenyBnk(kfr)%numProtons + &
               & ( progenyBnk(kfr)%numBaryons - progenyBnk(kfr)%numProtons )   )
          write(gsmObj%io%message, 2100)
          call gsmObj%io%print(3, 3, gsmObj%io%message)
       end if

       if ( progenyBnk(kfr)%origin == preequilibriumProgeny ) then
          prodMech = 100.0_real64   ! Fragment from preequlibrium emission
       else if ( progenyBnk(kfr)%origin == fermiBreakUpProgenyType ) then
          prodMech = 1500.0_real64   ! Fragment from Fermi Break-up of residual
       else
          ! Error - Invalid origin flag
          write(gsmObj%io%message, 1100) progenyBnk(kfr)%origin
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          prodMech = 100.0_real64
       end if

       ! Bank particle:
       results%numProgeny = results%numProgeny + 1
       results%progenyBnk(results%numProgeny)%numBaryons = progenyBnk(kfr)%numBaryons
       results%progenyBnk(results%numProgeny)%numProtons = progenyBnk(kfr)%numProtons
       results%progenyBnk(results%numProgeny)%kinEnergy  = progenyBnk(kfr)%kinEnergy / thousand
       results%progenyBnk(results%numProgeny)%restMass   = progenyBnk(kfr)%restMass / thousand
       results%progenyBnk(results%numProgeny)%phi        = progenyBnk(kfr)%phi
       results%progenyBnk(results%numProgeny)%theta      = progenyBnk(kfr)%theta
       results%progenyBnk(results%numProgeny)%sinTheta   = sin( progenyBnk(kfr)%theta )
       results%progenyBnk(results%numProgeny)%cosTheta   = cos( progenyBnk(kfr)%theta )
       results%progenyBnk(results%numProgeny)%typeID     = particleID
       results%progenyBnk(results%numProgeny)%prodMech   = prodMech
    end do progenyLoop

    if(results%tallySim) then
       ! Obtain distribution of nuclei prior to equilibration
       results%info%residual(1)%numBaryons = preeqRes%initResidual%numBaryons
       results%info%residual(1)%numProtons = preeqRes%initResidual%numProtons
       results%info%residual(1)%energy     = preeqRes%initResidual%thermEnergy
       results%info%residual(1)%linearMomentum = &
            & sqrt(preeqRes%initResidual%linearMomX**2 + &
            &        preeqRes%initResidual%linearMomY**2 + &
            &        preeqRes%initResidual%linearMomZ**2 ) * thousand
       results%info%residual(1)%angMom = preeqRes%initResidual%angMomFlag
    end if

    ! Obtain statistics for stable nuclei
    if ( preeqRes%residual%state == residualState ) then
       ! Preeq. process did not finish due to not enough progeny storage available
       !    (error addressed above)
    else if ( preeqRes%residual%state == compoundState ) then
       ! Preeq. process finished and a compound nucleus exists (undergo evaporation)
    else if ( preeqRes%residual%state == stableState ) then
       ! Accumulate statistics for stable nuclei before decay (no progeny created)
       if (results%tallySim) then
          if (preeqRes%numProgeny == 0) then
             call gsmObj%ststcs (preeqRes%residual%numBaryons, preeqRes%residual%numProtons, &
                  & preeqRes%residual%thermEnergy, preeqRes%residual%angMomFlag, &
                  & preeqRes%residual%fissBarr, 2)
          end if

          ! Accumulate statistics on stable nucleus and store into secondary arrays
          call gsmObj%ststcs (preeqRes%residual%numBaryons, preeqRes%residual%numProtons, &
               & preeqRes%residual%recEnergy, preeqRes%residual%angMomFlag, &
               & preeqRes%residual%fissBarr, 4)
       end if

       call gsmObj%restor1 (preeqRes%residual%numBaryons, preeqRes%residual%numProtons, &
            & preeqRes%residual%recEnergy, preeqRes%residual%restMass, &
            & preeqRes%residual%normSpeed, results)

    else if ( preeqRes%residual%state == fermiBreakupState ) then
       ! Nucleus underwent Fermi Break-up; set flags accordingly
       results%info%fusion = zro
       results%info%wf = zro
       ! Increment uses of Fermi Break-up
       gsmRxn%outData%ifermi = gsmRxn%outData%ifermi + 1

    else
       ! Unknown state flag - write error and assume compoundState
       write(gsmObj%io%message, 3950)
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       write(gsmObj%io%message, 3955)
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       preeqRes%residual%state = compoundState
    end if


    ! Store compound data locally
    gsmNucleus%ap    = preeqRes%residual%numBaryons
    gsmNucleus%zp    = preeqRes%residual%numProtons
    ! (Energies)
    gsmNucleus%up    = preeqRes%residual%kinEnergy
    gsmNucleus%ue    = preeqRes%residual%thermEnergy
    gsmNucleus%trec  = preeqRes%residual%recEnergy
    gsmNucleus%rstMass = preeqRes%residual%restMass
    ! (linear momentum)
    gsmNucleus%pnx   = preeqRes%residual%linearMomX
    gsmNucleus%pny   = preeqRes%residual%linearMomY
    gsmNucleus%pnz   = preeqRes%residual%linearMomZ
    ! (angular momentum)
    gsmNucleus%elx   = preeqRes%residual%angularMom(1)
    gsmNucleus%ely   = preeqRes%residual%angularMom(2)
    gsmNucleus%elz   = preeqRes%residual%angularMom(3)
    gsmNucleus%ln    = preeqRes%residual%angMomFlag
    ! (fission barrier height)
    gsmNucleus%bf0   = preeqRes%residual%fissBarr
    ! (state of nucleus)
    gsmNucleus%state = preeqRes%residual%state

    return
! ======================================================================
1100 format("Invalid fragment origin flag (", i3, ") was detected during ", &
          & "preequilibrium emission.")
2000 format("The GSM progeny bank was filled during the preequilibrium simulation.")
2010 format("   Unable to bank the remaining ", i3, " particles.")
2100 format("Unexpected particle found. Approximation using ZZZNNN identifier.")
3000 format("Not enough memory was given for progeny storage to allow ", &
          & "the preequilibrium simulation to complete.")
3600 format("The preequilibrium object wasn't constructed.")
3700 format("The preequilibrium data object wasn't constructed.")
3800 format("The preequilibrium results object wasn't constructed.")
3900 format("An unknown error (", i3, ") was flagged during the ", &
          & "preequilibrium simulation.")
3950 format("An unknown state (", i3, ") was flagged for the residual nucleus.")
3955 format("   Assuming residual is a compound nucleus...")
3999 format("   Preequilibrium simulation failed. Stopping...")
! ======================================================================
  end subroutine preequilibriumInterface
