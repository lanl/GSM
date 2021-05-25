
  subroutine standardDCMInterface( gsmObj, rxnData, sDCMProj, &
       & results)

! ==============================================================================
!
! Interface to the Standard DCM Class
!
!
! Written by CMJ, XCP-3 (01/2019)!
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro, hlf, twpi
    use standardDCMClass, only: StandardDCM, &   ! Class type and interface
         & sDCMProgeny, &
         & StandardDCMResults, newStandardDCMResults, &   ! Objects for class (results)
         & sDCMProjectile                                 ! Target and Projectile types

    implicit none
    class(GSM),            intent(inout) :: gsmObj
    type(rxnSpecificData), intent(inout) :: rxnData
    type(sDCMProjectile),  intent(inout) :: sDCMProj
    type(GSMResults),      intent(inout) :: results

    integer(int32) :: i, particleID
    real(real64)   :: phi, theta

    ! Standard DCM results object(s):
    type(sDCMProgeny), dimension( results%maxProgeny ) :: progenyBnk   ! Make memory for progeny bank
    type(StandardDCMResults) :: sDCMResults

! ==============================================================================

    ! NOTE: The Standard DCM Object is constructed elsewhere

    ! Initialize a results object (give access to progeny array)
10  continue
    progenyBnk(:) = sDCMProgeny()   ! Reset progeny array (when re-doing event)
    sDCMResults = newStandardDCMResults(progenyBnk)


    ! Perform simulation
    call gsmObj%genModels%sDCM%interact( sDCMProj, rxnData%sDCMTargetData, &
         & sDCMResults )


    ! Verify that simulation didn't fail
    if ( sDCMResults%simState /= 0 ) then
       write( gsmObj%io%message, 2100) sDCMResults%simState
       call gsmObj%io%print(3, 3, gsmObj%io%message)
       write(gsmObj%io%message, 3000)
       call gsmObj%io%print(3, 3, gsmObj%io%message)

       ! Special handling for various error states:
       if ( sDCMResults%simState == 1 ) then
          ! An unphysical residual nucleus was created. The sDCM returns the
          ! pre-error residual nucleus and progeny particle bank (before including
          ! the error-causing progeny).
          ! No action necessary (sDCM handles error)
          !    NOTE: GSM COULD RESTART EVENT IF DESIRED.
          if ( sDCMResults%residual%numProtons >= sDCMresults%residual%numBaryons ) then
             write(gsmObj%io%message, 2150)
             call gsmObj%io%print(3, 3, gsmObj%io%message)
             go to 10
          end if
       else if ( sDCMResults%simState == 2 ) then
          ! Flags an error obtaining angular distributions during INC
       else if ( sDCMResults%simState >= 5 ) then
          ! Fatal error that ONLY a developer should see
          stop
       end if
    end if


    ! Increment the number of times the sDCM had to restart the simulation
    results%modelUsage%numSDCMRestarts = &
         & results%modelUsage%numSDCMRestarts + &
         & sDCMResults%numRestarts(1) + &
         & sDCMResults%numRestarts(2)


   ! Check for elastic events
    results%numElasticEvents = results%numElasticEvents + sDCMResults%numElastic
    if ( sDCMResults%numElastic > 0 ) return



    ! Store the created progeny into GSM progeny bank
    tallyProgeny: do i = 1, sDCMResults%numProgeny

       ! Check for full particle bank:
       if ( results%numProgeny > results%maxProgenyM1 ) then
          write(gsmObj%io%message, 2200) sDCMResults%numProgeny - i + 1
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          exit tallyProgeny
       end if

       ! Obtain interim results for the progeny:
       theta = atan2( progenyBnk(i)%sinTheta, &
            & progenyBnk(i)%cosTheta )
       if (theta < zro) theta = twpi + theta
       phi = atan2( sDCMResults%progenyBnk(i)%sinPhi, &
            & sDCMResults%progenyBnk(i)%cosPhi )
       if (phi < zro) phi = twpi + phi
       particleID = 0
       if ( sDCMResults%progenyBnk(i)%restMass <= hlf ) then
          if ( sDCMResults%progenyBnk(i)%numProtons == -1 ) then
             particleID = 7   ! Pi-
          else if ( sDCMResults%progenyBnk(i)%numProtons == 0 ) then
             particleID = 8   ! Pi0
          else if ( sDCMResults%progenyBnk(i)%numProtons == 1 ) then
             particleID = 9   ! Pi+
          end if
       else
          if ( sDCMResults%progenyBnk(i)%numProtons == 0 ) then
             particleID = 1   ! Neutron
          else if ( sDCMResults%progenyBnk(i)%numProtons == 1 ) then
             particleID = 2   ! Proton
          end if
       end if
       if ( particleID <= 0 ) then
          ! Unexpected particle produced during the SDCM simulation. 
          write(gsmObj%io%message, 2300)
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          write(gsmObj%io%message, 2310)
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          particleID = 1   ! Assume ID is a neutron (keep actual characteristics)
       end if

       ! Tally the progeny:
       results%numProgeny = results%numProgeny + 1
       results%progenyBnk( results%numProgeny )%numBaryons = progenyBnk(i)%numBaryons
       results%progenyBnk( results%numProgeny )%numProtons = progenyBnk(i)%numProtons
       results%progenyBnk( results%numProgeny )%kinEnergy  = progenyBnk(i)%kinEnergy
       results%progenyBnk( results%numProgeny )%restMass   = progenyBnk(i)%restMass
       results%progenyBnk( results%numProgeny )%phi        = phi
       results%progenyBnk( results%numProgeny )%theta      = theta
       results%progenyBnk( results%numProgeny )%sinTheta   = progenyBnk(i)%sinTheta
       results%progenyBnk( results%numProgeny )%cosTheta   = progenyBnk(i)%cosTheta
       results%progenyBnk( results%numProgeny )%typeID     = particleID
       results%progenyBnk( results%numProgeny )%prodMech   = progenyBnk(i)%index

    end do tallyProgeny

    ! Residual nucleus is as follows:
    results%targRes%numBaryons = sDCMResults%residual%numBaryons
    results%targRes%numProtons = sDCMResults%residual%numProtons
    results%targRes%kinEnergy  = sDCMResults%residual%kinEnergy
    do i = 1, 3
       results%targRes%linearMom(i)  = sDCMResults%residual%linearMom(i)
       results%targRes%angularMom(i) = sDCMResults%residual%angularMom(i)
    end do
    

    ! Exciton data is as follows:
    results%targExc%numTotal    = sDCMResults%excitons%numExcProt + &
         & sDCMResults%excitons%numExcNeut + sDCMResults%excitons%numExcHoles
    results%targExc%numProtons  = sDCMResults%excitons%numExcProt
    results%targExc%numNeutrons = sDCMResults%excitons%numExcNeut
    results%targExc%numHoles    = sDCMResults%excitons%numExcHoles

    return
! ==============================================================================
2100 format("GSM's Standard DCM simulation for an event", &
          & " ended in an error state (=", i1, ").")
2150 format("   The event will be restarted.")
2200 format("The particle bank in GSM has been exceeded by ", i4, " particle(s).")
2300 format("An unexpected progeny was created during the Standard DCM simulation.")
2310 format("   The progeny will be arbitrarily identified as a neutron.")
3000 format("   Be wary of simulation results.")
! ==============================================================================
  end subroutine standardDCMInterface
