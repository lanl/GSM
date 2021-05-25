
  subroutine constructGeneralData( gsmObj )

! ==============================================================================
!
! This function constructs the general physics data used by GSM in the generalPhysicsData
! data type
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use molnixClass, only: newMolnix
    use fissionBarrierClass, only: newFissionBarrier
    use evaporationDataClass, only: newEvaporationData
 
    implicit none
    class(GSM), intent(inout) :: gsmObj

! ==============================================================================

    gsmObj%genData%numConstructionErrors = 0


    ! Molnix object:
    gsmObj%genData%molnixEnergies = newMolnix( &
         & gsmObj%options%molnixOpts, &
         & gsmObj%io%print )
    if ( .not.gsmObj%genData%molnixEnergies%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) "The Molnix object failed to properly construct."
       call gsmObj%io%print(1, 1, gsmObj%io%message)
       gsmObj%genData%numConstructionErrors = gsmObj%genData%numConstructionErrors + 1
    end if
    gsmObj%options%molnixOpts = gsmObj%genData%molnixEnergies%queryOptions()


    ! Fission Barrier object:
    gsmObj%genData%fissBarrier = newFissionBarrier( &
         & gsmObj%genData%molnixEnergies, &
         & gsmObj%options%fissBarOpts, gsmObj%io%print )
    if ( .not.gsmObj%genData%fissBarrier%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) "The Fission Barrier object failed to properly construct."
       call gsmObj%io%print(1, 1, gsmObj%io%message)
       gsmObj%genData%numConstructionErrors = gsmObj%genData%numConstructionErrors + 1
    end if
    gsmObj%options%fissBarOpts = gsmObj%genData%fissBarrier%queryOptions()


    ! Evaporation Data object:
    gsmObj%genData%evap = newEvaporationData( &
         & gsmObj%options%evapDataOpts, &
         & gsmObj%io%print )
    if ( .not.gsmObj%genData%evap%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) "The EvaporationData object failed to properly construct."
       call gsmObj%io%print(1, 1, gsmObj%io%message)
       gsmObj%genData%numConstructionErrors = gsmObj%genData%numConstructionErrors + 1
    end if
    gsmObj%options%evapDataOpts = gsmObj%genData%evap%queryOptions()


    ! Check how many errors occurred:
    if ( gsmObj%genData%numConstructionErrors /= 0 ) then
       write(gsmObj%io%message, 1100) gsmObj%genData%numConstructionErrors
       call gsmObj%io%print(1, 1, gsmObj%io%message)
    end if

    return
! ==============================================================================
1000 format(A)
1100 format(i2, "  data objects failed to construct within GSM.")
! ==============================================================================
  end subroutine constructGeneralData
