
  subroutine constructSpecificData( gsmObj, projObj, targObj, rxnData )

! ==============================================================================
!
! This function constructs the reaction-specific physics data used by
! GSM in the generalPhysicsData data type
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use gsm_params, only: thousand
    use standardDCMDataClass,    only: newStandardDCMData
    use coalescenceClass,        only: newCoalescenceData
    use preequilibriumDataClass, only: newPreequilibriumData

    implicit none
    class(GSM),            intent(inout) :: gsmObj
    type(GSMProjectile),   intent(in   ) :: projObj
    type(GSMTarget),       intent(in   ) :: targObj
    type(rxnSpecificData), intent(inout) :: rxnData

! ==============================================================================

    rxnData%numConstructionErrors = 0_int32

    ! Standard DCM Data:
    rxnData%sDCMTargetData = newStandardDCMData( &
         & targObj%numBaryons, &
         & targObj%numProtons, &
         & gsmObj%options%sDCMDataOpts, &
         & gsmObj%io%print )
    if ( .not.rxnData%sDCMTargetData%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) "The StandardDCMData object failed to construct."
       call gsmObj%io%print(1, 2, gsmObj%io%message)
       rxnData%numConstructionErrors = rxnData%numConstructionErrors + 1
    end if
    gsmObj%options%sDCMDataOpts = rxnData%sDCMTargetData%queryOptions()

    ! Coalescence Data:
    rxnData%coales = newCoalescenceData( &
         & projObj%kinEnergy, &
         & gsmObj%options%coalesDataOpts )
    if ( .not.rxnData%coales%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) "The CoalescenceData object failed to construct."
       call gsmObj%io%print(1, 2, gsmObj%io%message)
       rxnData%numConstructionErrors = rxnData%numConstructionErrors + 1
    end if
    gsmObj%options%coalesDataOpts = rxnData%coales%queryOptions()


    ! Preequilibrium Data:
    rxnData%preeq = newPreequilibriumData( &
         & projObj%numBaryons, &
         & projObj%numProtons, &
         & targObj%numBaryons, &
         & targObj%numProtons, &
         & (thousand * projObj%kinEnergy), &
         & gsmObj%genData%fissBarrier, &
         & gsmObj%io%print )
    if ( .not.rxnData%preeq%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) "The PreequilibriumData object failed to construct."
       call gsmObj%io%print(1, 2, gsmObj%io%message)
       rxnData%numConstructionErrors = rxnData%numConstructionErrors + 1       
    end if


    ! Check how many errors occurred:
    if ( rxnData%numConstructionErrors /= 0 ) then
       write(gsmObj%io%message, 1100) rxnData%numConstructionErrors
       call gsmObj%io%print(1, 2, gsmObj%io%message)
    end if

    return
! ==============================================================================
1000 format(A)
1100 format(i2, "  reaction-specific data objects failed to construct within GSM.")
! ==============================================================================
  end subroutine constructSpecificData
