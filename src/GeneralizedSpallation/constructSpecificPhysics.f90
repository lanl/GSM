
  subroutine constructSpecificPhysics( gsmObj, rxnData, rxnPhysics )

! ==============================================================================
!
! This function constructs the reaction-specific physics data used by
! GSM in the generalPhysicsData data type
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use coalescenceClass,    only: newCoalescence
    use preequilibriumClass, only: newPreequilibrium


    implicit none
    class(GSM),              intent(inout) :: gsmObj
    type(rxnSpecificData),   intent(in   ) :: rxnData
    type(rxnSpecificModels), intent(  out) :: rxnPhysics

! ==============================================================================

    rxnPhysics%numConstructionErrors = 0_int32


    ! Coalescence:
    rxnPhysics%coales = newCoalescence( &
         & rxnData%coales, &
         & gsmObj%options%coalesOpts, &
         & gsmObj%io%print )
    if ( .not.rxnPhysics%coales%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) "The Coalescence object failed to construct."
       call gsmObj%io%print(1, 2, gsmObj%io%message)
       rxnPhysics%numConstructionErrors = rxnPhysics%numConstructionErrors + 1
    end if
    gsmObj%options%coalesOpts = rxnPhysics%coales%queryOptions()


    ! Preequilibrium:
    if ( .not.gsmObj%usePhotonEmission ) then
       rxnPhysics%preeq = newPreequilibrium( &
            & rxnData%preeq, &
            & gsmObj%rang, &
            & gsmObj%genModels%fbu, &
            & gsmObj%options%preeqOpts, &
            & clientIO = gsmObj%io%print )
    else
       rxnPhysics%preeq = newPreequilibrium( &
            & rxnData%preeq, &
            & gsmObj%rang, &
            & gsmObj%genModels%fbu, &
            & gsmObj%options%preeqOpts, &
            & preeqGammaCascade, &
            & gsmObj%io%print ) 
   end if
    if ( .not.rxnPhysics%preeq%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) "The Preequilibrium object failed to construct."
       call gsmObj%io%print(1, 2, gsmObj%io%message)
       rxnPhysics%numConstructionErrors = rxnPhysics%numConstructionErrors + 1
    end if
    gsmObj%options%preeqOpts = rxnPhysics%preeq%queryOptions()


    ! Check how many errors occurred:
    if ( rxnPhysics%numConstructionErrors /= 0 ) then
       write(gsmObj%io%message, 1100) rxnPhysics%numConstructionErrors
       call gsmObj%io%print(1, 2, gsmObj%io%message)
    end if

    return
! ==============================================================================
1000 format(A)
1100 format(i2, "  reaction-specific physics objects failed to construct within GSM.")
! ==============================================================================
  end subroutine constructSpecificPhysics
