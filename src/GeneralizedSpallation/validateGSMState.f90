
  subroutine validateGSMState(gsmObj, results)

! ============================================================
!
! Validates the state of the GSM and results objects
!
!
! Written by CMJ, XCP-3 (07/2019)
!
! ============================================================

    use generalizedSpallationData, only: gsmDataInitialized, &
         & initializeGSMData

    implicit none
    class(GSM),        intent(inout) :: gsmObj
    class(gsmResults), intent(inout) :: results

! ============================================================

    ! Set default state of the results object:
    results%simState = successfulSingleEvent

    ! Verify GSM is constructed:
    if(.not.gsmObj%constructed) then
       ! Assume construction, then attempt to resolve errors:
       gsmObj%constructed = .TRUE.
       results%simState = successfulSingleEvent
       write(gsmObj%io%message, 1000) &
            & "The GSM object is in an incomplete state."
       call gsmObj%io%print(3, 3, gsmObj%io%message)

       ! Try to initialize data:
       if(.not.gsmDataInitialized) then
          write(gsmObj%io%message, 1000) &
               & "Model data for GSM failed to initialize."
          call gsmObj%io%print(2, 2, gsmObj%io%message)
          write(gsmObj%io%message, 1000) &
               & "   Attempting to initialize..."
          call gsmObj%io%print(2, 4, gsmObj%io%message)
          call initializeGSMData()
          if(.not.gsmDataInitialized) then
             write(gsmObj%io%message, 1000) "      Failed."
             call gsmObj%io%print(2, 3, gsmObj%io%message)
             results%simState = noDataConstruction
             gsmObj%constructed = .FALSE.
          else
             write(gsmObj%io%message, 1000) "      Success."
             call gsmObj%io%print(2, 4, gsmObj%io%message)
          end if
       end if

       ! Try to initialize general objects:
       call gsmObj%constructGeneralModels()
       if(gsmObj%genModels%numConstructionErrors > 0 .or. &
            & gsmObj%genData%numConstructionErrors > 0) &
            & results%simState = noObjectConstruction
    end if

    ! Verify the results object was constructed:
    if(.not.results%constructed) then
       results%simState = noResultsConstruction
       write(gsmObj%io%message, 1000) &
            & "The results object was not properly established."
       call gsmObj%io%print(1, 1, gsmObj%io%message)
    end if

    ! Warn user if not able to validate the objects' states:
    if(results%simState /= successfulSingleEvent) then
       write(gsmObj%io%message, 1000) &
            & "   Cannot simulate spallation physics."
       call gsmObj%io%print(1, 1, gsmObj%io%message)
    end if

    return
! ============================================================
1000 format(A)
! ============================================================
  end subroutine validateGSMState
