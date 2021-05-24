
  function setupReaction(gsmObj, proj, targ, &
       & results, output) result(reaction)

! ====================================================================
!
! This procedure constructs all reaction-specific objects needed by
! the sub-models and setups up all other misc. quantities for a
! specific reaction
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    use OutputDataMod, only: OutputData, newOutputData

    implicit none
    class(GSM),               intent(inout) :: gsmObj
    class(GSMProjectile),     intent(in   ) :: proj
    class(GSMTarget),         intent(in   ) :: targ
    class(GSMResults),        intent(inout), target :: results
    class(GSMOutput),         intent(inout), target, optional :: output
    type(gsmReaction) :: reaction

! ====================================================================

    ! Establish reaction object:
    reaction%results     => results
    if (present(output)) then
       reaction%output => output
    end if
    if (.not.associated(reaction%outData)) allocate(reaction%outData)
    reaction%outData = newOutputData()
    reaction%constructed = .TRUE.

    ! Construct physics data and models:
    call gsmObj%constructSpecificData(proj, targ, reaction%data)
    call gsmObj%constructSpecificPhysics(reaction%data, &
         & reaction%models)

     ! Obtain geometric cross section for elastic/inelastic estimation
     reaction%outData%sigom = reaction%data%sDCMTargetData%geomCrossSection()


    ! Check for errors regarding the setup:
    if ( reaction%data%numConstructionErrors > 0 .or. &
         & reaction%models%numConstructionErrors > 0 ) then
       ! Warn user:
       write(gsmObj%io%message, 1999)
       call gsmObj%io%print(1, 1, gsmObj%io%message)
       ! Flag bad construction:
       reaction%constructed = .FALSE.
       results%simState = noObjectConstruction
    end if

    return
! ====================================================================
1999 format("   Unable to construct reaction-specific objects.")
! ====================================================================
  end function setupReaction
