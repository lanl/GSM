
  subroutine simulateEvent(gsmObj, proj, targ, results)

! ====================================================================
!
! Simulates a single event in GSM
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, int64, real64
    use standardDCMClass, only: sDCMProjectile

    implicit none
    class(GSM),           intent(inout) :: gsmObj
    class(GSMProjectile), intent(inout) :: proj
    class(GSMTarget),     intent(inout) :: targ
    class(GSMResults),    intent(inout) :: results

    type(GSMReaction)       :: gsmRxn
    type(sDCMProjectile), target :: sDCMProj

! ====================================================================

    ! Verify state of the GSM object:
    call gsmObj%validateGSMState(results)
    if(results%simState /= successfulSingleEvent) return

    ! Validate nuclei:
    call gsmObj%formNuclei(proj, targ)
    call gsmObj%setMDCMReaction(proj, targ)
    results%initialProj = proj
    results%initialTarg = targ

    ! Setup the sDCM projectile:
    sDCMProj%numBaryons  = proj%numBaryons
    sDCMProj%numProtons  = proj%numBaryons
    sDCMProj%kinEnergy   = proj%kinEnergy
    sDCMProj%restMass    = proj%restMass
    sDCMProj%decayNumber = proj%decayNumber
    sDCMProj%gammaFlag   = proj%particleFlag

    ! Establish reaction-specific data:
    gsmRxn = gsmObj%setupReaction(proj, targ, results)
    gsmRxn%incFlag = gsmObj%inquireINCModel(proj, targ)
    gsmRxn%sDCMProj => sDCMProj

    ! Simulate an event:
    call gsmObj%eventLoop(gsmRxn)

    ! Unswap nuclei if applicable
    if (proj%system == antilabSystem) &
        & call gsmObj%swapNuclei(proj, targ)

    return
! ====================================================================
  end subroutine simulateEvent
