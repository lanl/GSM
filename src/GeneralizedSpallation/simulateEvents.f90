
  subroutine simulateEvents(gsmObj, gsmRxn, numDesiredEvents, maxEventAttempts)

! ====================================================================
!
! Simulates as many events as requested by the end-user/client
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, int64, real64
    use standardDCMClass, only: sDCMProjectile

    implicit none
    class(GSM),        intent(inout) :: gsmObj
    type(GSMReaction), intent(inout) :: gsmRxn
    integer(int64),    intent(in   ) :: numDesiredEvents
    integer(int64),    intent(in   ) :: maxEventAttempts

    integer(int64) :: numEventAttempts = 0_int64
    integer(int64), target :: eventNum = 0_int64
    type(sDCMProjectile), target :: sDCMProj

! ====================================================================

    ! Setup projectile:
    sDCMProj%numBaryons  = gsmRxn%results%initialProj%numBaryons
    sDCMProj%numProtons  = gsmRxn%results%initialProj%numProtons
    sDCMProj%kinEnergy   = gsmRxn%results%initialProj%kinEnergy
    sDCMProj%restMass    = gsmRxn%results%initialProj%restMass
    sDCMProj%decayNumber = gsmRxn%results%initialProj%decayNumber
    if ( gsmRxn%results%initialProj%particleFlag == nucleusProjFlag .or. &
       & gsmRxn%results%initialProj%particleFlag == pionProjFlag ) then
       sDCMProj%gammaFlag = 0_int32
    else
       sDCMProj%gammaFlag = 1_int32
    end if

    ! Reset values:
    eventNum = 0_int64

    ! Point reaction to interim results:
    gsmRxn%sDCMProj => sDCMProj
    gsmRxn%eventNum => eventNum

    !>>> PRIVATE GSMRXN MAY REQUIRE ADDITIONAL SETUP OF SUB-OBJECTS
    eventLoop: do numEventAttempts = 1, maxEventAttempts

       ! Exit the loop when the max. number of desired events has
       !    successfully been simulated:
       gsmRxn%eventAttempt = numEventAttempts
       if ( eventNum >= numDesiredEvents ) exit eventLoop


       ! Simulate an event until a valid event occurs:
       call gsmObj%eventLoop(gsmRxn)


       ! Tally results (no results with elastic events)
!       if ( results%numElasticEvents == 0 ) then
          ! Set information regarding fission:
!          if ( fusion > zro .and. gsmRxn%outData%idel > 0 ) nfis = nfis + 1
!          sfu = sfu + wf

          ! Tally event results:
!          call gsmObj%vlobd()
!          if ( ipisa /= 0 ) call gsmObj%pisaSpectra()
!       end if

    end do eventLoop
    gsmRxn%outData%limcc = eventNum + gsmRxn%outData%intel

    ! Unassociate pointers:
    gsmRxn%sDCMProj => NULL()
    gsmRxn%eventNum => NULL()


    return
! ====================================================================
  end subroutine simulateEvents
