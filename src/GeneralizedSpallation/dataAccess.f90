
  logical function properlyConstructed( gsmObj )

! ====================================================================
!
! This function returns to the client a boolean flag that indicates
! if the GSM class was properly constructed
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ====================================================================

    implicit none
    class(GSM), intent(in   ) :: gsmObj

! ====================================================================

    properlyConstructed = gsmObj%constructed

    return
! ====================================================================
  end function properlyConstructed

  function reactionProperlyConstructed( gsmRxn ) result(constructed)

! ====================================================================
!
! Returns the construction flag of the gsmReaction object (indicates
! if errors occurred during setup)
!
! ====================================================================

    implicit none
    class(GSMReaction), intent(in   ) :: gsmRxn
    logical :: constructed

! ====================================================================

    constructed = gsmRxn%constructed

    return
! ====================================================================
  end function reactionProperlyConstructed


  function queryOptions( gsmObj ) result(options)

! ====================================================================
!
! Returns the GSM options to the client
!
! ====================================================================

    implicit none
    class(GSM), intent(in   ) :: gsmObj
    type(GSMOptions) :: options

! ====================================================================

    options = gsmObj%options
    return
! ====================================================================
  end function queryOptions


  function queryRNG( gsmObj ) result(rng)

! ====================================================================
!
! Returns the GSM RNG to the client
!
! ====================================================================

    implicit none
    class(GSM), intent(in   ) :: gsmObj
    procedure(RANDOM), pointer :: rng => NULL()

! ====================================================================

    rng => gsmObj%rang
    return
! ====================================================================
  end function queryRNG


  function queryPhotoEmissionUse( gsmObj ) result(useEmission)

! ====================================================================
!
! Returns the logic flag for whether or not GSM is using photon emission
!
! ====================================================================

    implicit none
    class(GSM), intent(in   ) :: gsmObj
    logical :: useEmission

! ====================================================================

    useEmission = gsmObj%usePhotonEmission
    return
! ====================================================================
  end function queryPhotoEmissionUse


  subroutine updateOptions( gsmObj, newOptions )

! ====================================================================
!
! Updates the options utilized by the GSM object
!
! ====================================================================

    implicit none
    class(GSM), intent(inout) :: gsmObj
    type(GSMOptions), intent(in   ) :: newOptions

! ====================================================================

    gsmObj%options = newOptions
    call gsmObj%validateOptions()

    return
! ====================================================================
  end subroutine updateOptions


  subroutine updateRNG( gsmObj, newRNG )

! ====================================================================
!
! Updates the RNG procedure pointer used by GSM
!
! ====================================================================

    implicit none
    class(GSM),        intent(inout) :: gsmObj
    procedure(RANDOM), intent(in   ), pointer :: newRNG

! ====================================================================

    ! Update client RNG, if it points to something:
    if ( associated(newRNG) ) then
       gsmObj%useDefaultRNG = .FALSE.
       gsmObj%rang => newRNG
    else
       write(gsmObj%io%message, 1200)
       call gsmObj%io%print(2, 3, gsmObj%io%message)
    end if

    return
! ====================================================================
1200 format("The client-provided RNG for GSM is not associated and ", &
          & "will not be used.")
! ====================================================================
  end subroutine updateRNG


  subroutine updateMessageHandler(gsmObj, newIO)

! ====================================================================
!
! Updates the message handling procedure pointer used by GSM
!
! ====================================================================

    implicit none
    class(GSM),           intent(inout) :: gsmObj
    procedure(IOHANDLER), intent(in   ), pointer :: newIO

! ====================================================================

    ! Update the IO handler IF it points to something:
    if ( associated(newIO) ) then
       gsmObj%io%print => newIO
    else
       write(gsmObj%io%message, 1100)
       call gsmObj%io%print(2, 3, gsmObj%io%message)
    end if

    return
! ====================================================================
1100 format("The client-provided message handler for GSM is not ", &
          & "associated and will not be used.")
! ====================================================================
  end subroutine updateMessageHandler
