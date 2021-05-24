
  function new_GSM(clientRNG, clientOptions, clientIO ) result(gsmObj)

! ======================================================================
!
! Constructor for the GSM Class
!
! USE:
!    myGSM = GSM(clientMolnix, clientEvapData, ...)
!
! REQUIRED ARGUMENTS:
! (1) A Molnix class is required to be passed in. It is used throughout
!     GSM and its various sub-models/classes
! (2) An evaporation data class is required to be passed in. It is passed
!     in to the Evaporation constructor and saves on computation time
!     significantly.
!
! OPTIONAL ARGUMENTS:
! (1)
!
!
! Written by CMJ, XCP-3, August 2018 - Decemeber 2018 (GSM Class Creation)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: real64
    use generalizedSpallationData, only: &
         & initializeGSMData, gsmDataInitialized, &
         & queryRNG

    implicit none
    procedure(RANDOM),        intent(in   ), pointer, optional :: clientRNG
    type(GSMOptions),         intent(in   ),          optional :: clientOptions
    procedure(IOHANDLER),     intent(in   ), pointer, optional :: clientIO
    type(GSM) :: gsmObj

! ======================================================================

    ! Flag object as being constructed:
    gsmObj%constructed = .TRUE.


    ! Establish GSM data:
    if ( .not. gsmDataInitialized ) then
       ! Warn user that an attempt to create data will be made:
       write(gsmObj%io%message, 1000) "The GSM object's data was not established."
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       write(gsmObj%io%message, 1000) "   Establishing data..."
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       call initializeGSMData(clientIO=gsmObj%io%print)
       ! Flag data initialization failure:
       if ( .not. gsmDataInitialized ) then
          write(gsmObj%io%message, 1000) "   Unable to establish GSM data."
          call gsmObj%io%print(1, 2, gsmObj%io%message)
          write(gsmObj%io%message, 1000) "   Cannot simulate spallation physics."
          call gsmObj%io%print(1, 2, gsmObj%io%message)
          gsmObj%constructed = .FALSE.
       else
          write(gsmObj%io%message, 1000) "   Done."
          call gsmObj%io%print(2, 3, gsmObj%io%message)
       end if
    end if


    ! Setup I/O Handler FIRST (all messages will be printed as desired):
    if ( present(clientIO) ) then
       call gsmObj%updateMessageHandler( clientIO )
    end if


    ! Point to client's various procedure pointer(s):
    gsmObj%useDefaultRNG = .TRUE.
    gsmObj%rang => queryRNG()
    if ( present(clientRNG) ) then
       call gsmObj%updateRNG( clientRNG )
    end if


    ! Set simulation options (if client specified)
    if ( present(clientOptions) ) then
       call gsmObj%updateOptions( clientOptions )
    end if


    ! Setup photon emission for both evaporation/preequilibrium physics:
    gsmObj%usePhotonEmission = useGammaCascade

    ! Construct general physics objects:
    call gsmObj%constructGeneralModels()

    return
! ======================================================================
1000 format(A)
! ======================================================================
  end function new_GSM
