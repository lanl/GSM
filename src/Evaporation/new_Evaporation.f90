
  function new_Evaporation(data, clientRNG, clientMolnix, clientFBU, clientOptions, &
       & photonEmission, clientIO) result(evapObj)

! ===================================================================================
!
! Constructor for the 'Evaporation' class.
!
! USE:
!    evaporationObject = Evaporation(fragmentBnk, evaporationData, [options], [clientIO] )
!
!
! REQUIRED ARGUMENTS: Constructor must be given the fragment bank array for the particle
! bank data structure (evaporationFragment) to store the fragment data in, AND data
! that used for the simulation, but that is based on a calculation option for the
! level density. The "errorFlag" variable will warn the user if construction of the object
! has failed.
!
! OPTIONAL ARGUMENTS:
! (1) Clients may include the 'evaporationOptions' type if desired. Doing so allows
!     client programs to control how the evaporation model simulates and what
!     parameterizations to use.
! (2) The 'clientIO' variable allows clients to control how the evaporaiton model
!     outputs messages. The evaporation model makes calls to this procedure whenever
!     a message of any importance is created, giving the procedure a flag stating
!     how important the message is and a flag for the type of message (comment, error, etc.)
!
!
! Written by CMJ, XCP-3, 8/2018 (Class creation)
!
! ===================================================================================

    use, intrinsic :: iso_fortran_env, only: int32
    use evaporationFissionData, only: evaporationDataEstablished, &
         & initializeEvaporationData
    use molnixClass, only: Molnix, newMolnix
    use fermiBreakUpClass, only: FermiBreakup, newFermiBreakUp, fermiBreakUpOptions
    use evaporationDataClass,   only: EvaporationData, newEvaporationData

    implicit none
    type(EvaporationData),     intent(in   ), target   :: data
    procedure(RANDOM),         intent(in   ), pointer  :: clientRNG
    type(Molnix),              intent(in   ), target   :: clientMolnix
    type(FermiBreakUp),        intent(in   ), optional :: clientFBU
    type(evaporationOptions),  intent(in   ), optional :: clientOptions
    procedure(PHOTOEMISSION),  intent(in   ), optional, pointer :: photonEmission
    procedure(IOHANDLER),      intent(in   ), optional, pointer :: clientIO
    type(Evaporation)                                  :: evapObj

    integer(int32) :: dummyFlag = 0
    type(fermiBreakUpOptions) :: fbuOptions

! ===================================================================================

    ! Object has been successfully constructed
    evapObj%constructed = .true.


    ! Set where output goes (if specified)
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          evapObj%io%print => clientIO
       else
          write(evapObj%io%message, 1000)
          call evapObj%io%print(2, 3, evapObj%io%message)
       end if
    end if


    ! Ensure data is establish
    if ( .not. evaporationDataEstablished ) then
       ! Attempt to initialze evaporation data using the module's default values:
       ! Print to user the attempt...
       write(evapObj%io%message, 4000)
       call evapObj%io%print(1, 2, evapObj%io%message)
       dummyFlag = initializeEvaporationData()

       if ( evaporationDataEstablished ) then
          ! Attempt succeeded - print comment to user
          write(evapObj%io%message, 4010)
          call evapObj%io%print(1, 3, evapObj%io%message)
       else
          ! Failed to establish the required data. Warn user/client and return
          write(evapObj%io%message, 4020)
          call evapObj%io%print(1, 2, evapObj%io%message)
          evapObj%constructed = .FALSE.
       end if
    end if


    ! Point to necessary data and ensure its validity:
    evapObj%data => data
    evapObj%evapMolnix => clientMolnix
    if ( .not.evapObj%evapMolnix%properlyConstructed() ) then
       evapObj%evapMolnix = newMolnix( clientIO = evapObj%io%print )
    end if
    if ( .not.evapObj%data%properlyConstructed() ) then
       evapObj%data = newEvaporationData( clientIO = evapObj%io%print )
    end if
    if ( associated(clientRNG) ) then
       evapObj%rang => clientRNG
    else
       write(evapObj%io%message, 1100)
       call evapObj%io%print(2, 3, evapObj%io%message)
       evapObj%constructed = .FALSE.
    end if


    ! Set Fermi Break-Up object:
    if ( present(clientFBU) ) then
       evapObj%fbuObj = clientFBU
       if ( .not. evapObj%fbuObj%properlyConstructed() ) then
          evapObj%fbuObj = newFermiBreakUp( evapObj%rang, clientFBU%queryOptions(), &
               & evapObj%io%print )
       end if
    else
       evapObj%fbuObj = newFermiBreakUp( evapObj%rang, fbuOptions, evapObj%io%print )
    end if


    ! Setup photon emission procedure if specified:
    if ( present(photonEmission) ) then
       if ( associated(photonEmission) ) then
          evapObj%usePhotoEmission = .TRUE.
          evapObj%photonEmission => photonEmission
       else
          write(evapObj%io%message, 1200)
          call evapObj%io%print(2, 3, evapObj%io%message)
       end if
    end if


    ! Set object simulation behavior/options
    if ( present(clientOptions) ) then
       evapObj%options = clientOptions
       call evapObj%validateOptions()
    end if
    ! Setup the class variable(s) used in various procedures w/ these options:
    if ( evapObj%options%inverseParameter == 10 ) then
       evapObj%rrVcoul = defaultRRVcoul
    else if ( evapObj%options%inverseParameter > 0 ) then
       evapObj%rrVcoul = evapObj%options%inverseParameter
    end if


    return
! ===================================================================================
1000 format("The I/O procedure given to the evaporation object is not associated ", &
          & "and will not be used.")
1100 format("The RNG procedure given to the evaporation object is not associated ", &
          & "and will not be used.")
1200 format("The photon emission procedure given to the evaporation object is not associated ", &
          & "and will not be used.")
4000 format("Evaporation data was not established. An attempt will be made...")
4010 format("   Success.")
4020 format("   Failure. Be suspect of results.")
! ===================================================================================
  end function new_Evaporation
