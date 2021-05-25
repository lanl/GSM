
  function new_FermiBreakUp(clientRNG, clientOptions, clientIO) &
       & result(fbuObj)

! ====================================================================
!
! This procedure is the main constructor for the Fermi BreakUp object.
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    use fermiBreakUpData, only: fermiDataInitialized

    implicit none
    procedure(RANDOM),         intent(in   ), pointer  :: clientRNG
    type(fermiBreakUpOptions), intent(in   ), optional :: clientOptions
    procedure(IOHANDLER),      intent(in   ), optional, pointer :: clientIO
    type(FermiBreakUp) :: fbuObj

! ====================================================================

    ! Flag object as being properly constructed
    fbuObj%constructed = .TRUE.


    ! Assign client's I/O if present and valid:
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          fbuObj%io%print => clientIO
       else
          write(fbuObj%io%message, 1100)
          call fbuObj%io%print(2, 3, fbuObj%io%message)
       end if
    end if


    ! Ensure data was established:
    if ( .not. fermiDataInitialized ) then
       write(fbuObj%io%message, 1000)
       call fbuObj%io%print(2, 3, fbuObj%io%message)
       call initializeFermiData()
       if ( fermiDataInitialized ) then
          write(fbuObj%io%message, 1010)
          call fbuObj%io%print(2, 3, fbuObj%io%message)
       else
          write(fbuObj%io%message, 1020)
          call fbuObj%io%print(2, 2, fbuObj%io%message)
          fbuObj%constructed = .FALSE.
       end if
    end if


    ! Assign RNG:
    fbuObj%rang => clientRNG
    if ( .not.associated(clientRNG) ) then
       write(fbuObj%io%message, 1200)
       call fbuObj%io%print(1, 1, fbuObj%io%message)
       fbuObj%constructed = .FALSE.
    end if


    ! Use and validate client options if present:
    if ( present(clientOptions) ) then
       fbuObj%options = clientOptions
       call fbuObj%validateOptions()
    end if


    return
! ====================================================================
1000 format("The Fermi Break-Up data was not initialized prior to ", &
          & "object construction.")
1010 format("   Data initialized.")
1020 format("   Data initialization failed. Unable to simulate Fermi ", &
          & "Break-Up physics.")
1100 format("The I/O handler given by the client is not associated ", &
          & "and will not be used.")
1200 format("The random number generator procedure given by the ", &
          & "client is not associated and will not be used.")
! ====================================================================
  end function new_FermiBreakUp
