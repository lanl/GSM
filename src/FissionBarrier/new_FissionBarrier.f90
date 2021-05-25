
  function new_FissionBarrier(clientMolnix, clientOptions, clientIO) result(fbObj)

! ===================================================================================
!
! Constructor for the FissionBarrier class
!
! USE:
!    fissionBarrierObject = FissionBarrier(newMolnix, [r0m])
!
! REQUIRED ARGUMENTS: A Molnix class must be passed in to this class for fission
!     barriers to be properly calculated.
!
! OPTIONAL ARGUMENTS:
! (1) r0m may be specified to use a different radius parameter to obtain fission
!     barrier values.
!
!
! Written by CMJ, XCP-3, 8/2018
!
! ===================================================================================

    use molnixClass, only: Molnix, newMolnix

    implicit none
    type(Molnix),                intent(in   ), target            :: clientMolnix
    type(fissionBarrierOptions), intent(in   ), optional          :: clientOptions
    procedure(IOHANDLER),        intent(in   ), optional, pointer :: clientIO
    type(FissionBarrier)                                          :: fbObj

! ===================================================================================

    ! Indicate that the object was properly constructed:
    fbObj%constructed = .TRUE.


    ! Set where messages go:
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          fbObj%io%print => clientIO
       else
          write(fbObj%io%message, 1100)
          call fbObj%io%print(2, 3, fbObj%io%message)
       end if
    end if

    ! Use the information of the 'newMolnix' class
    fbObj%fbMolnix => clientMolnix
    if ( .not.fbObj%fbMolnix%properlyConstructed() ) then
       fbObj%fbMolnix = newMolnix( clientIO = fbObj%io%print )
    end if


    ! Use r0 for radius parameter
    if ( present(clientOptions) ) then

       ! Use client's options:
       fbObj%options = clientOptions

       ! Verify the options specified are valid:
       if ( fbObj%options%r0m < 0 ) then
          write(fbObj%io%message, 1000) fbObj%options%r0m
          call fbObj%io%print(2, 3, fbObj%io%message)
          write(fbObj%io%message, 1010) defaultR0m
          call fbObj%io%print(2, 3, fbObj%io%message)
          fbObj%options%r0m = defaultR0m
       end if
    end if

    return
! ===================================================================================
1000 format("An invalid radius multiplier (", f6.3, " [fm]) was specified in the ", &
          & "Fission Barrier object.")
1010 format("   A multiplier of ", f6.3, " [fm] will instead be used.")
1100 format("The provided I/O procedure for the Fission Barrier object is not ", &
          & "associated and will not be used.")
! ===================================================================================
  end function new_FissionBarrier
