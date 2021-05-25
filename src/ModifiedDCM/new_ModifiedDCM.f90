
  function new_ModifiedDCM( clientRNG, &  ! Required arguments
       & clientOptions, clientIO ) &   ! Optional arguments
       & result( mDCM )

! ====================================================================
!
! Main constructor for the mDCMResults object (points to progeny array
! and sets internal-progeny related information)
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ====================================================================

    use modifiedDCMData, only: mDCMDataInitialized, &
         & initializeModifiedDCMData

    implicit none
    procedure(RANDOM),    intent(in   ), pointer  :: clientRNG
    type(mDCMOptions),    intent(in   ), optional :: clientOptions
    procedure(IOHANDLER), intent(in   ), optional, pointer :: clientIO
    type(ModifiedDCM) :: mDCM

! ====================================================================

    ! Flag object as constructed:
    mDCM%constructed = .TRUE.

    ! Setup I/O if provided:
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          mDCM%io%print => clientIO
       else
          write(mDCM%io%message, 1000)
          call mDCM%io%print(2, 3, mDCM%io%message)
       end if
    end if

    ! Check if model data was constructed:
    if ( .not. mDCMDataInitialized ) then
       call initializeModifiedDCMData( clientIO = mDCM%io%print)
       if ( .not. mDCMDataInitialized ) then
          write(mDCM%io%message, 2000)
          call mDCM%io%print(1, 2, mDCM%io%message)
          write(mDCM%io%message, 2010)
          call mDCM%io%print(1, 2, mDCM%io%message)
          mDCM%constructed = .FALSE.
       end if
    end if


    ! Validate and use the RNG:
    if ( associated(clientRNG) ) then
       mDCM%rang => clientRNG
    else
       write(mDCM%io%message, 1100)
       call mDCM%io%print(1, 2, mDCM%io%message)
       mDCM%constructed = .FALSE.
    end if


    ! Set options, if provided:
    if ( present(clientOptions) ) then
       call mDCM%setOptions( clientOptions )
    end if

    return
! ====================================================================
1000 format("The I/O procedure given to the mDCM object is not ", &
          & "associated and will not be used.")
1100 format("The RNG procedure given to the mDCM object is not ", &
          & "associated and will not be used.")
2000 format("The mDCM data could not be initialized during object ", &
          & "construction.")
2010 format("   Object will be constructed.")
! ====================================================================
  end function new_ModifiedDCM
