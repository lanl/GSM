
  subroutine initializeModifiedDCMData( &
       & clientPhotoFile, &
       & clientDecayFile, &
       & clientIO )

! ====================================================================
!
! Initializes all non-parameterized data from the mDCM data module
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    character(len=*),     intent(in   ), optional :: clientPhotoFile
    character(len=*),     intent(in   ), optional :: clientDecayFile
    procedure(IOHANDLER), intent(in   ), optional, pointer :: clientIO

! ====================================================================

    character(len=128) :: photoFile = defaultPhotoFile
    character(len=128) :: decayFile = defaultDecayFile

! ====================================================================

    mDCMDataInitialized = .TRUE.

    if ( present(clientIO ) ) then
       if ( associated(clientIO) ) then
          dataIO%print => clientIO
       else
          write(dataIO%message, 1000)
          call dataIO%print(2, 3, dataIO%message)
       end if
    end if

    ! Apply all client-specified arguments:
    if( present(clientPhotoFile) ) photoFile = clientPhotoFile
    if( present(clientDecayFile) ) decayFile = clientDecayFile

    ! Set the effective files:
    effectiveDecayFile  = decayFile

    ! Initialize all non-parameterized data here and beyond:
    call readPhotonData(photoFile)
    call initam()

    return
! ====================================================================
1000 format("The provided I/O procedure for the mDCM data is not ", &
          & "associated and will not be used.")
! ====================================================================
  end subroutine initializeModifiedDCMData
