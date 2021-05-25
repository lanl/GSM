
  subroutine initializeGSMData( &
       & clientRndmType, &
       & gammaFile, &
       & photoFile, decayFile, &
       & massFile, levelFile, shellFile, &
       & clientIO)

! ==============================================================================
!
! Initializes all GSM data
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32
    use randomNumberGenerator, only: numRNGenerators, standard_generator, &
         & RN_init_problem, io
    use fund_data, only: radius_rms
    use standardDCMData, only: &
         & defaultGammaFile
    use modifiedDCMData, only: &
         & defaultPhotoFile, defaultDecayFile
    use evaporationFissionData, only: &
         & defaultMassFile,  &
         & defaultLevelFile, &
         & defaultShellFile

    implicit none
    ! For evaporation data initialization:
    integer(int32),   intent(in   ), optional :: clientRndmType
    character(LEN=*), intent(in   ), optional :: gammaFile
    character(LEN=*), intent(in   ), optional :: photoFile
    character(LEN=*), intent(in   ), optional :: decayFile
    character(LEN=*), intent(in   ), optional :: massFile
    character(LEN=*), intent(in   ), optional :: levelFile
    character(LEN=*), intent(in   ), optional :: shellFile
    procedure(IOHANDLER), intent(in   ), optional, pointer :: clientIO

    ! Default values for all optional arguments:
    integer(int32)     :: rndmType      = defaultRndmType
    character(LEN=128) :: gammaFilename = defaultGammaFile
    character(LEN=128) :: photoFilename = defaultPhotoFile
    character(LEN=128) :: decayFilename = defaultDecayFile
    character(LEN=128) :: massFilename = defaultMassFile
    character(LEN=128) :: levelFilename = defaultLevelFile
    character(LEN=128) :: shellFilename = defaultShellFile

! ==============================================================================

    ! If data was already established successfully, exit:
    if ( gsmDataInitialized ) return

    ! Flag data as being established:
    gsmDataInitialized = .TRUE.

    ! Use client-provided I/O handler if appropriate:
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          dataIO%print => clientIO
       else
          write(dataIO%message, 1000)
          call dataIO%print(3, 3, dataIO%message)
       end if
    end if
    io%print => dataIO%print

    ! Use client specified file names/units if included:
    if ( present(clientRndmType) ) rndmType = clientRndmType
    if ( present(gammaFile) ) gammaFileName = gammaFile
    if ( present(photoFile) ) photoFileName = photoFile
    if ( present(decayFile) ) decayFileName = decayFile
    if ( present(massFile) ) massFileName = massFile
    if ( present(levelFile) ) levelFileName = levelFile
    if ( present(shellFile) ) shellFileName = shellFile


    ! Initialize Random Number Generator and point to it:
    if ( rndmType < 1 .or. rndmType > numRNGenerators ) then
       write(dataIO%message, 1100) rndmType
       call dataIO%print(1, 3, dataIO%message)
       rndmType = defaultRndmType
       write(dataIO%message, 1110) rndmType
       call dataIO%print(1, 3, dataIO%message)
    end if
    call RN_init_problem( rndmType, defaultSeed, defaultStride, &
         & defaultNumSeedAdvance, defaultPrintInfo)

    ! Set r_rms multipliers
    call radius_rms()

    ! Initialize sub-model data:
    call initializeSubModelData( &
         & gammaFileName, &
         & photoFileName, decayFileName, &
         & massFileName, levelFileName, shellFileName, &
         & gsmDataInitialized)
    

    return
! ==============================================================================
1000 format("The GSM data provided I/O handler is not associated and will ", &
          & "not be used.")
1100 format("The random number generator of type ", i1, " does not exit.")
1110 format("   Defaulting to the random number generator of type ", i1, ".")
! ==============================================================================
  end subroutine initializeGSMData
