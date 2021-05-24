
  subroutine initializeSubModelData( &
       & gFile, &                 ! for SDCM Data
       & pFile, dFile, &          ! for MDCM Data
       & mFile, lFile, sFile, &   ! For Evaporation data
       & dataInitialized)

! ====================================================================
!
! Initializes all sub-model data variables
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32
    use standardDCMData,        only: initializeStandardDCMData, &
         & sDCMDataEstablished
    use modifiedDCMData,        only: initializeModifiedDCMData, &
         & mDCMDataInitialized
    use fermiBreakUpClass,      only: initializeFermiData
    use preequilibriumClass,    only: preequilibriumInit
    use evaporationFissionData, only: initializeEvaporationData, &
         & evaporationDataEstablished

    implicit none
    character(LEN=128), intent(in   ) :: gFile   ! Photon file name
    character(LEN=128), intent(in   ) :: pFile   ! Photon file name (mDCM)
    character(LEN=128), intent(in   ) :: dFile   ! Decay file name (mDCM)
    character(LEN=128), intent(in   ) :: mFile   ! Mass file name
    character(LEN=128), intent(in   ) :: lFile   ! Level file name
    character(LEN=128), intent(in   ) :: sFile   ! Shell file name
    logical,            intent(inout) :: dataInitialized

    integer(int32) :: errorFlag

! ====================================================================

    !>>> ASSIGN EACH CALL BELOW AN OPENMP TASK TO READ IN PHYSICS DATA

    ! sDCM Data:
    errorFlag = initializeStandardDCMData(gFile, dataIO%print)
    if ( .not.sDCMDataEstablished ) then
       ! Failure:
       write(dataIO%message, 1000) "Standard"
       call dataIO%print(2, 2, dataIO%message)
       dataInitialized = .FALSE.
    end if


    ! mDCM Data:
    call initializeModifiedDCMData(pFile, dFile, dataIO%print)
    if ( .not.mDCMDataInitialized ) then
       ! Failure:
       write(dataIO%message, 1000) "Modified"
       call dataIO%print(2, 2, dataIO%message)
       dataInitialized = .FALSE.
    end if


    ! Fermi Break-Up Data:
    call initializeFermiData()


    ! Preequilibrium Data:
    call preequilibriumInit()


    ! Evaporation Data:
    errorFlag = initializeEvaporationData(mFile, lFile, &
         & sFile)
    if ( .not.evaporationDataEstablished ) then
       ! Failure:
       write(dataIO%message, 1300)
       call dataIO%print(2, 2, dataIO%message)
       dataInitialized = .FALSE.
    end if


    return
! ====================================================================
1000 format("The ", A, " DCM models' data parameterizations failed ", &
          & "to initialize.")
1300 format("The Evaporation model's data parameterizations failed ", &
          & "to initialize.")
! ====================================================================
  end subroutine initializeSubModelData
