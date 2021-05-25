
  function initializeStandardDCMData( &
       & clientGammaFile, &
       & clientIO) result(errorFlag)

! ====================================================================
!
! Initializes data utilized by the standard Dubna Cascade Model (S-DCM)
!
!
! Written by CMJ, XCP-3, 12/2018
!
! ====================================================================

    use, intrinsic :: iso_fortran_env, only: int32

    implicit none
    character(LEN=*), intent(in   ), optional :: clientGammaFile
    procedure(IOHANDLER), pointer, intent(in   ), optional :: clientIO
    integer(int32) :: errorFlag

    character(LEN=64) :: gammaFile = defaultGammaFile

! ====================================================================

    ! Error flags for the associated data initialization routines
    integer(int32) :: qintsError = 0_int32
    integer(int32) :: gammaError = 0_int32

! ====================================================================

    errorFlag = sDCMDataSetup

    ! If client wants to control printing, change message handling procedure
    if ( present(clientIO) ) then
       dataIO => clientIO
    end if


    ! Intialize data for the "qints" function
    qintsError = initializeQintsData()
    if ( qintsError == sDCMQintError ) then
       write(message,1000) "The cross section interpolation data failed to initialize."
       call dataIO(0, 2, message)
       errorFlag = errorFlag + qintsError
    end if


    ! Initialize photon cross section data
    ! (use file name and/or unit specified)
    if ( present(clientGammaFile) ) gammaFile = trim(clientGammaFile)
    gammaError = initializeGammaData(gammaFile)
    if ( gammaError == sDCMGammaError ) then
       write(message,1000) "Photon cross section data was not properly initialized."
       call dataIO(0, 2, message)
       errorFlag = errorFlag + gammaError
    end if


    
    if ( errorFlag == sDCMDataSetup ) sDCMDataEstablished = .TRUE.

    return
! ====================================================================
1000 format(A)
! ====================================================================
  end function initializeStandardDCMData
