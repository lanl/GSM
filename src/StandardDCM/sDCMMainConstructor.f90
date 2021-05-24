
  function sDCMMainConstructor(clientRNG, clientMolnix, &
       & clientOptions, clientIO) result(sDCM)

! ====================================================================
!
! Returns a class for the standard DCM
!
!
! Written by CMJ, XCP-3, 12/2018
!
! ====================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use standardDCMData, only: sDCMDataEstablished, &
         & initializeStandardDCMData
    use molnixClass, only: Molnix, newMolnix

    implicit none
    procedure(RANDOM),    intent(in   ), pointer  :: clientRNG
    type(Molnix),         intent(in   ), target   :: clientMolnix
    type(sDCMOptions),    intent(in   ), optional :: clientOptions
    procedure(IOHANDLER), intent(in   ), optional, pointer :: clientIO
    type(standardDCM)                             :: sDCM

    integer(int32) :: errorFlag = 0_int32

! ====================================================================

    ! Flag proper construction
    sDCM%constructed = .TRUE.


    ! Point to client's I/O handler if present
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          sDCM%io%print => clientIO
       else
          write(sDCM%io%message, 1000)
          call sDCM%io%print(2, 3, sDCM%io%message)
       end if
    end if


    ! Verify Standard DCM data was properly initialized (use defaults)
    if ( .not. sDCMDataEstablished ) then
       ! Note: warning for print is embedded in the called interface (don't include here)
       errorFlag = initializeStandardDCMData( clientIO = sDCM%io%print )
    end if


    ! Setup class options
    if ( present(clientOptions) ) then
       sDCM%options = clientOptions
    end if


    ! Point to client's RNG
    if ( associated(clientRNG) ) then
       sDCM%rang => clientRNG
    else
       write(sDCM%io%message, 1100)
       call sDCM%io%print(1, 2, sDCM%io%message)
    end if


    ! Point to the necessary physics
    sDCM%molnixE => clientMolnix
    if ( .not.sDCM%molnixE%properlyConstructed() ) then
       sDCM%molnixE = newMolnix( clientIO = sDCM%io%print )
    end if


    return
! ====================================================================
1000 format("The I/O procedure given to the sDCM object is not ", &
          & "associated and will not be used.")
1100 format("The RNG procedure given to the sDCM object is not ", &
          & "associated and will not be used.")
! ====================================================================
  end function sDCMMainConstructor
