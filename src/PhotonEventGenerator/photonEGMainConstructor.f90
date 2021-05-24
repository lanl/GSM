
  function photonEGMainConstructor ( thetaValues, numElements, &
       & clientIO ) result(photonEG)

! ======================================================================
!
! Constructs a photon event generator object (client I/O) and
! sets up data for the "qintxs" function (i.e. angle bins for interpolation)
!
!
! Written by CMJ, XCP-3, 12/2018
! Edited  by CMJ, XCP-3, 01/2019
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    real(real64),   intent(in   ) :: thetaValues(:)
    integer(int32), intent(in   ) :: numElements
    procedure(IOHANDLER), intent(in   ), pointer, optional :: clientIO
    type(PhotonEventGenerator)    :: photonEG

! ======================================================================
 
    ! Point to client's I/O if it is being utilized
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          photonEG%io%print => clientIO
       else
          write(photonEG%io%message, 1000)
          call photonEG%io%print(2, 3, photonEG%io%message)
          write(photonEG%io%message, 1010)
          call photonEG%io%print(2, 3, photonEG%io%message)
       end if
    endif


    ! Set up data for the cross section interpolation function
    call photonEG%setupQintsData( thetaValues, numElements )


    ! Flag that the class was constructed
    photonEG%constructed = .TRUE.

    return
! ======================================================================
1000 format("The message handling procedure provided to the photon ", &
          & "event generator object")
1010 format("   is not associated and will not be used.")
! ======================================================================
  end function photonEGMainConstructor
