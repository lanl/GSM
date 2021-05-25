
  function sDCMDataMainConstructor(aTarg, zTarg, clientOptions, clientIO) result(dataObj)

! ==============================================================================
!
! Returns a data class for the standard DCM
!
!
! Written by CMJ, XCP-3, 12/2018
!
! ==============================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use standardDCMDataParams, only: zro, one, ato3rd

    implicit none
    real(real64),   intent(in   ) :: aTarg   ! Number of baryons in the target
    real(real64),   intent(in   ) :: zTarg   ! Number of protons in the target
    type(sDCMDataOptions), intent(in   ), optional          :: clientOptions
    procedure(IOHANDLER),  intent(in   ), optional, pointer :: clientIO
    type(StandardDCMData) :: dataObj

    integer(int32) :: targetError = 0
    real(real64)    :: targetA, targetZ

! ==============================================================================

    ! State that object was constructed
    dataObj%objectConstructed = .TRUE.
    dataObj%constructed = properlyConstructed


    ! Set data object's I/O handling
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          dataObj%io%print => clientIO
       else
          write(dataObj%io%message, 1200)
          call dataObj%io%print(0, 2, dataObj%io%message)
       end if
    end if


    ! Check for physicality
    targetA = aTarg
    targetZ = zTarg
    if ( (targetA < 1 .or. targetA > 300) .or. (targetZ > targetA) ) then
       ! Unphysical target object, signal failed construction
       write(dataObj%io%message, 1000) nint(targetA), nint(targetZ)
       call dataObj%io%print(1, 2, dataObj%io%message)
       dataObj%constructed = dataObj%constructed + invalidTargetFlag

       ! Approximate the nucleus and approximate:
       if ( targetA <   1 ) targetA =   1
       if ( targetA > 300 ) targetA = 300
       if ( targetZ > targetA ) targetA = targetZ
       write(dataObj%io%message, 1050) nint(targetA), nint(targetZ)
       call dataObj%io%print(1, 2, dataObj%io%message)
    end if
    dataObj%target%numBaryons  = targetA
    dataObj%target%numProtons  = targetZ
    dataObj%target%aTargThrd = ato3rd( nint(targetA) )


    if ( present(clientOptions) ) then
       ! Use and validate all client-specified options:
       dataObj%options = clientOptions
       call dataObj%validateOptions()
    end if


    ! Construct target object now:
    targetError = dataObj%setupTarget()
    if ( targetError == targetInitFailed ) then
       ! An error occurred while setting up the target data
       write(dataObj%io%message, 1200) targetA, targetZ
       call dataObj%io%print(0, 2, dataObj%io%message)
       dataObj%constructed = dataObj%constructed + targetInitFailed
    end if


    if ( dataObj%constructed /= properlyConstructed ) then
       dataObj%objectConstructed = .FALSE.
    end if


    return
! ==============================================================================
1000 format("The target object cannot be constructed due to an unphysical ", &
          & "target mass number (A=", i3, ", Z=", i3, ").")
1050 format("   The target will be approximated (A=", i3, ", Z=", i3, ").")
1200 format("Target nucleus properties in the Standard DCM Data object could ", &
          & "not be established (A=", f5.1, ", Z=", f5.1, ").")
! ==============================================================================
  end function sDCMDataMainConstructor
