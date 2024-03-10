
  function mDCMDataMainConstructor(clientOptions, clientIO) result(dataObj)

! ==============================================================================
!
! Returns a data class for the modified DCM
!
!
! Written by CMJ, XCP-3, 02/2024
!
! ==============================================================================

    implicit none
    class(mDCMDataOptions), intent(in   ), optional          :: clientOptions
    procedure(IOHANDLER),   intent(in   ), optional, pointer :: clientIO
    type(ModifiedDCMData) :: dataObj

! ==============================================================================

    ! State that object was constructed
    dataObj%objectConstructed = .TRUE.
    dataObj%constructed = properlyConstructed


    ! Set data object's I/O handling
    if ( present(clientIO) .and. associated(clientIO) ) then
       dataObj%io%print => clientIO
    end if


    ! Check for physicality of parameters
    ! <nothing yet>


    ! Set options if they were specified
    if ( present(clientOptions) ) then
       ! Use and validate all client-specified options:
       dataObj%options = clientOptions
       call dataObj%validateOptions()
    end if


    ! If construction in error, flag as not constructed
    if ( dataObj%constructed /= properlyConstructed ) then
       dataObj%objectConstructed = .FALSE.
    end if


    return
! ==============================================================================
1000 format("The target object cannot be constructed due to an unphysical ", &
          & "target mass number (A=", i3, ", Z=", i3, ").")
1050 format("   The target will be approximated (A=", i3, ", Z=", i3, ").")
1200 format("Target nucleus properties in the Modified DCM Data object could ", &
          & "not be established (A=", f5.1, ", Z=", f5.1, ").")
! ==============================================================================
  end function mDCMDataMainConstructor
