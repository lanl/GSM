
  logical function properlyConstructed(mDCMObj)

! ====================================================================
!
! This function returns to the client a boolean flag that indicates
! if the Modified DCM class was properly constructed
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ====================================================================

    implicit none
    class(ModifiedDCM), intent(in) :: mDCMObj

! ====================================================================

    properlyConstructed = mDCMObj%constructed

    return
! ====================================================================
  end function properlyConstructed


  function queryOptions( mDCMObj ) result(options)

! ====================================================================
!
! Returns the 'options' object to the client
!
! ====================================================================

    implicit none
    class(ModifiedDCM), intent(in) :: mDCMObj
    type(mDCMOptions) :: options

! ====================================================================

    options = mDCMObj%options
    return
! ====================================================================
  end function queryOptions


  function queryRNG( mDCMObj ) result(rng)

! ====================================================================
!
! Returns the RNG procedure pointer to the client
!
! ====================================================================

    implicit none
    class(ModifiedDCM), intent(in) :: mDCMObj
    procedure(RANDOM), pointer :: rng => NULL()

! ====================================================================

    rng => mDCMObj%rang
    return
! ====================================================================
  end function queryRNG


  subroutine setOptions( mDCM, options )

! ====================================================================
!
! Sets and validates the options utilized by the object for its simulations
!
! ====================================================================

    implicit none
    class(ModifiedDCM), intent(inout) :: mDCM
    type(mDCMOptions),  intent(in   ) :: options

! ====================================================================

    mDCM%options = options
    call mDCM%validateOptions()

    return
! ====================================================================
  end subroutine setOptions


  subroutine validateOptions( mDCM )

! ====================================================================
!
! Validates all options contained by the mDCM object
!
! ====================================================================

    implicit none
    class(ModifiedDCM), intent(inout) :: mDCM

! ====================================================================

    ! Validate RM:
    if ( mDCM%options%rm <= 0.0_real64 ) then
       write(mDCM%io%message, 1000) mDCM%options%rm
       call mDCM%io%print(2, 3, mDCM%io%message)
       mDCM%options%rm = defaultRm
       write(mDCM%io%message, 1999) mDCM%options%rm
       call mDCM%io%print(2, 3, mDCM%io%message)
    end if

    ! Validate DELTA:
    if ( mDCM%options%delta <= 0.0_real64 ) then
       write(mDCM%io%message, 1100) mDCM%options%delta
       call mDCM%io%print(2, 3, mDCM%io%message)
       mDCM%options%delta = defaultDelta
       write(mDCM%io%message, 1999) mDCM%options%delta
       call mDCM%io%print(2, 3, mDCM%io%message)
    end if

    return
! ====================================================================
1000 format("An invalid RM option was utilized (", f7.4, ").")
1100 format("An invalid DELTA option was specified (", f7.4, ").")
1999 format("   Using the default value (", f7.4, ").")
! ====================================================================
  end subroutine validateOptions

