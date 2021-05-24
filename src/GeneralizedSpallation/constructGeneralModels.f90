
  subroutine constructGeneralModels(gsmObj)

! ======================================================================
!
! Constructs the general data objects and the general physics models
! used by the GSM object
!
!
! Written by CMJ, XCP-3, August 2018 - Decemeber 2018 (GSM Class Creation)
!
! ======================================================================

    implicit none
    class(GSM), intent(inout) :: gsmObj

! ======================================================================

    ! Construct the various general objects used by GSM:
    call gsmObj%constructGeneralData()
    call gsmObj%constructGeneralPhysics()

    ! Print failure to construct if any models didn't construct properly:
    if ( gsmObj%genModels%numConstructionErrors > 0 .or. &
         & gsmObj%genData%numConstructionErrors > 0 ) then
       write(gsmObj%io%message, 1999)
       call gsmObj%io%print(1, 1, gsmObj%io%message)
       gsmObj%constructed = .FALSE.
    end if

    return
! ======================================================================
1999 format("   Unable to construct GSM object.")
! ======================================================================
  end subroutine constructGeneralModels
