
  function properlyConstructedObject(dataObj) result(constructed)

! ==============================================================================
!
! Returns flag stating if the data object is constructed or not
!
! ==============================================================================

    implicit none
    class(ModifiedDCMData), intent(in   ) :: dataObj
    logical :: constructed

! ==============================================================================

    constructed = dataObj%objectConstructed

    return
! ==============================================================================
  end function properlyConstructedObject


  function queryOptions( dataObj ) result(options)

! ==============================================================================
!
! Returns the options contained by the mDCM Data object
!
! ==============================================================================

    implicit none
    class(ModifiedDCMData), intent(in   ) :: dataObj
    type(mDCMDataOptions) :: options

! ==============================================================================

    options = dataObj%options
    return
! ==============================================================================
  end function queryOptions
