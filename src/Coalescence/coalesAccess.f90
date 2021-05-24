
! ====================================================================
!
! This file provides the various querying methods by which users can
! query the coalescence object
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

  function properlyConstructedObject(coalesObj) result(constructed)

! ====================================================================
!
! Returns if the class is constructed or not
!
! ====================================================================

    implicit none
    class(Coalescence), intent(in   ) :: coalesObj
    logical :: constructed

! ====================================================================

    constructed = coalesObj%constructed

    return
! ====================================================================
  end function properlyConstructedObject


  function queryOptions(coalesObj) result(options)

! ====================================================================
!
! Returns the options of the class
!
! ====================================================================

    implicit none
    class(Coalescence), intent(in   ) :: coalesObj
    type(coalescenceOptions) :: options

! ====================================================================

    options = coalesObj%options

    return
! ====================================================================
  end function queryOptions


  function queryData(coalesObj) result(dataObj)

! ====================================================================
!
! Returns the data object of the class
!
! ====================================================================

    implicit none
    class(Coalescence), intent(in   ) :: coalesObj
    type(CoalescenceData) :: dataObj

! ====================================================================

    dataObj = coalesObj%data

    return
! ====================================================================
  end function 
