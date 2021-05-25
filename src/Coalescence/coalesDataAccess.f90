
! ====================================================================
!
! This file provides the various querying methods by which users can
! query the coalescence data object
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

  function properlyConstructedData(dataObj) result(constructed)

! ====================================================================
!
! Returns if the class is constructed or not
! NOTE: The class can NOT be constructed and still operate correctly.
!
! ====================================================================

    implicit none
    class(CoalescenceData), intent(in   ) :: dataObj
    logical :: constructed

! ====================================================================

    constructed = dataObj%constructed

    return
! ====================================================================
  end function properlyConstructedData


  function queryDataOptions( dataObj ) result(options)

! ====================================================================
!
! Returns the options contained by the data object
!
! ====================================================================

    implicit none
    class(CoalescenceData), intent(in   ) :: dataObj
    type(coalescenceDataOptions) :: options

! ====================================================================

    options = dataObj%options
    return

! ====================================================================
  end function queryDataOptions


  function coalesRadiiDeut(dataObj) result(coalRadius)

! ====================================================================
!
! Returns the deuteron coalescnece radius
!
! ====================================================================

    implicit none
    class(CoalescenceData), intent(in   ) :: dataObj
    real(real64) :: coalRadius

! ====================================================================

    coalRadius = dataObj%radii%coalesRadiiDeut

    return
! ====================================================================
  end function coalesRadiiDeut


  function coalesRadiiTrit(dataObj) result(coalRadius)

! ====================================================================
!
! Returns the tritium coalescnece radius
!
! ====================================================================

    implicit none
    class(CoalescenceData), intent(in   ) :: dataObj
    real(real64) :: coalRadius

! ====================================================================

    coalRadius = dataObj%radii%coalesRadiiTrit

    return
! ====================================================================
  end function coalesRadiiTrit


  function coalesRadiiAlpha(dataObj) result(coalRadius)

! ====================================================================
!
! Returns the alpha (He-4) coalescnece radius
!
! ====================================================================

    implicit none
    class(CoalescenceData), intent(in   ) :: dataObj
    real(real64) :: coalRadius

! ====================================================================

    coalRadius = dataObj%radii%coalesRadiiAlpha

    return
! ====================================================================
  end function coalesRadiiAlpha


  function coalesRadiiLFrag(dataObj) result(coalRadius)

! ====================================================================
!
! Returns the light fragment coalescnece radius
!
! ====================================================================

    implicit none
    class(CoalescenceData), intent(in   ) :: dataObj
    real(real64) :: coalRadius

! ====================================================================

    coalRadius = dataObj%radii%coalesRadiiLFrag

    return
! ====================================================================
  end function coalesRadiiLFrag


  function coalesRadiiDeutSqrd(dataObj) result(coalRadius)

! ====================================================================
!
! Returns the deuteron coalescnece radius
!
! ====================================================================

    implicit none
    class(CoalescenceData), intent(in   ) :: dataObj
    real(real64) :: coalRadius

! ====================================================================

    coalRadius = dataObj%radii%coalesRadiiDeutSqrd

    return
! ====================================================================
  end function coalesRadiiDeutSqrd


  function coalesRadiiTritSqrd(dataObj) result(coalRadius)

! ====================================================================
!
! Returns the tritium coalescnece radius
!
! ====================================================================

    implicit none
    class(CoalescenceData), intent(in   ) :: dataObj
    real(real64) :: coalRadius

! ====================================================================

    coalRadius = dataObj%radii%coalesRadiiTritSqrd

    return
! ====================================================================
  end function coalesRadiiTritSqrd


  function coalesRadiiAlphaSqrd(dataObj) result(coalRadius)

! ====================================================================
!
! Returns the alpha (He-4) coalescnece radius
!
! ====================================================================

    implicit none
    class(CoalescenceData), intent(in   ) :: dataObj
    real(real64) :: coalRadius

! ====================================================================

    coalRadius = dataObj%radii%coalesRadiiAlphaSqrd

    return
! ====================================================================
  end function coalesRadiiAlphaSqrd


  function coalesRadiiLFragSqrd(dataObj) result(coalRadius)

! ====================================================================
!
! Returns the light fragment coalescnece radius
!
! ====================================================================

    implicit none
    class(CoalescenceData), intent(in   ) :: dataObj
    real(real64) :: coalRadius

! ====================================================================

    coalRadius = dataObj%radii%coalesRadiiLFragSqrd

    return
! ====================================================================
  end function coalesRadiiLFragSqrd
