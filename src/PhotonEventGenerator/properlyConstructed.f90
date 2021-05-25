
  function properlyConstructed( photonEG ) result(constructed)

! ====================================================================
!
! This function returns the logical flag stating if the photonEG
! object was constructed or not.
!
! ====================================================================

    implicit none
    class(PhotonEventGenerator), intent(in   ) :: photonEG
    logical                                    :: constructed

! ====================================================================

    constructed = photonEG%constructed
    return

! ====================================================================
  end function properlyConstructed
