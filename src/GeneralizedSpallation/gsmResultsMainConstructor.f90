
  function gsmResultsMainConstructor(clientProgeny) result(results)

! ==============================================================================
!
! This function returns a results object to the client for GSM to use (client
! can track as they desire this way)
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    implicit none
    type(GSMProgeny), intent(in   ), dimension(:), target :: clientProgeny
    type(GSMResults)                                      :: results

! ==============================================================================

    ! Point to particle bank client wants to use
    results%progenyBnk => clientProgeny


    ! Determine maximum allowed number of progeny in the bank
    results%maxProgeny = size( results%progenyBnk )
    results%maxProgenyM1 = results%maxProgeny - 1

    ! Flag that results object was properly initialized
    results%constructed = .TRUE.
    results%simState = successfulSingleEvent

    return
! ==============================================================================
  end function gsmResultsMainConstructor
