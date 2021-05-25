
  function new_StandardDCMResults(clientProgeny) result(results)

! ==============================================================================
!
! This function returns a standard DCM results data type (stores progeny array
! into results type)
!
!
! Written by CMJ, XCP-3, 12/2018
!
! ==============================================================================

    implicit none
    type(sDCMProgeny), intent(in   ), dimension(:), target :: clientProgeny
    type(StandardDCMResults)                               :: results

! ==============================================================================

    ! Point to the desired progeny array
    results%progenyBnk => clientProgeny


    ! Set progeny information
    results%maxProgeny = size( results%progenyBnk )
    results%maxProgenyM3 = results%maxProgeny - 3 

    ! Flag that results object was properly initialized
    results%initialized = .TRUE.

    return
! ==============================================================================
  end function new_StandardDCMResults
