
  function new_Results( progenyBnk ) result(results)

! ====================================================================
!
! Main constructor for the mDCMResults object (points to progeny array
! and sets internal-progeny related information)
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ====================================================================

    implicit none
    type(mDCMProgeny), intent(in   ), dimension(:), target :: progenyBnk
    type(mDCMResults) :: results

! ====================================================================

    ! Point to progeny bank:
    results%progenyBnk => progenyBnk

    ! Set progeny information:
    results%maxProgeny = size( results%progenyBnk )
    results%maxProgenyM1 = results%maxProgeny - 1

    ! Flag object as constructed:
    results%constructed = .TRUE.

    return
! ====================================================================
  end function new_Results
