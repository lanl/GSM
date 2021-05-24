
  function new_FermiBreakUpResults( progenyBnk ) result(fbuResults)

! ====================================================================
!
! This procedure constructs the Fermi BreakUp results object that
! will be returned to clients upon simulation of FBU physics
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    implicit none
    type(fermiBreakUpProgeny), intent(in   ), dimension(:), target :: progenyBnk
    type(fermiBreakUpResults) :: fbuResults

! ====================================================================

    ! Assign progeny information given the particle bank:
    fbuResults%progenyBnk => progenyBnk
    fbuResults%maxProgeny = size( fbuResults%progenyBnk )

    ! Flag as properly constructed:
    fbuResults%constructed = .TRUE.

    return
! ====================================================================
  end function new_FermiBreakUpResults
