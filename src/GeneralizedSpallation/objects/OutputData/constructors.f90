! =============================================================================
!
!> \file
!> \brief  Contains the OutputData object constructors
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

! =============================================================================
!
!> \fn    constructorMain
!> \brief Main OutputData constructor
!
! ARGUMENTS:
!> NONE
!
! =============================================================================
 function constructorMain( &
         & ) &
         & result(this)

     implicit none
     type(OutputData) :: this

! =============================================================================

     ! Interface to establish internal values of the object
     ! NOTE: This is where private data will be established/initially set.
     !       This will include counting variables, however there likely won't
     !       be many of these.
     this = OutputData()
     return
! =============================================================================
end function constructorMain


