! =============================================================================
!
!> \file
!> \fn
!> \brief  Contains the require function for pre-condition checks
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================
subroutine require(cond)
   implicit none
   logical, intent(in) :: cond

!   Require(cond)
   if (.not. (cond)) then
      print *, &
           & 'Pre-condition check failed in ', &
           & __FILE__, ', line ', __LINE__
      error stop 'Program execution halted.'; 
   end if

   return
end subroutine require
