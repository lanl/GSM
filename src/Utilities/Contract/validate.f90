! =============================================================================
!
!> \file
!> \brief  Contains the Validate subroutine
!> \author CMJ (XCP-3; LANL)
!
!
! =============================================================================
subroutine validate(cond, msg)
   implicit none
   logical,      intent(in   ) :: cond
   character(*), intent(in   ) :: msg

   Insist(cond, msg)
   return
end subroutine validate
