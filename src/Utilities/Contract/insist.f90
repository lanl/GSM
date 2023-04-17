! =============================================================================
!
!> \file
!> \brief  Contains the Contracts module
!> \author CMJ (XCP-3; LANL)
!
!> The Contracts module contains the various functions and macros for the
!> design-by-contract specifications. This module relies heavily on the Fortran
!> preprocessor.
!
! =============================================================================
subroutine insist(cond, msg)
   implicit none
   logical,      intent(in   ) :: cond
   character(*), intent(in   ) :: msg

   !Insist(cond, msg)
   return
end subroutine insist
