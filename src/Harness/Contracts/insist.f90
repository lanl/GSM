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
    use iso_fortran_env, only: error_unit
   implicit none
   logical,      intent(in   ) :: cond
   character(*), intent(in   ) :: msg

   if (.not.(cond)) then
       error stop "Insist failure: " // msg
   end if

   return
end subroutine insist
