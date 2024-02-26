! =============================================================================
!
!> \file
!> \brief  Contains the insist method
!> \author CMJ (XCP-3; LANL)
!
! Defines the "insist" contract to ensure a condition is met.
! Insist is ALWAYS on.
!
! =============================================================================
subroutine insist(cond, msg, fileName, line)
    use iso_fortran_env, only: int32
   implicit none
   logical,        intent(in   ) :: cond
   character(*),   intent(in   ) :: msg
   character(*),   intent(in   ), optional :: fileName
   integer(int32), intent(in   ), optional :: line

   character(:), allocatable :: insistMsg

   if (.not.(cond)) then
       insistMsg = "Insist failure: " // msg
       call stop_execution(insistMsg, fileName, line)
   end if

   return
end subroutine insist
