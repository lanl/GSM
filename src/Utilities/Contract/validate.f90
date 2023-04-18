! =============================================================================
!
!> \file
!> \brief  Contains the Validate subroutine
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================
subroutine validate(cond, msg, file, line)
   use iso_fortran_env, only: int32, error_unit

   implicit none
   logical,        intent(in   ) :: cond
   character(*),   intent(in   ) :: msg
   character(*),   intent(in   ) :: file
   integer(int32), intent(in   ) :: line

   character(len = 6) :: lineStr

   if (.not.(cond)) then
       write(lineStr, "(A)") line

       write(error_unit, "(A)") "Validation failed: " // trim(adjustl(msg))
       if (showLine) then
           write(error_unit, "(A)") "^^^ at " // &
               & file // ":" // trim(adjustl(lineStr))
       end if

       error stop "Execution stopped."
   end if

   return
end subroutine validate
