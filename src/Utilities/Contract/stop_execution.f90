! =============================================================================
!
!> \file
!> \fn
!> \brief  Contains the "stop_execution" method.
!> \author CMJ (XCP-3; LANL)
!
! Contains implementation to stop the program's current execution
!
! =============================================================================
subroutine stop_execution(msg, fileName, line)
    use iso_fortran_env, only: int32
   implicit none
   character(:),   intent(inout), allocatable :: msg
   character(*),   intent(in   ),              optional :: fileName
   integer(int32), intent(in   ),              optional :: line

   call append_location(msg, fileName, line)
   error stop msg;

   return
end subroutine stop_execution
