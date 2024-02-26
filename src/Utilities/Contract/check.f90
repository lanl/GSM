! =============================================================================
!
!> \file
!> \fn
!> \brief  Contains the check function for intra-scope checks
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================
subroutine check(cond, fileName, line)
    use iso_fortran_env, only: int32
   implicit none
   logical,        intent(in   ) :: cond
   character(*),   intent(in   ), optional :: fileName
   integer(int32), intent(in   ), optional :: line

#if CONTRACTS_LEVEL >= 2
   character(:), allocatable :: msg
   if (.not.(cond)) then
       msg = "Intra-scope check failed"
       call stop_execution(msg, fileName, line)
   end if
#endif

   return
end subroutine check
