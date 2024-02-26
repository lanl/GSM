! =============================================================================
!
!> \file
!> \fn
!> \brief  Specifies that a method or feature has not been implemented.
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================
subroutine not_implemented(feature, fileName, line)
    use iso_fortran_env, only: int32
   implicit none
   character(*),   intent(in   ) :: feature
   character(*),   intent(in   ), optional :: fileName
   integer(int32), intent(in   ), optional :: line

   character(:), allocatable :: msg
   msg = "The feature '" // feature // "' is not implemented"
   call stop_execution(msg, fileName, line)

   return
end subroutine not_implemented
