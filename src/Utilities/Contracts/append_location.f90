! =============================================================================
!
!> \file
!> \brief  Definition of appending line information to contract message
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================
subroutine append_location(msg, fileName, line)
   use iso_fortran_env, only: int32
   implicit none
   character(:),   intent(inout), allocatable :: msg
   character(*),   intent(in   ),              optional :: fileName
   integer(int32), intent(in   ),              optional :: line

   character(len = 10) :: lineStr
   character(:), allocatable :: fileStr
   character(:), allocatable :: lineFmt

   if (showLine) then
       ! Dynamically set filename, if not provided
      if (present(fileName)) then
         fileStr = fileName
      else
         fileStr = "(not specified)"
      end if

      ! Dynamically set line number as a string
      if (present(line)) then
         write(lineStr, '(I0)') line
         lineFmt = ": " // lineStr
      else
         lineFmt = ""
      end if

      ! Modify the message
      msg = msg // newLine // "   ^^^ at " // fileStr // lineFmt

   end if
   return
end subroutine append_location
