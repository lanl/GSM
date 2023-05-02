! =============================================================================
!
!> \file
!> \brief  Contains the Validate macro
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

#ifndef Validate
#define Validate(COND, MSG) \
    if (.not. (COND)) validate(COND, MSG, __FILE__, __LINE__)
!   if (.not.(COND)) then \
!       character(len=:), allocatable :: error_line \
!       write(error_line, "(I5)") __LINE__ \
!       character(len=:), allocatable :: error_str \
!       error_str = "Validation failed: " // trim(adjustl(MSG) \
!       error_str = error_str // & \
!           & NEWLINE // "^^^ at " // __FILE__ // ":" // error_line \
!           \
!       error stop error_str \
!   endif
#endif
