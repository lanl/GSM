! =============================================================================
!>
!> \file
!> \brief  Contains the procedures for the FException object
!> \author CMJ
!>
! =============================================================================

  ! ===========================================================================
  ! Return file name the exception occurred on
  pure function location(this) result(str)
    implicit none
    class(FException), intent(in) :: this
    character(:), allocatable :: str
    str = this%errorLocation
  end function location

  ! ===========================================================================
  ! Return exception text
  pure function what(this) result(errText)
    implicit none
    class(FException), intent(in) :: this
    character(:), allocatable :: errText
    errText = this%errorText
  end function what
  
  ! ===========================================================================
  ! Create and throw an exception
  subroutine throw(this, errFile, errLine, errText)
    implicit none
    class(FException), intent(inout) :: this
    character(*),      intent(in   ) :: errFile
    integer(gsmInt8),  intent(in   ) :: errLine
    character(*),      intent(in   ), optional :: errText
    
    ! Set exception location
    write(this%errorLocation, "(A, i0)") &
         & newLine // " ^^^ Error at '" // errFile // "', line ", errLine
    
    ! Build the message text
    if (present(errText)) then
       call this%buildMessage(errText)
    else
       call this%buildMessage("")
    end if
    
    ! Add the location to the message if desired
    if (includeLocation) then
       this%errorText = this%errorText // this%errorLocation
    end if
    
    ! Now "throw" the error
    call throwException(this)
  end subroutine throw
  
  ! ===========================================================================
  !> \brief Throws a general exception message
  subroutine throwExceptionBase(fExcept)
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    implicit none
    class(FException), intent(inout) :: fExcept
    write(ERROR_UNIT, "(A)") fExcept%what()
    stop 1
  end subroutine throwExceptionBase
  
  ! ===========================================================================
  !> Sets the method used to "throw" and display exceptions
  subroutine setExceptionThrower (throwPtr)
    implicit none
    procedure(throwExceptionBase) :: throwPtr
    throwException => throwPtr
  end subroutine setExceptionThrower
  
