! =============================================================================
!
!> \file
!> \brief  Contains the Logger's general procedures
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================


! =============================================================================
!
!> \fn    cleanString
!> \brief Cleans a string (i.e. allocatable character) by removing whitespace
!
!> Cleans an allocatable character by removing the leading and trailing
!> whitespace of an allocatable character.
!
! ARGUMENTS:
!> \param[inout] string  The string to be cleaned
!
! =============================================================================
subroutine cleanString(string)
    implicit none
    character(:),   intent(inout), allocatable:: string

    string = trim(adjustl(string))
    return
end subroutine cleanString


! =============================================================================
!
!> \fn    prepend
!> \brief Modifies a message with the level type preprended to it
!
! ARGUMENTS:
!> \param[in   ] level   The level of the message
!> \param[inout] message The message to be prepended to
!
! =============================================================================
subroutine prepend(level, message)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    integer(int32), intent(in   ):: level
    character(:),   intent(inout), allocatable:: message

    ! Prepend text as appropriate:
    select case (level)
       case (QUIETLog); return
       case (CRITICALLog); message = preCritical // message
       case (ERRORLog); message = preError // message
       case (WARNINGLog); message = preWarning // message
       case (INFOLog); message = preInfo // message
       case (DEBUGLog); message = preDebug // message
       case (UNSETLog); message = preUnset // message
       case default; message = preCritical // message
    end select

    return
end subroutine prepend


! =============================================================================
!
!> \fn    prepareMessage
!> \brief Cleans and prepends message type to the message
!
! ARGUMENTS:
!> \param[in   ] level   The level of the message
!> \param[inout] message The message being printed
!
! =============================================================================\
subroutine prepareMessage(level, message)
    implicit none
    integer(int32), intent(in   ):: level
    character(:),   intent(inout), allocatable:: message
    
    call cleanString(message)
    call prepend(level, message)

    return
end subroutine prepareMessage


! =============================================================================
!
!> \fn    baseLog
!> \brief Prints a message to the Logger's stream
!
!> Base procedure that prints a given message to the Logger object's stream. The
!> calling procedures simply perform any message pre- and post-processing.
!
!> \todo Modify this procedure to limit the frequency of buffer printing, such
!> as creating an array of allocatable messages and printing several at once
!> (culminating the messages here).
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!> \param[in   ] level   The level of the message
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine baseLog(this, level, message)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger),  intent(in   ):: this
    integer(int32), intent(in   ):: level
    character(*),   intent(in   ):: message

    ! Check if the Logger prints this type of message and print if so
    if (this%printLevel(level)) &
        & write(this%stream(), stringFormat) message

    return
end subroutine baseLog


! =============================================================================
!
!> \fn    critical
!> \brief Prints a critical message to the Logger's stream
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine critical(this, message)
    implicit none
    class(Logger),  intent(in   ):: this
    character(*),   intent(in   ):: message

    ! Print the message
    call this%baseLog(CRITICALLog, message)

    ! Kill the thread here:
    !> \todo Add exception module for handling this error; something like:
    !> call criticalError("A critical error was encountered: " // string, __FILE__)
    error stop "A critical error was encountered!"

    return
end subroutine critical


! =============================================================================
!
!> \fn    error
!> \brief Prints an error message to the Logger's stream
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine error(this, message)
    implicit none
    class(Logger),  intent(in   ):: this
    character(*),   intent(in   ):: message

    ! Print the message
    call this%baseLog(ERRORLog, message)

    return
end subroutine error


! =============================================================================
!
!> \fn    warning
!> \brief Prints a warning message to the Logger's stream
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine warning(this, message)
    implicit none
    class(Logger),  intent(in   ):: this
    character(*),   intent(in   ):: message
    
    ! Print the message
    call this%baseLog(WARNINGLog, message)

    return
end subroutine warning


! =============================================================================
!
!> \fn    info
!> \brief Prints an info message to the Logger's stream
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine info(this, message)
    implicit none
    class(Logger),  intent(in   ):: this
    character(*),   intent(in   ):: message

    ! Print the message
    call this%baseLog(INFOLog, message)

    return
end subroutine info


! =============================================================================
!
!> \fn    debug
!> \brief Prints a debug message to the Logger's stream
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine debug(this, message)
    implicit none
    class(Logger),  intent(in   ):: this
    character(*),   intent(in   ):: message

    ! Print the message
    call this%baseLog(DEBUGLog, message)

    return
end subroutine debug


! =============================================================================
!
!> \fn    print
!> \brief Prints an unspecified message type to the Logger's stream
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine print(this, message)
    implicit none
    class(Logger),  intent(in   ):: this
    character(*),   intent(in   ):: message

    ! Print the message
    call this%baseLog(UNSETLog, message)

    return
end subroutine print

