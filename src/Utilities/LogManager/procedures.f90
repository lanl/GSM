! =============================================================================
!
!> \file
!> \brief  Contains the LogManager's general procedures
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

!> @{
!> \brief Interfaces the LogManager's message to the appropriate printing
!> procedure
subroutine criticalExtended(this)
    class(LogManager), intent(in):: this
    call this%criticalBase(this%message)
    return
end subroutine criticalExtended
subroutine errorExtended(this)
    class(LogManager), intent(in):: this
    call this%errorBase(this%message)
    return
end subroutine errorExtended
subroutine warningExtended(this)
    class(LogManager), intent(in):: this
    call this%warningBase(this%message)
    return
end subroutine warningExtended
subroutine infoExtended(this)
    class(LogManager), intent(in):: this
    call this%infoBase(this%message)
    return
end subroutine infoExtended
subroutine debugExtended(this)
    class(LogManager), intent(in):: this
    call this%debugBase(this%message)
    return
end subroutine debugExtended
subroutine printExtended(this)
    class(LogManager), intent(in):: this
    call this%printBase(this%message)
    return
end subroutine printExtended
!> @}


! =============================================================================
!
!> \fn    prepareMessage
!> \brief Prepends the message type and removes trailing/leading spaces
!
! ARGUMENTS:
!> \param[in   ] level   The level of the message
!> \param[in   ] message The message being printed
!
! =============================================================================
function prepareMessage(level, message) result(prepared)
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: prepLogMessage => prepareMessage
    implicit none
    integer(int32), intent(in   ):: level
    character(*),   intent(in   ):: message
    character(:), allocatable:: prepared

    prepared = message
    call prepLogMessage(level, prepared)

    return
end function prepareMessage


! =============================================================================
!
!> \fn    criticalBase
!> \brief Prints a critical message to the LogManager's Loggers
!
! ARGUMENTS:
!> \param[in   ] this    The LogManager object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine criticalBase(this, message)
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: CRITICALLog
    implicit none
    class(LogManager), intent(in   ):: this
    character(*),      intent(in   ):: message
    character(:), allocatable:: string
    integer(int32):: indx

    string = prepareMessage(CRITICALLog, message)
    logLoop: do indx = 1_int32, this%numLoggers()
        call this%d_loggers(indx)%critical(string)
    end do logLoop

    ! Kill the thread here
    ! \todo Add macro call to obtian a line and file/"exception handling"
    error stop "A critical error was encountered in the LogManager!"

    return
end subroutine criticalBase


! =============================================================================
!
!> \fn    errorBase
!> \brief Prints an error message to the LogManager's Loggers
!
! ARGUMENTS:
!> \param[in   ] this    The LogManager object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine errorBase(this, message)
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: ERRORLog
    implicit none
    class(LogManager), intent(in   ):: this
    character(*),      intent(in   ):: message
    character(:), allocatable:: string
    integer(int32):: indx

    string = prepareMessage(ERRORLog, message)
    logLoop: do indx = 1_int32, this%numLoggers()
        call this%d_loggers(indx)%error(string)
    end do logLoop

    return
end subroutine errorBase


! =============================================================================
!
!> \fn    warningBase
!> \brief Prints a warning message to the LogManager's Loggers
!
! ARGUMENTS:
!> \param[in   ] this    The LogManager object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine warningBase(this, message)
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: WARNINGLog
    implicit none
    class(LogManager), intent(in   ):: this
    character(*),      intent(in   ):: message
    character(:), allocatable:: string
    integer(int32):: indx

    string = prepareMessage(WARNINGLog, message)
    logLoop: do indx = 1_int32, this%numLoggers()
        call this%d_loggers(indx)%warning(string)
    end do logLoop

    return
end subroutine warningBase


! =============================================================================
!
!> \fn    infoBase
!> \brief Prints an info message to the LogManager's Loggers
!
! ARGUMENTS:
!> \param[in   ] this    The LogManager object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine infoBase(this, message)
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: INFOLog
    implicit none
    class(LogManager), intent(in   ):: this
    character(*),      intent(in   ):: message
    character(:), allocatable:: string
    integer(int32):: indx

    string = prepareMessage(INFOLog, message)
    logLoop: do indx = 1_int32, this%numLoggers()
        call this%d_loggers(indx)%info(string)
    end do logLoop

    return
end subroutine infoBase


! =============================================================================
!
!> \fn    debugBase
!> \brief Prints a debug message to the LogManager's Loggers
!
! ARGUMENTS:
!> \param[in   ] this    The LogManager object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine debugBase(this, message)
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: DEBUGLog
    implicit none
    class(LogManager), intent(in   ):: this
    character(*),      intent(in   ):: message
    character(:), allocatable:: string
    integer(int32):: indx

    string = prepareMessage(DEBUGLog, message)
    logLoop: do indx = 1_int32, this%numLoggers()
        call this%d_loggers(indx)%debug(string)
    end do logLoop

    return
end subroutine debugBase


! =============================================================================
!
!> \fn    printBase
!> \brief Prints a general message to the LogManager's Loggers
!
! ARGUMENTS:
!> \param[in   ] this    The LogManager object
!> \param[in   ] message The message being printed
!
! =============================================================================
subroutine printBase(this, message)
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: UNSETLog
    implicit none
    class(LogManager), intent(in   ):: this
    character(*),      intent(in   ):: message
    character(:), allocatable:: string
    integer(int32):: indx

    string = prepareMessage(UNSETLog, message)
    logLoop: do indx = 1_int32, this%numLoggers()
        call this%d_loggers(indx)%print(string)
    end do logLoop

    return
end subroutine printBase

