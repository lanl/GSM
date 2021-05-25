! =============================================================================
!
!> \file
!> \brief  Tests a Logger object's level implementation
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================



! =============================================================================
!
!> \fn    tstLoggerLevel
!> \brief Tests the a Logger's level implementation
!
! =============================================================================
subroutine tstLoggerLevel
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: Logger, newLogger, &
        & WARNINGLog, UNSETLog
    implicit none

    ! Construct logger with an invalid level
    type(Logger):: log
    if (log%level() /= UNSETLog) &
        & error stop "Initial level failed to be default!"
   
   
    log = newLogger(24)
    if(log%level() == 24) &
        & error stop "Logger failed at construction!"

    ! Set logger to a valid level:
    call log%setLevel(WARNINGLog)
    if(log%level() /= WARNINGLog) &
        & error stop "Logger failed to set a valid level!"

    ! Set logger to an invalid level:
    call log%setLevel(-WARNINGLog)
    if(log%level() /= UNSETLog) &
        & error stop "Logger failed to recognize error with an invalid level!"

    return
end subroutine tstLoggerLevel
