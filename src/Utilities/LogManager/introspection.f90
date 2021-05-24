! =============================================================================
!
!> \file
!> \brief  Contains the LogManager introspection (state getters) procedures
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================



! =============================================================================
!
!> \fn    getNumLoggers
!> \brief Returns the number of Loggers contained in the LogManager object
!
! ARGUMENTS:
!> \param[in   ] this    The LogManager object
!
! =============================================================================
function getNumLoggers(this) result(numLoggers)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(LogManager), intent(in   ):: this
    integer(int32):: numLoggers

    numLoggers = this%d_numLoggers
    return
end function getNumLoggers


! =============================================================================
!
!> \fn    getMaxLoggers
!> \brief Returns the size of the logger array in the LogManager object
!
! ARGUMENTS:
!> \param[in   ] this    The LogManager object
!
! =============================================================================
function getMaxLoggers(this) result(maxLoggers)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(LogManager), intent(in   ):: this
    integer(int32):: maxLoggers

    maxLoggers = this%d_maxLoggers
    return
end function getMaxLoggers


! =============================================================================
!
!> \fn    getLogger
!> \brief Returns the i^th logger in the LogManager object
!
! ARGUMENTS:
!> \param[in   ] this    The LogManager object
!> \param[in   ] logIndx The index of the Logger object
!
! =============================================================================
function getLogger(this, logIndx) result(log)
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: Logger
    implicit none
    class(LogManager), intent(in   ):: this
    integer(int32),    intent(in   ):: logIndx
    type(Logger):: log

    if (this%numLoggers() <= 0) then
        ! Do nothing - no logger exists
    else if (logIndx <= 0) then
        ! Return the first Logger
        log = this%d_loggers(1)
    else if (logIndx >= this%numLoggers()) then
        ! Return the specified Logger
        log = this%d_loggers(logIndx)
    else
        ! Return the last Logger
        log = this%d_loggers(this%numLoggers())
    end if

    return
end function getLogger


! =============================================================================
!
!> \fn    describe
!> \brief Calls to each of the LogManager's describe functions
!
! ARGUMENTS:
!> \param[in   ] this    The LogManager object
!
! =============================================================================
subroutine describe(this)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(LogManager), intent(in   ):: this

    integer(int32):: indx
    do indx = 1, this%numLoggers()
        call this%d_loggers(indx)%describe()
    end do

    return
end subroutine describe
