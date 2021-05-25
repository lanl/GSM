! =============================================================================
!
!> \file
!> \brief  Contains the Logger introspection (state getters) procedures
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================



! =============================================================================
!
!> \fn    getLevel
!> \brief Returns the Logger's maximum message level printed
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!
! =============================================================================
function getLevel(this) result(level)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger), intent(in   ):: this
    integer(int32) :: level

    level = this%b_level
    return
end function getLevel


! =============================================================================
!
!> \fn    getStream
!> \brief Returns the Logger's stream that messages get printed to
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!
! =============================================================================
function getStream(this) result(stream)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger), intent(in   ):: this
    integer(int32) :: stream

    stream = this%b_stream
    return
end function getStream


! =============================================================================
!
!> \fn    getNumFilters
!> \brief Returns the number of the Logger's message type filters
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!
! =============================================================================
function getNumFilters(this) result(numFilters)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger), intent(in   ):: this
    integer(int32):: numFilters

    numFilters = this%b_numFilters
    return
end function getNumFilters


! =============================================================================
!
!> \fn    getFilters
!> \brief Returns the Logger's message type filters
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!
! =============================================================================
function getFilters(this) result(filters)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger), intent(in   ):: this
    integer(int32), dimension(:), allocatable :: filters

    filters = this%b_filters
    return
end function getFilters


! =============================================================================
!
!> \fn    filter
!> \brief Returns a boolean to flag if the message type ought to be filtered
!
! ARGUMENTS:
!> \param[in   ] this    The Logger object
!> \param[in   ] level   The message level
!
! =============================================================================
function filter(this, level) result(filterMessage)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger),  intent(in   ):: this
    integer(int32), intent(in   ):: level
    logical:: filterMessage

    integer(int32):: indx

    filterMessage = .FALSE.
    filterLoop: do indx = 1, this%b_numFilters
       if (level == this%b_filters(indx)) then
           filterMessage = .TRUE.
           exit filterLoop
       end if
    end do filterLoop

    return
end function filter


! =============================================================================
!
!> \fn    printLevel
!> \brief Returns a boolean to flag for if a message level ought to be printed
!
! ARGUMENTS:
!> \param[in   ] this           The Logger object
!> \param[in   ] messageLevel   The message level
!
! =============================================================================
function printLevel(this, messageLevel) result(shouldPrint)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger),  intent(in   ):: this
    integer(int32), intent(in   ):: messageLevel
    logical:: shouldPrint

    ! Initialize default and check if the message's level is below that of the
    ! Logger
    shouldPrint = .FALSE.
    if ((messageLevel <= this%level()) .and. &
        & (.not.this%filter(messageLevel))) shouldPrint = .TRUE.

    return
end function printLevel


! =============================================================================
!
!> \fn    describe
!> \brief Describes the Logger object to its own stream
!
! ARGUMENTS:
!> \param[in   ] this           The Logger object
!
! =============================================================================
subroutine describe(this)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger), intent(in   ):: this
    integer(int32):: indx
    integer(int32), dimension(this%b_numFilters):: filters

    ! Write  basic info:
    write(this%stream(), 1000) this%stream() &
        &, this%level() &
        &, this%numFilters()

    ! Write out filters if they exist:
    if (this%numFilters() > 0) then
        filters = this%filters()
        write(this%stream(), 1100, advance="no")
        do indx = 1, this%numFilters()
            write(this%stream(), 1110, advance="no") filters(indx)
        end do
        write(this%stream(), 1200)
    end if

    return
! =============================================================================
1000 format(/, &
          & /5x, "LOGGER DESCRIPTION:", &
          & /5x, "==================================", &
          & /5x, "Stream       = ", i5, &
          & /5x, "Output level = ", i5, &
          & /5x, "Num. Filters = ", i5)
1100 format(/5x, "     =")
1110 format(" ", i2)
1200 format(/)
! =============================================================================
end subroutine describe

