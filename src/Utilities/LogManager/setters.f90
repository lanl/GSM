! =============================================================================
!
!> \file
!> \brief  Contains the LogManager setter and validation procedures
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================



! =============================================================================
!
!> \fn    resize
!> \brief Resizes the LogManager's logger array
!
! ARGUMENTS:
!> \param[inout] this    The LogManager object
!> \param[in   ] newSize The new size of the logger array
!
! =============================================================================
subroutine resize(this, newSize)
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: Logger
    implicit none
    class(LogManager),  intent(inout):: this
    integer(int32),     intent(in   ):: newSize

    type(Logger), dimension(this%d_numLoggers):: oldLoggers

    ! Only resize when the new size is valid
    if (newSize > 0) then
        ! Copy previous array if Loggers existed
        if (this%numLoggers() > 0) oldLoggers(1:this%numLoggers()) = &
            & this%d_loggers(1:this%numLoggers())

        ! Resize the array
        if (allocated(this%d_loggers)) deallocate(this%d_loggers)
        allocate(this%d_loggers(newSize))
        this%d_maxLoggers = newSize

        ! Copy previous array of Loggers back to the array
        this%d_numLoggers = min(newSize, this%numLoggers())
        if (this%numLoggers() > 0) &
            & this%d_loggers(1:this%numLoggers()) = &
            & oldLoggers(1:this%numLoggers())

    else

        ! Desired size will yield array without contents - deallocate entirely
        if (allocated(this%d_loggers)) deallocate(this%d_loggers)
        this%d_numLoggers = 0_int32
        this%d_maxLoggers = 0_int32

    end if

    return
end subroutine resize


! =============================================================================
!
!> \fn    addLogger
!> \brief Adds a Logger object to the array
!
! ARGUMENTS:
!> \param[inout] this      The LogManager object
!> \param[in   ] newLogger The Logger to be added
!
! =============================================================================
subroutine addLogger(this, newLogger)
    use, intrinsic:: iso_fortran_env, only: int32
    use Loggers, only: Logger
    implicit none
    class(LogManager),  intent(inout):: this
    class(Logger),      intent(in   ):: newLogger

    ! Increase the array size if appropriate
    if ((this%numLoggers() + 1_int32) > this%maxLoggers()) &
        & call this%resize (this%numLoggers() + 1_int32)

    ! Add the logger to the array
    this%d_numLoggers = this%d_numLoggers + 1
    this%d_loggers(this%d_numLoggers) = newLogger
    
    return
end subroutine addLogger

