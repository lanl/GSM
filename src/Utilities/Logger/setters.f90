! =============================================================================
!
!> \file
!> \brief  Contains the Logger setter and validation procedures
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================



! =============================================================================
!
!> \fn    setLevel
!> \brief Sets the Logger's maximum message type to be printed
!
! ARGUMENTS:
!> \param[inout] this    The Logger object
!> \param[in   ] level   The new level
!
! =============================================================================
subroutine setLevel(this, level)
    use, intrinsic:: iso_C_binding, only: c_int
    implicit none
    class(Logger),  intent(inout):: this
    integer(c_int), intent(in   ):: level

    this%b_level = this%validateLevel(level)
    return
end subroutine setLevel


! =============================================================================
!
!> \fn    validateLevel
!> \brief Returns a valid message type level
!
! ARGUMENTS:
!> \param[inout] this    The Logger object
!> \param[in   ] level   The new level
!
! =============================================================================
function validateLevel(this, level) result(validLevel)
    use, intrinsic:: iso_C_binding, only: c_int
    implicit none
    class(Logger),  intent(in   ):: this
    integer(c_int), intent(in   ):: level
    integer(c_int):: validLevel

    select case (level)
       case (QUIETLog); validLevel = QUIETLog
       case (CRITICALLog); validLevel = CRITICALLog
       case (ERRORLog); validLevel = ERRORLog
       case (WARNINGLog); validLevel = WARNINGLog
       case (INFOLog); validLevel = INFOLog
       case (DEBUGLog); validLevel = DEBUGLog
       case (UNSETLog); validLevel = UNSETLog
       case default;
           validLevel = defaultLevel
           call this%warning("An invalid message level was provided to the Logger.")
           call this%warning("   The default will be used.")
    end select
    return
end function validateLevel


! =============================================================================
!
!> \fn    setStream
!> \brief Sets the stream of the Logger object
!
! ARGUMENTS:
!> \param[inout] this    The Logger object
!> \param[in   ] stream  The stream to be printed to
!
! =============================================================================
subroutine setStream(this, stream)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger),  intent(inout):: this
    integer(int32), intent(in   ):: stream

    this%b_stream = validateStream(stream)
    return
end subroutine setStream


! =============================================================================
!
!> \fn    validateStream
!> \brief Returns a valid stream unit
!
! ARGUMENTS:
!> \param[in   ] stream  The stream to be printed to
!
! =============================================================================
function validateStream(stream) result(validStream)
    use, intrinsic:: iso_fortran_env, only: int32, &
        & input_unit, output_unit, error_unit
    implicit none
    integer(int32), intent(in   ):: stream
    integer(int32):: validStream

    logical:: exists

    ! Set the default and ensure positive
    validStream = stream
    if (validStream <= 0_int32 .or. validStream == input_unit) &
        & validStream = defaultStream

    ! If output- or error-unit, don't check for unit existence
    if (validStream /= output_unit .and. validStream /= error_unit) then
        !> \todo Revisit special handling for printing to files, i.e. does
        !> the Logger manage the file (open/close) and assign a name, or do we
        !> throw our hands in the air and let the client do it if they want the
        !> unit printed to?
        inquire(validStream, exist=exists)
        if (.not.exists) validStream = defaultStream
    end if

    return
end function validateStream


! =============================================================================
!
!> \fn    resizeFilters
!> \brief Resizes the filters used by the Logger
!
! ARGUMENTS:
!> \param[inout] this    The Logger object
!> \param[in   ] length  The new length of the filters
!
! =============================================================================
subroutine resizeFilters(this, length)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger),  intent(inout):: this
    integer(int32), intent(in   ):: length

    integer(int32), dimension(this%b_numFilters):: tempFilters

    ! Set the number of filters based on the given argument
    if (length < 0) then
        this%b_numFilters = 0
    else if (length > numMessageLevels) then
        this%b_numFilters = numMessageLevels
    else
        this%b_numFilters = length
    end if

    ! Temporarily store the current filters
    tempFilters = this%filters()

    ! Resize now
    if (allocated(this%b_filters)) deallocate(this%b_filters)
    allocate(this%b_filters(this%b_numFilters))

    ! Restore the previous values now
    this%b_filters(1:min(this%b_numFilters, size(tempFilters))) = &
        & tempFilters(1:min(this%b_numFilters, size(tempFilters)))

    return
end subroutine resizeFilters


! =============================================================================
!
!> \fn    filterExists
!> \brief Returns a boolean flag if a filter already exists
!
! ARGUMENTS:
!> \param[inout] this    The Logger object
!> \param[in   ] filter  The filter inquired about
!
! =============================================================================
function filterExists(this, filter) result(exists)
    use, intrinsic::iso_fortran_env, only: int32
    implicit none
    class(Logger),  intent(in   ):: this
    integer(int32), intent(in   ):: filter
    logical:: exists

    integer(int32):: indx

    ! Set default and search through the filters to see if one exists
    filterLoop: do indx = 1_int32, this%b_numFilters
       if (filter == this%b_filters(indx)) then
           exists = .TRUE.
           exit filterLoop
       end if
    end do filterLoop

    return
end function filterExists


! =============================================================================
!
!> \fn    addFilters
!> \brief Adds message filters to the Logger
!
! ARGUMENTS:
!> \param[inout] this    The Logger object
!> \param[in   ] filters The filters to be applied
!
! =============================================================================
subroutine addFilters(this, filters)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger),  intent(inout):: this
    integer(int32), intent(in   ), dimension(:):: filters

    integer(int32):: indx

    ! Interface to addFilter
    filterLoop: do indx = 1_int32, size(filters)
       call this%addFilter(filters(indx))
    end do filterLoop

    return
end subroutine addFilters


! =============================================================================
!
!> \fn    addFilter
!> \brief Adds a single message filter to the Logger
!
! ARGUMENTS:
!> \param[inout] this    The Logger object
!> \param[in   ] filter  The filter to be applied
!
! =============================================================================
subroutine addFilter(this, filter)
    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(Logger),  intent(inout):: this
    integer(int32), intent(in   ):: filter

    logical:: validTypeFilter
    integer(int32), dimension(this%b_numFilters + 1_int32):: newFilters

    ! Ensure the filter is a valid message type:
    select case(filter)
       case (CRITICALLog); validTypeFilter = .TRUE.
       case (ERRORLog); validTypeFilter = .TRUE.
       case (WARNINGLog); validTypeFilter = .TRUE.
       case (INFOLog); validTypeFilter = .TRUE.
       case (DEBUGLog); validTypeFilter = .TRUE.
       case default; validTypeFilter = .FALSE.
    end select

    ! Add the filter is its valid and doesn't exist)
    if (validTypeFilter) then
        if (.not.this%filterExists(filter)) then
            ! Resize the filter array
            call this%resizeFilters(size(newFilters))

            ! Add the new value to the end
            this%b_filters(this%b_numFilters) = filter
        end if
    end if

    return
end subroutine addFilter

