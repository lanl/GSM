
! ==============================================================================
!
!> \file
!> \brief   Contains the implementation of \c getValidIOStream
!> \author  CMJ, XCP-3 (LANL)
!
!  PROCEDURE DESCRIPTION:
!
!> \fn     getValidIOStream
!> \brief  Returns a valid I/O stream from the provided one
!
!> Verifies that the provided I/O stream is available to write to. This
!> verifies the stream is valid (\c >0) and that the stream is open. In the case
!> of the stream being \c output_unit or \c error_unit, it will simply use that
!> stream (not validation necessary).
!> An invalid stream will return a default stream of \c output_unit.
!
!  ARGUMENTS:
!> \param[in   ] ioStream     the stream desired being written to
!> \return       newStream    the valid stream (coming from \c ioStream)
!
! ==============================================================================

  function getValidIOStream(ioStream) result(newStream)

      use, intrinsic:: iso_fortran_env, only: int32, &
           & input_unit, output_unit, error_unit

      implicit none
      integer(int32), intent(in   ) :: ioStream
      integer(int32) :: newStream

      logical :: isValid = .true.    !< Flags if the given stream is valid
      logical :: isOpen  = .false.   !< Flags if the file unit is open

! ==============================================================================

      ! Set default values
      newStream = output_unit

      ! Check if stream is positive
      if (ioStream < 0 .or. ioStream == input_unit) isValid = .false.

      ! Validate the unit is open IF ioStream is not the output/error_unit
      if (isValid .and. ioStream /= output_unit .and. ioStream /= error_unit) then
         inquire(ioStream, opened = isOpen)
         if (.not.isOpen) isValid = .false.
      end if

      ! Set new stream
      if (isValid) then
         newStream = ioStream
      end if

      return
! ==============================================================================
  end function getValidIOStream
