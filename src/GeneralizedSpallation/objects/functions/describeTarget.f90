
! ==============================================================================
!
!> \file
!> \brief   Contains the implementation of \c describeTarget
!> \author  CMJ, XCP-3 (LANL)
!
!  PROCEDURE DESCRIPTION:
!
!> \fn     describeTarget
!> \brief  Writes the attributes of the target to an output stream
!
!> Writes the attributes of the target object to an output stream. Clients
!> may provide an output stream to the procedure, however it will be validated
!> prior to being utilized.
!> The default output stream is to terminal (\c output_unit)
!
!  ARGUMENTS:
!> \param[in   ] targ          the target to be described
!> \param[in   ] ioStream      the unit number (int32) to write to (opt.)
!
! ==============================================================================

  subroutine describeTarget(targ, ioStream)

      use, intrinsic:: iso_fortran_env, only: int32, real64, output_unit

      implicit none
      class(gsmTarget), intent(in   ) :: targ
      integer(int32),   intent(in   ), optional :: ioStream

      integer(int32) :: stream = output_unit   !< The actual stream variable 

! ==============================================================================

      ! Set stream
      if (present(ioStream)) then
         stream = getValidIOStream(ioStream)
      end if

      ! Describe attributes:
      write(stream, 1000) &
         & trim(adjustl(targ%particleName)), &
         & targ%numBaryons, &
         & targ%numProtons, &
         & targ%restMass, &
         & targ%afMultiplier, &
         & targ%czMultiplier

      return
! ==============================================================================
1000 format(/, &
          & /5x, "TARGET '", A, "':", &
          & /5x, "--------------------------------------", &
          & /5x, "A    = ", f5.1, &
          & /5x, "Z    = ", f5.1, &
          & /5x, "E_0  = ", f8.3, " [GeV/c**2]", &
          & /5x, "A_f  = ", f6.3, &
          & /5x, "C_z  = ", f6.3)
! ==============================================================================
  end subroutine describeTarget

