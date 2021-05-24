
! ==============================================================================
!
!> \file
!> \brief   Contains the implementation of \c describeProjectile
!> \author  CMJ, XCP-3 (LANL)
!
!  PROCEDURE DESCRIPTION:
!
!> \fn     describeProjectile
!> \brief  Writes the attributes of the projectile to an output stream
!
!> Writes the attributes of the projectile object to an output stream. Clients
!> may provide an output stream to the procedure, however it will be validated
!> prior to being utilized.
!> The default output stream is to terminal (\c output_unit)
!
!  ARGUMENTS:
!> \param[in   ] proj          the projectile to be described
!> \param[in   ] ioStream      the unit number (int32) to write to (opt.)
!
! ==============================================================================

  subroutine describeProjectile(proj, ioStream)

      use, intrinsic:: iso_fortran_env, only: int32, real64, output_unit

      implicit none
      class(gsmProjectile), intent(in   ) :: proj
      integer(int32),       intent(in   ), optional :: ioStream

      integer(int32) :: stream = output_unit   !< The actual stream variable 
      character(20)  :: type = ""              !< String for particle type
      character(20)  :: system = ""            !< String for system type

! ==============================================================================

      ! Set stream
      if (present(ioStream)) then
         stream = getValidIOStream(ioStream)
      end if

      ! Obtain particle type flag
      if (proj%particleFlag == pionProjFlag) then
         type = "pion"
      else if (proj%particleFlag == photonProjFlag) then
         type = "photon (mono.)"
      else if (proj%particleFlag == bremsProjFlag) then
         type = "photon (brems.)"
      else
         type = "nucleus"
      end if

      ! Write the system utilized
      if (proj%system == labSystem) then
         system = "lab. system"
      else if (proj%system == antilabSystem) then
         system = "anti. lab. system"
      else
         system = "unspecified"
      end if

      ! Describe attributes:
      write(stream, 1000) &
         & trim(adjustl(proj%particleName)), &
         & proj%numBaryons, &
         & proj%numProtons, &
         & proj%decayNumber, &
         & trim(type), &
         & proj%restMass, &
         & proj%kinEnergy, &
         & proj%dKinEnergy, &
         & proj%kinEnergyMax, &
         & proj%afMultiplier, &
         & proj%czMultiplier, &
         & trim(system)

      return
! ==============================================================================
1000 format(/, &
          & /5x, "PROJECTILE '", A, "':", &
          & /5x, "--------------------------------------", &
          & /5x, "A         = ", i3, &
          & /5x, "Z         = ", i3, &
          & /5x, "decay no. = ", i3, &
          & /5x, "type      = ", A, &
          & /5x, "E_0       = ", f8.3, "   [GeV/c**2]", &
          & /5x, "E_k       = ", es10.4, " [GeV]", &
          & /5x, "  dE_k    = ", f8.3, "   [MeV]", &
          & /5x, "  E_k,max = ", es10.4, " [GeV]", &
          & /5x, "A_f       = ", f6.3, &
          & /5x, "C_z       = ", f6.3, &
          & /5x, "System    = ", A)
! ==============================================================================
  end subroutine describeProjectile
