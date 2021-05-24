! =============================================================================
!
!> \file
!> \brief  Contains the BremsPhoton object introspection functions
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

! =============================================================================
!
!> \fn    incrementTEqv
!> \brief Increments the equivalent KE of the photon by some value
!
! ARGUMENTS:
!> \param[inout] self        The BremsPhoton object
!> \param[in   ] val         The value to increment by
!
! =============================================================================
  subroutine incrementTEqv(self, val)

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(BremsPhoton), intent(inout) :: self
    real(real64),       intent(in   ) :: val

    self%tEquiv = self%tEquiv + val
    return
! =============================================================================
  end subroutine incrementTEqv


! =============================================================================
!
!> \fn    incrementSXAbs
!> \brief Increments the absolute scattering cross section of the photon by
!>        some value
!
! ARGUMENTS:
!> \param[inout] self        The BremsPhoton object
!> \param[in   ] val         The value to increment by
!
! =============================================================================
  subroutine incrementSXAbs(self, val)

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(BremsPhoton), intent(inout) :: self
    real(real64),       intent(in   ) :: val

    self%absSX = self%absSX + val
    return
! =============================================================================
  end subroutine incrementSXAbs


! =============================================================================
!
!> \fn    sampleEnergy
!> \brief Returns the sampled energy of the brems. photon via Schiff's spectrum
!
! ARGUMENTS:
!> \param[inout] self        The BremsPhoton object
!> \param[in   ] rndm        A random number
!
! =============================================================================
  function sampleEnergy(self, rndm) result(energy)

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(BremsPhoton), intent(in   ) :: self
    real(real64),       intent(in   ) :: rndm
    real(real64) :: energy

    energy = self%tMin() * (self%tMax() / self%tMin())**(rndm)
    
    return
! =============================================================================
  end function sampleEnergy

