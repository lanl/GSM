! =============================================================================
!
!> \file
!> \brief  Contains the BremsPhoton object setters and reset functions
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

! =============================================================================
!
!> \fn    setTMin
!> \brief Sets the minimum brems. photon energy
!
! ARGUMENTS:
!> \param[inout] this        The BremsPhoton object
!> \param[in   ] minE        Minimum energy of the brems. photon
!
! =============================================================================
  subroutine setTMin(this, minE)

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(BremsPhoton), intent(inout) :: this
    real(real64),       intent(in   ) :: minE

    this%tgmin = minE
    return
! =============================================================================
  end subroutine setTMin


! =============================================================================
!
!> \fn    setTMax
!> \brief Sets the maximum brems. photon energy
!
! ARGUMENTS:
!> \param[inout] this        The BremsPhoton object
!> \param[in   ] maxE        Maximum energy of the brems. photon
!
! =============================================================================
  subroutine setTMax(this, maxE)

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(BremsPhoton), intent(inout) :: this
    real(real64),       intent(in   ) :: maxE

    this%tgmax = maxE
    return
! =============================================================================
  end subroutine setTMax


! =============================================================================
!
!> \fn    resetTEqv
!> \brief Resets the equivalent brems. photon KE
!
! ARGUMENTS:
!> \param[inout] this        The BremsPhoton object
!
! =============================================================================
  subroutine resetTEqv(this)

    implicit none
    class(BremsPhoton), intent(inout) :: this

    this%tEquiv = 0.0_real64
    return
! =============================================================================
  end subroutine resetTEqv


! =============================================================================
!
!> \fn    resetSXAbs
!> \brief Resets the photon scattering cross section
!
! ARGUMENTS:
!> \param[inout] this        The BremsPhoton object
!
! =============================================================================
  subroutine resetSXAbs(this)

    implicit none
    class(BremsPhoton), intent(inout) :: this

    this%absSX = 0.0_real64
    return
! =============================================================================
  end subroutine resetSXAbs

