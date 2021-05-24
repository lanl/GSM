! =============================================================================
!
!> \file
!> \brief  Contains the BremsPhoton object getters (queries object state)
!> \author CMJ (XCP-3; LANL)
!
! =============================================================================

! =============================================================================
!
!> \fn    getTMin
!> \brief Gets the minimum brems. photon energy
!
! ARGUMENTS:
!> \param[in   ] this        The BremsPhoton object
!
! =============================================================================
  function getTMin(this) result(tMin)

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(BremsPhoton), intent(in   ) :: this
    real(real64) :: tMin

    tMin = this%tgmin
    return
! =============================================================================
  end function getTMin


! =============================================================================
!
!> \fn    getTMax
!> \brief Gets the maximum brems. photon energy
!
! ARGUMENTS:
!> \param[in   ] this        The BremsPhoton object
!
! =============================================================================
  function getTMax(this) result(tMax)

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(BremsPhoton), intent(in   ) :: this
    real(real64) :: tMax

    tMax = this%tgmax
    return
! =============================================================================
  end function getTMax


! =============================================================================
!
!> \fn    getTEqv
!> \brief Gets the brems. photon equivalent KE
!
! ARGUMENTS:
!> \param[in   ] this        The BremsPhoton object
!
! =============================================================================
  function getTEqv(this) result(tEqv)

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(BremsPhoton), intent(in   ) :: this
    real(real64) :: tEqv

    tEqv = this%tEquiv
    return
! =============================================================================
  end function getTEqv


! =============================================================================
!
!> \fn    getSXAbs
!> \brief Gets the brems. photon absolute scattering cross section
!
! ARGUMENTS:
!> \param[in   ] this        The BremsPhoton object
!
! =============================================================================
  function getSXAbs(this) result(sxAbs)

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(BremsPhoton), intent(in   ) :: this
    real(real64) :: sxAbs

    sxAbs = this%absSX
    return
! =============================================================================
  end function getSXAbs


