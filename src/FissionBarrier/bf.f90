
  function bf (fbObj, a, z, ln, egs0)

! ======================================================================
!
!   Function to give fission barrier in Mev, given A and Z.
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Modified for CEM97a, December, 1997.
!   Modified by AJS, March, 1999.
!   "Last" change: 12-AUG-2003 by NVMokhov
!   Edited by A. J. Sierk, LANL T-16, October, 2003.
!   Modified by AJS, February, 2005.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fissionBarrierParams, only: zro

    implicit none
    class(FissionBarrier), intent(inout) :: fbObj
    real(real64),          intent(in   ) :: a
    real(real64),          intent(in   ) :: z
    integer(int32),        intent(in   ) :: ln
    real(real64),          intent(  out) :: egs0
    real(real64)                         :: bf

    real(real64) :: bf0, daz

! ======================================================================

!***********************************************************************
!  Sierk [Phys. Rev. C 33, 2039 (1986)] global fit for macroscopic Bf. *
!  Yukawa-plus-exponential nuclear energy; diffuse-surface Coulomb     *
!  energy; diffuse-matter moments of inertia.                          *
!***********************************************************************
!    Bf0 = bf without shell and odd-even corrections to the g.s. mass
    call fbObj%barfit (a, z, ln, bf0, egs0)

    ! Apply shell correction
    daz = fbObj%fbMolnix%shellEnergy(a, z)
    bf = bf0 - daz

    ! Ensure physically allowable
    bf = max(bf, zro)


    return
! ======================================================================
  end function bf
