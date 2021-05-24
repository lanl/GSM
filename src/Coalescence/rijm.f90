
  subroutine  rijm (coalObj, results, j1, j2, rmin)

! ======================================================================
!
!    Distance of closest approach of 2 particles?????
!    NOT USED!!!   SGM, July, 2005.
! 
!   "Last" change: 13-AUG-2003 by NVMokhov
!   Edited by A. J. Sierk, LANL T-16, October, 2003.
!
! ======================================================================

    use iso_fortran_env, only: int32, real64
    use coalescenceParams, only : zro, avg_mass

    implicit none
    class(Coalescence),       intent(inout) :: coalObj
    type(coalescenceResults), intent(inout) :: results
    integer(int32),           intent(in   ) :: j1
    integer(int32),           intent(in   ) :: j2
    real(real64),             intent(  out) :: rmin

    integer(int32) :: k
    real(real64)   :: s, r2, v2
    real(real64), dimension(3) :: r12 = zro, v12 = zro

! ======================================================================

    ! (positions)
!    r12(1) = results%partBnk(j1)%xCoord - results%partBnk(j2)%xCoord
!    r12(2) = results%partBnk(j1)%yCoord - results%partBnk(j2)%yCoord
!    r12(3) = results%partBnk(j1)%zCoord - results%partBnk(j2)%zCoord

    ! (momenta)
    v12(1) = results%partBnk(j1)%linearMomX/(results%partBnk(j1)%kinEnergy + avg_mass) - &
         & results%partBnk(j2)%linearMomX/(results%partBnk(j2)%kinEnergy + avg_mass)
    v12(2) = results%partBnk(j1)%linearMomY/(results%partBnk(j1)%kinEnergy + avg_mass) - &
         & results%partBnk(j2)%linearMomY/(results%partBnk(j2)%kinEnergy + avg_mass)
    v12(3) = results%partBnk(j1)%linearMomZ/(results%partBnk(j1)%kinEnergy + avg_mass) - &
         & results%partBnk(j2)%linearMomZ/(results%partBnk(j2)%kinEnergy + avg_mass)

    s = r12(1)*v12(1) + r12(2)*v12(2) + r12(3)*v12(3)
    r2 = r12(1)**2 + r12(2)**2 + r12(3)**2
    v2 = v12(1)**2 + v12(2)**2 + v12(3)**2
    rmin = 100.d0
    if (s <= zro) then
       rmin = sqrt(abs(r2 - s**2/v2))
       rmin = sqrt(r2)
    endif
    return

! ======================================================================
  end subroutine rijm
