
  subroutine  codir (coalObj, results, ind, nn)

! ======================================================================
!
!    Coalescence "direction".  Average momentum of the nn coalesced
!    particles; placed into the results%partBnk array with the index of the first
!    particle.
!
!   "Last" change: 12-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by CMJ, LANL XCP-3, July 2018 (Coalescence class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use coalescenceParams, only : zro, one, avg_mass

    implicit none
    class(Coalescence),       intent(inout) :: coalObj
    type(coalescenceResults), intent(inout) :: results
    integer(int32),           intent(in   ) :: ind(7)
    integer(int32),           intent(in   ) :: nn

    integer(int32) :: i, i1, k
    real(real64)   :: ps(3), psm_sq, pmm

! ======================================================================

    ps(1) = zro
    ps(2) = zro
    ps(3) = zro
    do  k = 1,nn
       i = ind(k)
       ps(1) = ps(1) + results%partBnk(i)%linearMomX
       ps(2) = ps(2) + results%partBnk(i)%linearMomY
       ps(3) = ps(3) + results%partBnk(i)%linearMomZ
    end do

! For some reason ps is, very rarely, huge (1e+154), causing a floating error when squared
    if (ps(1)>1.d100 .or. ps(2)>1.d100 .or. ps(3)>1.d100) then
       ps(:) = zro
       write(coalObj%io%message, 1000)
       call coalObj%io%print(4, 3, coalObj%io%message)
    end if

    psm_sq = ps(1)**2 + ps(2)**2 + ps(3)**2
    i1 = ind(1)
    pmm = dble(nn)*avg_mass
    results%partBnk(i1)%kinEnergy    = sqrt(psm_sq + pmm*pmm) - pmm
    results%partBnk(i1)%restMass     = pmm
    results%partBnk(i1)%linearMomX   = ps(1)
    results%partBnk(i1)%linearMomY   = ps(2)
    results%partBnk(i1)%linearMomZ   = ps(3)
    results%partBnk(i1)%coalesceFlag = one

    return
! ======================================================================
1000 format("Infinite momentum detected in ", &
          & "'codir.f90', line 944. Error will be resolved; results ", &
          & "may be suspect.")
! ======================================================================
  end subroutine codir
