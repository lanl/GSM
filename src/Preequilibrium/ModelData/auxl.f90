
  subroutine auxl (preeqData, a, z, angmom, bf0, ln, erotev, delu)

! ======================================================================
!
!   This subroutine takes care of some bookkeeping related to
!   angular momentum, rotational energy and fission-barrier heights.
!   It was extracted from the old PRECOF routine.
!
!   CALLS: BF 
!
!    Written by A. J. Sierk, LANL T-16, October, 2003.
!    Modified to remove filling of eb and egs arrays; AJS, Feb, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!    Edited by CMJ, XCP-3, Aug. 2018 (Preeq. class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use preequilibriumParams, only: zro, one

    implicit none
    class(PreequilibriumData), intent(inout) :: preeqData
    real(real64),              intent(in   ) :: a           ! Residual mass number
    real(real64),              intent(in   ) :: z           ! Residual's charge (atomic number)
    real(real64),              intent(inout) :: angmom(3)   ! Angular Momentum
    real(real64),              intent(  out) :: bf0         ! Fission Barrier [MeV]
    integer(int32),            intent(  out) :: ln          ! Angular Momentum (quantum number)
    real(real64),              intent(  out) :: erotev      ! Rotational Energy [MeV]
    real(real64),              intent(  out) :: delu        ! Excess Energy [MeV]

    integer(int32) :: ia1, iz1, ln1
    real(real64)   :: egs0, ratio, temp, tempZ, um2

! ======================================================================

    real(real64), parameter :: erotevLimit = one   ! Limits rotational energy  [MeV/A]
    real(real64)            :: erotevMax   = zro   ! = aResidual * erotevLimit [MeV]

! ======================================================================

    ! Initialize variables
    delu = zro
    um2 = angmom(1)**2 + angmom(2)**2 + angmom(3)**2
    um2 = max(um2, zro)        ! Total angular momentum
    ln = nint(sqrt(um2))       ! Angular momentum quantum number
    iz1 = nint(preeqData%compound%numProtons - z) + 1  ! Number of protons
    ia1 = nint(preeqData%compound%numBaryons - a) + 1  ! Number of neutrons
    ln1 = ln + 1

  ! Obtain ground state energies and fission barriers
    if (iz1 > 0 .and. ia1 > 0) then
       if (iz1 <= dataEgsEbDim1 .and. ia1 <= dataEgsEbDim2 .and. ln1 <= dataEgsEbDim3) then
          bf0 = preeqData%compEnergy%eb(iz1, ia1, ln1)
          egs0 = preeqData%compEnergy%egs(iz1, ia1, ln1)
       else
          tempZ = z
          if ( tempZ < 0.01 ) then
             ! Check for when rare neutrons clusters reach this point: Temporary fix! (CMJ)
             write(preeqData%io%message,2000) a, z
             call preeqData%io%print(1, 3, preeqData%io%message)
             tempZ = 1
          end if
          bf0 = preeqData%fissBarr%bf (a, tempZ, ln, egs0)
       endif
    else
       bf0 = preeqData%fissBarr%bf (a, z, ln, egs0)
    endif
    erotev = egs0   ! Ground state energy [MeV]


!   For very light nuclei, limit the rotational energy to 1 MeV
!   per nucleon; renormalize L and excitation energy.
    erotevMax = erotevLimit * a
    if (erotev > erotevMax) then
       temp = abs(erotevMax/erotev)   ! Renormalization
       if ( temp < zro ) then
          temp = 0.1d0
          write(preeqData%io%message,1100) "71"
          call preeqData%io%print(4, 3, preeqData%io%message)
       end if
       ratio = sqrt(temp)         ! Ratio for the renormalization
       delu = erotev - erotevMax   ! Excess energy (lost due to renormalization)
       erotev = erotevMax          ! Set rotational energy [MeV]

       ! Re-obtain angular momentum
       angmom(1) = ratio*angmom(1)
       angmom(2) = ratio*angmom(2)
       angmom(3) = ratio*angmom(3)
       temp = um2
       if ( temp < zro ) then
          temp = 0.1d0
          write(preeqData%io%message,1000) "82"
          call preeqData%io%print(4, 3, preeqData%io%message)
       end if
       ln = nint(ratio*sqrt(temp))   ! Re-obtain angular momentum quantum number
    endif

    return
! ======================================================================
1000 format("Divide by zero error prevented in 'auxl.f90', line ", A)
1100 format("Square root error prevented in 'auxl.f90', line ", A)
2000 format("Large neutron cluster (A = ", f4.1, ", Z = ", &
          & f4.1, ") exceeds array bounds in ", &
          & "'auxl' routine.")
! ======================================================================
  end subroutine auxl
