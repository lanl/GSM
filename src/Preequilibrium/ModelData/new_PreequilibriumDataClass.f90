
  function new_PreequilibriumData(projBaryons, projProtons, &
       & targBaryons, targProtons, kinEnergyMeV, &
       & clientFissBarr, clientIO) &
       & result(preeqData)

! ======================================================================
!
! Constructor for the PreequilibriumData class.
!
! USE:
!    preeqData = Preequilibrium(projBaryons, projProtons, targBaryons, &
!         & targProtons, kinEnergyMeV, fissBarrObj, [clientIO])
!
!
! REQUIRED ARGUMENTS: Users are required to input the projectile information
! (A, Z, Kinetic energy [kinEnergyMeV]) and the target information (A, Z only).
!
! OPTIONAL ARGUMENTS:
! (1) clientIO may be speicifed so that users may control where I/O
!     handling in this class goes. Currently, I/O is restricted mostly
!     to the "auxl" routine for errors and warnings.
!
!
!
! Written by CMJ, XCP-3, 8/2018 (Preequil. class creation)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use preequilibriumParams, only: zro, hlf, thousand
    use fissionBarrierClass,  only: FissionBarrier

    implicit none
    integer(int32),       intent(in   ) :: projBaryons
    integer(int32),       intent(in   ) :: projProtons
    real(real64),         intent(in   ) :: targBaryons
    real(real64),         intent(in   ) :: targProtons
    real(real64),         intent(in   ) :: kinEnergyMeV
    type(FissionBarrier), intent(in   ), target :: clientFissBarr
    procedure(IOHANDLER), intent(in   ), pointer, optional :: clientIO

    type(PreequilibriumData)   :: preeqData


    integer(int32) :: ia1, iz1, ln1, ln
    real(real64)   :: aa, bf0, bf, daz, egs0, z

! ======================================================================

    ! Set where I/O goes
    if ( present(clientIO) ) then
       if ( associated(clientIO) ) then
          preeqData%io%print => clientIO
       else
          write(preeqData%io%message, 1000)
          call preeqData%io%print(2, 3, preeqData%io%message)
       end if
    end if


    ! Point to necessary data:
    preeqData%fissBarr  => clientFissBarr                     ! Point to appropriate fission barrier object
    preeqData%molEnergy => preeqData%fissBarr%queryMolnix()   ! Obtains the "molnix" object


    ! Establish ground state and fission barrier energies for all possible nuclei
    preeqData%compound%numBaryons = targBaryons + dble(projBaryons)
    preeqData%compound%numProtons = targProtons + dble(projProtons)
    preeqData%compound%kinEnergy  = kinEnergyMeV / thousand
    do iz1 = 1,dataEgsEbDim1
       z = preeqData%compound%numProtons - dble(iz1 - 1)
       if (z > hlf) then
          ! Repeats for total A, Z greater than 0
          do ia1 = 1,dataEgsEbDim2
             aa = preeqData%compound%numBaryons - dble(ia1 - 1)
             if (aa > hlf) then
                do ln1 = 1,dataEgsEbDim3
                   ln = ln1 - 1  ! Total Angular Momentum
                   call preeqData%fissBarr%barfit (aa, z, ln, bf0, egs0)
                   daz = preeqData%molEnergy%shellEnergy(aa, z)
                   bf = bf0 - daz
                   bf = max(bf, zro)

                   ! Store ground state and fission barrier energies
                   preeqData%compEnergy%eb(iz1, ia1, ln1) = bf
                   preeqData%compEnergy%egs(iz1, ia1, ln1) = egs0
                end do
             endif
          end do
       endif
    end do


    ! Establish F_j array now
    preeqData%gammaJObj = GammaJ(kinEnergyMeV, targBaryons, projBaryons, projProtons)


    preeqData%constructed = .TRUE.

    return
! ======================================================================
1000 format("The provided I/O handler to the preeq. data object is ", &
          & "not associated and will not be used.")
! ======================================================================
  end function new_PreequilibriumData
