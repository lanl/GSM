
  subroutine setMDCMReaction(gsmObj, proj, targ)

! ====================================================================
!
! Sets the projectile/target for the mDCM calculation (interfaces to
! mDCM common blocks).
!
! NOTE: ONLY nuclei information is transferred - the mDCM is utilized
!       to fill in the remaining the common block's variables
!
!
! Written by CMJ, XCP-3 (05/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMClass, only: setupMDCMReaction

    implicit none
    class(GSM),           intent(inout) :: gsmObj
    class(gsmProjectile), intent(in   ) :: proj
    class(gsmTarget),     intent(in   ) :: targ

    real(real64) :: tempA

! ====================================================================

    ! mDCM Common Blocks:
    ! (proj/targ nuclei)
    real(real64) :: anucl1, anucl2,znucl1,znucl2,t0,eps1,eps2,vpi,a1,a2, &
         & c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    common/hcasc/anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2,vpi,a1,a2, &
         & c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2
    ! (rest mass/decay quantum number)
    real(real64) :: stin, amin
    common /stin/ stin,amin
    ! (anti-lab system)
    integer(int32) :: kobr
    real(real64) :: blab, glab
    common/bglab/blab,glab,kobr   ! Anti-lab system 

! ====================================================================

    ! Allow only one thread to enter this region at a time
    !$OMP critical

    ! Transfer projectile information:
    anucl1 = dble(proj%numBaryons)
    znucl1 = dble(proj%numProtons)
    stin = dble(proj%decayNumber)
    if ( proj%system == labSystem ) then
       amin = proj%restMass
    else if ( proj%system == antiLabSystem ) then
       amin = targ%restMass
    else
       ! Default to lab system:
       amin = proj%restMass
    end if

    ! Transfer target information:
    anucl2 = targ%numBaryons
    znucl2 = targ%numProtons

    ! Set incident energy based on the system (ensure usage of correct incident energy):
    if ( proj%system == labSystem ) then
       tempA = dble(proj%numBaryons)
    else if ( proj%system == antiLabSystem ) then
       tempA = targ%numBaryons
    else
       ! Default to lab system:
       tempA = dble(proj%numBaryons)
    end if
    if ( tempA > 1 ) then
       t0 = proj%kinEnergy / tempA
    else
       t0 = proj%kinEnergy
    end if

    ! Set lab system (lab=0, antiLab=1):
    kobr = proj%system

    ! Setup all internals for mDCM:
    call setupMDCMReaction()

    !$OMP end critical

    return
! ====================================================================
  end subroutine setMDCMReaction
