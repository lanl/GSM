
  subroutine formNuclei( gsmObj, proj, targ )

! ====================================================================
!
! Establishes/sets the variables in the projectile/target nuclei for a
! given event
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use BremsPhotonMod, only: newBremsPhoton
    use gsm_params, only: emprot, emneut, massPiPM, massPi0

    implicit none
    class(GSM),           intent(inout) :: gsmObj
    class(GSMProjectile),  intent(inout) :: proj
    class(GSMTarget),      intent(inout) :: targ

    real(real64)      :: temp
    character(len= 4) :: pname
    character(len=10) :: name

! ====================================================================

    ! Set projectile defaults
    if ( proj%numBaryons <  0 ) proj%numBaryons = defaultProjBaryons
    if ( proj%numProtons < -1 ) proj%numProtons = defaultProjProtons

    ! Ensure all values are physically valid (use default if not):
    if ( proj%numBaryons <= proj%numProtons .and. proj%numBaryons > 1 ) then
       ! A is greater than Z for a nucleus; use defaults:
       write(gsmObj%io%message, 1100) "projectile", proj%numBaryons, proj%numProtons
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       proj%numBaryons = defaultProjBaryons
       proj%numProtons = defaultProjProtons
    end if
    if ( proj%kinEnergy  <  0 ) proj%kinEnergy  = defaultKinEnergy
    if ( targ%numBaryons <  0 ) targ%numBaryons = defaultTargBaryons
    if ( targ%numProtons < -1 ) targ%numProtons = defaultTargProtons
    if ( targ%numBaryons <= targ%numProtons .and. targ%numBaryons > 1 ) then
       ! A is greater than Z for a nucleus; use defaults:
       write(gsmObj%io%message, 1100) "target", nint(targ%numBaryons), nint(targ%numProtons)
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       targ%numBaryons = defaultTargBaryons
       targ%numProtons = defaultTargProtons
    end if

    ! Check if anti-lab system is to be used:
    if ( proj%numBaryons > targ%numBaryons ) then
       call gsmObj%swapNuclei(proj, targ)
    else
       proj%system = labSystem
    end if

    ! Populate projectile information:
    pname = trim(proj%particleName)
    select case( pname )
       case( projNames(1) )
          ! Proton
          proj%numBaryons = 1
          proj%numProtons = 1
          proj%restMass   = emprot
          proj%particleFlag = nucleusProjFlag
       case( projNames(2) )
          ! Neutron
          proj%numBaryons = 1
          proj%numProtons = 0
          proj%restMass   = emneut
          proj%particleFlag = nucleusProjFlag
       case( projNames(3) )
          ! pi+
          proj%numBaryons = 0
          proj%numProtons = 1
          proj%restMass   = massPiPM
          proj%particleFlag = pionProjFlag
       case( projNames(4) )
          ! Pi-
          proj%numBaryons = 0
          proj%numProtons = -1
          proj%restMass   = massPiPM
          proj%particleFlag = pionProjFlag
       case( projNames(5) )
          ! Pi0
          proj%numBaryons = 0
          proj%numProtons = 0
          proj%restMass   = massPi0
          proj%particleFlag = pionProjFlag
       case( projNames(6) )
          ! Mono-energetic photon:
          proj%numBaryons = 0
          proj%numProtons = 0
          proj%restMass   = 0.0
          proj%particleFlag = photonProjFlag
       case( projNames(7) )
          ! Brems. photon:
          proj%numBaryons = 0
          proj%numProtons = 0
          proj%restMass   = 0
          proj%particleFlag = bremsProjFlag

          ! Bremsstrahlung spectrum f(E)dE=const*dE/E, tgmin < E < tgmax
          if (.not.associated(proj%brems)) allocate(proj%brems)
          proj%brems = newBremsPhoton(proj%kinEnergy, proj%kinEnergyMax)
       case default
          ! Light- or Heavy-ion:
          ! Use the provided A/Z flags
          if ( proj%restMass <= 0 ) proj%restMass = &
               & determineRestMass( proj%numBaryons, proj%numProtons )
          proj%particleFlag = nucleusProjFlag
    end select

    ! Populate target information:
    if ( targ%restMass <= 0 ) targ%restMass = &
         & determineRestMass( nint(targ%numBaryons), nint(targ%numProtons) )

     ! Modify Af/Cz multipliers for the evaporation simulation
     call gsmObj%updateLevelDensities(proj, targ)


    ! Check energy range (warn client if outside recommendation)
    if ( 1000 * proj%kinEnergy < minRecommendedEnergy ) then
       ! Below sDCM threshold:
       write(gsmObj%io%message, 1200) 1000*proj%kinEnergy, "MeV", "below"
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       write(gsmObj%io%message, 1210) minRecommendedEnergy, "MeV"
       call gsmObj%io%print(2, 3, gsmObj%io%message)
    else if ( proj%kinEnergy / 1000 > maxRecommendedEnergy ) then
       ! Above mDCM threshold:
       write(gsmObj%io%message, 1200) proj%kinEnergy/1000, "TeV", "above"
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       write(gsmObj%io%message, 1210) maxRecommendedEnergy, "TeV"
       call gsmObj%io%print(2, 3, gsmObj%io%message)
    end if

    return
! ====================================================================
1100 format("The ", A, " is invalid (A=", i3, "Z=", i3, "). Using defaults.")
1200 format("The projectile has incident energy (", f8.3, " ", A, ") ", A)
1210 format("   the recommendation for this model (", f8.3, " ", A, ").")
! ====================================================================
  end subroutine formNuclei
