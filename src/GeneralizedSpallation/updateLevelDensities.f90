
  subroutine updateLevelDensities(gsmObj, proj, targ)

! ====================================================================
!
! This procedure is used to update the level density multipliers of
! the projectile and target objects based on the A/Z of the projectile
! and target as well as the incident energy [GeV] of the projectile.
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro
    use evaporationClass, only: setLevelDensity

    implicit none
    class(GSM),           intent(inout) :: gsmObj
    class(GSMProjectile), intent(inout) :: proj
    class(GSMTarget),     intent(inout) :: targ

    integer(int32) :: afCzFlag = 0_int32
    real(real64)   :: afMult, czMult

! ====================================================================

    if ( gsmObj%options%useSDCMAfCzMultipliers ) then
       afCzFlag = 0_int32
    else
       afCzFlag = 1_int32
    end if

    ! Set Af and Cz multipliers, if applicable, for evap/fission physics:
    ! PROJECTILE:
    call setLevelDensity( dble(proj%numBaryons), &
         & dble(proj%numProtons), &
         & proj%kinEnergy, afMult, czMult, afCzFlag )
    if ( proj%afMultiplier < zro ) &
       proj%afMultiplier = afMult
    if ( proj%czMultiplier < zro ) &
         & proj%czMultiplier = czMult

    ! TARGET:
    call setLevelDensity( targ%numBaryons, &
         & targ%numProtons, &
         & proj%kinEnergy, afMult, czMult, afCzFlag )
    if ( targ%afMultiplier < zro ) &
         & targ%afMultiplier = afMult
    if ( targ%czMultiplier < zro ) &
         & targ%czMultiplier = czMult

    return
! ====================================================================
  end subroutine updateLevelDensities
