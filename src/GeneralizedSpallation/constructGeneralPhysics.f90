
  subroutine constructGeneralPhysics( gsmObj )

! ==============================================================================
!
! This function constructs the sub-models used by GSM in the generalPhysicsModels
! data type
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

    use standardDCMClass,      only: newStandardDCM
    use fermiBreakUpClass,     only: newFermiBreakup
   use evaporationClass,      only: newEvaporation
 
    implicit none
    class(GSM), intent(inout) :: gsmObj

! ==============================================================================

    gsmObj%genModels%numConstructionErrors = 0


    ! Standard DCM:
    gsmObj%genModels%sDCM = newStandardDCM( &
         & gsmObj%rang, &
         & gsmObj%genData%molnixEnergies, &
         & gsmObj%options%sDCMOpts, &
         & gsmObj%io%print)
    if ( .not.gsmObj%genModels%sDCM%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) &
            & "The Standard DCM model failed to properly construct."
       call gsmObj%io%print(1, 1, gsmObj%io%message )
       gsmObj%genModels%numConstructionErrors = gsmObj%genModels%numConstructionErrors + 1
    end if
    gsmObj%options%sDCMOpts = gsmObj%genModels%sDCM%queryOptions()


    ! Fermi Break-Up:
    gsmObj%genModels%fbu = newFermiBreakUp( &
         & gsmObj%rang, &
         & gsmObj%options%fbuOpts, &
         & gsmObj%io%print )
    if ( .not.gsmObj%genModels%fbu%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) &
            & "The Fermi Break-Up model failed to properly construct."
       call gsmObj%io%print(1, 1, gsmObj%io%message )
       gsmObj%genModels%numConstructionErrors = gsmObj%genModels%numConstructionErrors + 1
    end if
    gsmObj%options%fbuOpts = gsmObj%genModels%fbu%queryOptions()


    ! Evaporation:
    if ( .not. gsmObj%usePhotonEmission ) then
       gsmObj%genModels%evap = newEvaporation( &
            & gsmObj%genData%evap, &
            & gsmObj%rang, &
            & gsmObj%genData%molnixEnergies, &
            & gsmObj%genModels%fbu, &
            & gsmObj%options%evapOpts, &
            & clientIO=gsmObj%io%print )
    else
       gsmObj%genModels%evap = newEvaporation( &
            & gsmObj%genData%evap, &
            & gsmObj%rang, &
            & gsmObj%genData%molnixEnergies, &
            & gsmObj%genModels%fbu, &
            & gsmObj%options%evapOpts, &
            & evapGammaCascade, &
            & gsmObj%io%print )
    end if
     if ( .not.gsmObj%genModels%evap%properlyConstructed() ) then
       write(gsmObj%io%message, 1000) &
            & "The Evaporation model failed to properly construct."
       call gsmObj%io%print(1, 1, gsmObj%io%message )
       gsmObj%genModels%numConstructionErrors = gsmObj%genModels%numConstructionErrors + 1
    end if
    gsmObj%options%evapOpts = gsmObj%genModels%evap%queryOptions()


    ! Check how many errors occurred:
    if ( gsmObj%genModels%numConstructionErrors /= 0 ) then
       write(gsmObj%io%message, 1100) gsmObj%genModels%numConstructionErrors
       call gsmObj%io%print(1, 1, gsmObj%io%message)
    end if


    return
! ==============================================================================
1000 format(A)
1100 format(i2, " general sub-models failed to construct within GSM.")
! ==============================================================================
  end subroutine constructGeneralPhysics
