
  function inquireINCModel( gsmObj, proj, targ ) result(incFlag)

! ====================================================================
!
! Returns an integer flag indicating which INC model to use
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    class(GSM),           intent(inout) :: gsmObj
    class(GSMProjectile), intent(inout) :: proj
    class(GSMTarget),     intent(inout) :: targ
    integer(int32) :: incFlag

! ====================================================================

    ! Assume the sDCM is used:
    incFlag = sDCMFlagged

    ! Check for anti-lab system: use mDCM if found
    if ( proj%numBaryons > targ%numBaryons ) call gsmObj%formNuclei(proj, targ)
    if ( proj%system == antiLabSystem ) then
       incFlag = mDCMFlagged
       return
    end if

    ! For incident light- and heavy-ions, use the mDCM (always):
    if ( proj%numBaryons > 1 ) then
       incFlag = mDCMFlagged
       return
    end if


    ! Determine use based on small projectiles:
    select case( trim(proj%particleName) )
       case( projNames(1), projNames(2) )
          ! Nucleon
          if ( proj%kinEnergy <= &
               & gsmObj%sampleEnergy( gsmObj%options%nucleonTransitionE, &
               & gsmObj%options%transitionWidth ) ) then
             incFlag = sDCMFlagged
          else
             incFlag = mDCMFlagged
          end if
       case( projnames(3), projNames(4), projNames(5) )
          ! Pion
          if ( proj%kinEnergy <= &
               & gsmObj%sampleEnergy( gsmObj%options%pionTransitionE, &
               & gsmObj%options%transitionWidth ) ) then
             incFlag = sDCMFlagged
          else
             incFlag = mDCMFlagged
          end if
       case( projNames(6) )
          ! Photon (mono-energetic)
          if ( proj%kinEnergy <= &
               & gsmObj%sampleEnergy( gsmObj%options%monoPhotonTransitionE, &
               & gsmObj%options%transitionWidth ) ) then
             incFlag = sDCMFlagged
          else
             incFlag = mDCMFlagged
          end if
       case( projNames(7) )
          ! Photon (brems.)
          if ( proj%kinEnergy <= &
               & gsmObj%sampleEnergy( gsmObj%options%bremPhotonTransitionE, &
               & gsmObj%options%transitionWidth ) ) then
             incFlag = sDCMFlagged
          else
             incFlag = mDCMFlagged
          end if
       case default
          ! Assume nucleon (sample accordingly)
          if ( proj%kinEnergy <= &
               & gsmObj%sampleEnergy( gsmObj%options%nucleonTransitionE, &
               & gsmObj%options%transitionWidth ) ) then
             incFlag = sDCMFlagged
          else
             incFlag = mDCMFlagged
          end if
    end select

    return
! ====================================================================
  end function inquireINCModel
