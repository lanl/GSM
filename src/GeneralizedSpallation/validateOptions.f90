
  subroutine validateOptions( gsmObj )

! ======================================================================
!
! Validates the GSM-specific options that the GSM class is to use
!
! ======================================================================

    implicit none
    class(GSM), intent(inout) :: gsmObj

! ======================================================================

    ! Transition width:
    if ( gsmObj%options%smoothTransition ) then
       if ( gsmObj%options%transitionWidth < 0 ) then
          write(gsmObj%io%message, 1000) "transition energy", gsmObj%options%transitionWidth
          call gsmObj%io%print(2, 3, gsmObj%io%message)
          gsmObj%options%transitionWidth = transitionWidthDefault
          write(gsmObj%io%message, 1050) gsmObj%options%transitionWidth
          call gsmObj%io%print(2, 3, gsmObj%io%message)
       end if
    end if

    ! Transition energies:
    ! NOTE: Clients may use these options to restrict sDCM/mDCM use for
    !       specific incident particles.
    !       i.e. nucleonTransition = 0 means ONLY mDCM usage (conversely for sDCM)
    ! (nucleons):
    if ( gsmObj%options%nucleonTransitionE < 0 ) then
       write(gsmObj%io%message, 1000) "nucleon transition energy", gsmObj%options%nucleonTransitionE
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       gsmObj%options%nucleonTransitionE = nucleonTransEDefault
       write(gsmObj%io%message, 1050) gsmObj%options%nucleonTransitionE
       call gsmObj%io%print(2, 3, gsmObj%io%message)
    end if
    ! (pions):
    if ( gsmObj%options%pionTransitionE < 0 ) then
       write(gsmObj%io%message, 1000) "pion transition energy", gsmObj%options%pionTransitionE
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       gsmObj%options%pionTransitionE = pionTransEDefault
       write(gsmObj%io%message, 1050) gsmObj%options%pionTransitionE
       call gsmObj%io%print(2, 3, gsmObj%io%message)
    end if
    ! (photons; mono-energetic):
    if ( gsmObj%options%monoPhotonTransitionE < 0 ) then
       write(gsmObj%io%message, 1000) "mono-energetic photon transition energy", &
            & gsmObj%options%monoPhotonTransitionE
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       gsmObj%options%monoPhotonTransitionE = monoPhotonTransEDefault
       write(gsmObj%io%message, 1050) gsmObj%options%monoPhotonTransitionE
       call gsmObj%io%print(2, 3, gsmObj%io%message)
    end if
    ! (photons; brems.):
    if ( gsmObj%options%bremPhotonTransitionE < 0 ) then
       write(gsmObj%io%message, 1000) "brems. photon transition energy", &
            & gsmObj%options%bremPhotonTransitionE
       call gsmObj%io%print(2, 3, gsmObj%io%message)
       gsmObj%options%bremPhotonTransitionE = bremPhotonTransEDefault
       write(gsmObj%io%message, 1050) gsmObj%options%bremPhotonTransitionE
       call gsmObj%io%print(2, 3, gsmObj%io%message)
    end if

    return
! ======================================================================
1000 format("An invalid ", A, " was detected (", es15.4,").")
1050 format("   Using default value (", es15.4, ")")
! ======================================================================
  end subroutine validateOptions
