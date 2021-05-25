
  subroutine checkMomentum(gsmObj, results, conserved)

! ======================================================================
!
! Performs the conservation of momentum check within GSM
!
!
! Written by CMJ, XCP-3 (03/2019); parsed from original "gsm03.f90" routine
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: hlf, one, two
    
    implicit none
    class(GSM),        intent(inout) :: gsmObj
    class(GSMResults), intent(in   ) :: results
    logical,           intent(  out) :: conserved

    real(real64)   :: maxMomentum, totalSystemMomentum, totParticleMomentum

    ! For per nucleon check:
    real(real64)   :: approxNumNucleons

    ! For system total check:
    integer(int32) :: pIndx
    real(real64)   :: pxTot, pyTot, pzTot, progenyMom

! ======================================================================

    ! Flag default value:
    conserved = .TRUE.

    ! Obtain total system momentum:
    totalSystemMomentum = sqrt((results%initialProj%kinEnergy + &
         & two * results%initialProj%restMass) * results%initialProj%kinEnergy)

    if ( .not.gsmObj%options%conserveTotalMomentum ) then

       ! Obtain an expected 'average' momentum per nucleon in the system
       approxNumNucleons = results%initialTarg%numBaryons + &
            & max(results%initialProj%restMass, dble(results%initialProj%numBaryons))
       maxMomentum = results%targRes%numBaryons * &
            & totalSystemMomentum/(approxNumNucleons) + hlf

       ! Total Momenta of Recoiled Particle (cannot exceed projectile's momentum)
       ! NOTE: The momentum check SHOULD include the projectile residual's momentum
       !       and ALL PROGENY as well, HOWEVER this is NOT advisable due to the
       !       renormalization routine!
       totParticleMomentum = results%targRes%totalLinearMomentum()
       ! totParticleMomentum = sqrt( (results%targRes%linearMom(1) + results%projRes%linearMom(1) )**2 + &
       !      & (results%targRes%linearMom(2) + results%projRes%linearMom(2) )**2 + &
       !      (results%targRes%linearMom(3) + results%projRes%linearMom(3) )**2 )

    else

       ! Check for total momentum conservation:
       maxMomentum = totalSystemMomentum + hlf
       pxTot = 0
       pyTot = 0
       pzTot = 0
       do pIndx = 1, results%numProgeny
          progenyMom = sqrt( results%progenyBnk(pIndx)%kinEnergy * &
               & (results%progenyBnk(pIndx)%kinEnergy + &
               & two * results%progenyBnk(pIndx)%restMass ) )
          pxTot = pxTot + progenyMom * results%progenyBnk(pIndx)%sinTheta * &
               & cos(results%progenyBnk(pIndx)%phi )
          pyTot = pyTot + progenyMom * results%progenyBnk(pIndx)%sinTheta * &
               & sin(results%progenyBnk(pIndx)%phi )
          pzTot = pzTot + progenyMom * results%progenyBnk(pIndx)%cosTheta
       end do
       pxTot = pxTot + results%targRes%linearMom(1) + results%projRes%linearMom(1)
       pyTot = pyTot + results%targRes%linearMom(2) + results%projRes%linearMom(2)
       pzTot = pzTot + results%targRes%linearMom(3) + results%projRes%linearMom(3)
       totParticleMomentum = sqrt( pxTot**2 + pytot**2 + pzTot**2 )

    end if

    ! Check if the momentum is greater than that allowed in the check:
    ! NOTE: If conserving on average, ensure NOT in anti-lab system, otherwise
    !       check for conservation.
    ! NOTE: For checking on total, anti-lab system does NOT need accounted for.
    if ( totParticleMomentum > maxMomentum .and. &
         & (results%initialProj%system /= antilabSystem .or. &
         & gsmObj%options%conserveTotalMomentum) ) then
       conserved = .FALSE.
    endif

    return
! ======================================================================
  end subroutine checkMomentum
