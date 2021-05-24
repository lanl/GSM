
  subroutine eventLoop(gsmObj, gsmRxn)

! ======================================================================
!
! Main physics loop for spallation events - simulates INC, coalescence,
! preequilibrium, evaporation, and Fermi Break-up emission where
! appropriate
!
!
! Written by CMJ, XCP-3 (03/2019); parsed from original "gsm03.f90" routine
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, int64, real64, error_unit
    use gsm_params, only: zro, hlf, one, two
    use randomNumberGenerator, only: RN_next_particle
    
    ! Some objects used (mostly for setup):
    use standardDCMData,  only: photonEG
    use evaporationClass, only : setLevelDensity

    ! For initializing things...
    use hist_mod, only: hist_tally

! ======================================================================

    implicit none
    class(GSM),        intent(inout) :: gsmObj
    type(GSMReaction), intent(inout), target :: gsmRxn

    integer(int32) :: ibadren, in, iz, ln
    real(real64)   :: bf0, delu, erotev, pevap, &
         & u, ue, wam
    logical :: inelasticEvent = .TRUE.
    logical :: momentumConserved = .TRUE.
    logical :: mDCMUsed = .FALSE.

    type(GSMResults),    pointer :: results => NULL()
    type(GSMProjectile), pointer :: proj    => NULL()
    type(GSMTarget),     pointer :: targ    => NULL()

! ======================================================================

    results => gsmRxn%results
    proj    => results%initialProj
    targ    => results%initialTarg
    mDCMUsed = (gsmRxn%incFlag == mDCMFlagged)

30  continue   ! for cem reference
35  continue
    inelasticEvent = .TRUE.

    ! Reset simulation results object:
    call results%resetEvent()

    ! ----- for bremss. photons only ----- (sample incident energy)
    !   Sampling the energy of gamma from Schiff's spectra (1/E) and
    !   preparing g + N --> pi + N angular distributions:
    !   [The latter moved to TYPINT  02/25/05]
    if (results%initialProj%particleFlag == bremsProjFlag) then
       ! Sample energy and update energy for sDCM/mDCM simulation
       proj%kinEnergy = proj%brems%sampleEnergy(gsmObj%rang())
       gsmRxn%sDCMProj%kinenergy = proj%kinEnergy
    endif


    ! Sample INC if smooth transitions are desired:
    if ( gsmObj%options%smoothTransition .and. proj%numBaryons <= 1 ) then
       gsmRxn%incFlag = gsmObj%inquireINCModel( proj, targ )
       mDCMUsed = (gsmRxn%incFlag == mDCMFlagged)
    end if

    ! Performing IntraNuclear Cascade (INC) calculation
    select case( gsmRxn%incFlag )
       case( sDCMFlagged )
          call gsmObj%standardDCMInterface (gsmRxn%data, gsmRxn%sDCMProj, results)
       case( mDCMFlagged )
          !>>> CRITICAL REGION
          call gsmObj%modifiedDCMInterface( gsmRxn )
          !>>> END CRITICAL REGION
       case default
          call gsmObj%standardDCMInterface (gsmRxn%data, gsmRxn%sDCMProj, results)
    end select


    ! THESE 'IF' BLOCKS CHECK FOR PHYSICALITY (IF NOT, 'RESETS' THE EVENT
    ! [GO TO 30])

    ! Performing Coalescence, Pre-Equil, Equil, and Fermi Break-Up OR
    ! Evaporation/Fission when appropriate; unallowed reaction re-runs cascade.
    if (results%targRes%kinEnergy <= zro) then
       inelasticEvent = .FALSE.

       ! Negative Exit energy (elastic event); run cascade again (increment num. elastic events with CEM cascade only)
       if ( mDCMUsed ) then
          ! Print warning to user (shouldn't reach here)
          write(gsmObj%io%message, 1000) "An unexpected elastic interaction was flagged."
          call gsmObj%io%print(5, 3, gsmObj%io%message)
       else
          ! CEM cascade allows elastic events in this manner
          if (associated(gsmRxn%outData)) gsmRxn%outData%intel = gsmRxn%outData%intel + 1
       end if
    else
       if (results%numElasticEvents > 0 .and. .not.mDCMUsed) then
          ! Elastic Event; run cascade again (no momentum imparted)
          if (associated(gsmRxn%outData)) gsmRxn%outData%intel = gsmRxn%outData%intel + 1
          inelasticEvent = .FALSE.
       else
          ! Inelastic Event occurred, obtaining data
          ! Rerun cascade if an unphysical particle nucleus was produced via CASCAD:
          if ( ( results%targRes%numBaryons < 4.d0 .or. &
               & results%targRes%numProtons < one  .or. &
               & results%targRes%numBaryons < results%targRes%numProtons) .and. &
               & .not.(mDCMUsed) ) then
             write(gsmObj%io%message,1100) &
                  & "The target size (A<=4) is not allowed when using the Standard DCM INC model in GSM", gsmRxn%eventNum
             call gsmObj%io%print(2, 3, gsmObj%io%message)
             go to 30
          endif

         
          ! ABOVE PHYSICS IS VALID - CONTINUE ON
          ! Renormalize residuals to ensure physicality
          if ( .not.mDCMUsed ) then
             ! NOTE: renorm() is NOT used with LAQGSM's INC becuase the recalculation of
             !       energy is not performed for the projectile, and in addition would
             !       require modification of the equations used. Mass/charge conservation
             !       is observed some of the time with the use of the LAQGSM INC, however
             !       not often. 'renorm()' has been updated to allow mass/charge conservation
             !       checking with the LAQGSM INC, however should not be used as it will
             !       provide incorrect results (due to A/Z conservation failure in 'CASCAW()')
             !>>> ASSIGN TO AN OPENMP TASK?
             !>>> {
             call gsmObj%renorm (results, ibadren)

             ! Bad Re-Normalization; conservation laws broken
             if ( ibadren > 0) then
                ! call gsmObj%printPartProp (proj, targ, outData%ncas, intel)  ! for debug.
                write(gsmObj%io%message,1100) "A bad renormalization occurred", gsmRxn%eventNum
                call gsmObj%io%print(4, 3, gsmObj%io%message)
                go to 30
             endif
             !>>> }
          endif

          !>>> ASSIGN TO AN OPENMP TASK?
          !>>> {
          ! Form compound nuclei from the emitted INC progeny
          call gsmObj%coalescenceInterface (gsmRxn%models, results, gsmRxn%outData)

          ! Check for momentum conservation:
          call gsmObj%checkMomentum(results, momentumConserved)
          if ( .not. momentumConserved ) then
             write(gsmObj%io%message,1100) "Conservation of momentum failed", gsmRxn%eventNum
             call gsmObj%io%print(4, 3, gsmObj%io%message)
             go to 30
          endif
          !>>> }

          ! Converting to MeV
          results%targRes%kinEnergy = results%targRes%kinEnergy*1000.d0
          results%projRes%kinEnergy = results%projRes%kinEnergy*1000.d0

          ! Checking for allowed residual particle; LAQGSM Checks for this above
          if ( ( results%targRes%numBaryons >= 4.d0 .and. &
               & results%targRes%numProtons >= one  .and. &
               & results%targRes%numBaryons >= results%targRes%numProtons) .or. &
               & (mDCMUsed) ) then

             iz = nint(results%targRes%numProtons)       ! number of protons/charged particles
             in = nint(results%targRes%numBaryons - results%targRes%numProtons)  ! number of neutrons/neutral particles
             pevap = gsmObj%genData%molnixenergies%defineenergy (iz, in, 3)

             ! Momenta (ln), Angular Momenta (amnucl), Energy, and Fission Barrier (bf) Book-keeping
             call gsmRxn%data%preeq%auxl (results%targRes%numBaryons, results%targRes%numProtons, &
                  & results%targRes%angularMom, bf0, ln, erotev, delu)

             u = results%targRes%kinEnergy + delu
             ue = u - pevap - erotev ! rotational energy (ev)

             !>>> CRITICAL REGION (FOR TIME BEING)
             !>>> {
             if(results%tallySim) then
                ! Statistics at Various Stages of the Reaction
                call gsmObj%ststcs (results%targRes%numBaryons, &
                     & results%targRes%numProtons, ue, ln, bf0, 1)

                if (gsmRxn%output%printHist) then
                   ! Tallying Particle Data
                   call hist_tally(results%targRes%linearMom, results%targRes%angularMom, &
                        & results%targRes%kinEnergy, results%targRes%numBaryons, &
                        & results%targRes%numProtons, results%targExc%numNeutrons, &
                        & results%targExc%numProtons, results%targExc%numHoles )
                end if

                ! Valid event occurred, count event and notify user/client
                gsmRxn%outData%ncas = gsmRxn%outData%ncas + 1
                gsmRxn%eventNum = gsmRxn%eventNum + 1
                if ( mod(gsmRxn%eventNum, gsmObj%options%printIncrement) == 0 ) then
                   if ( gsmRxn%eventNum < 99999999 ) then
                      write (gsmObj%io%message, 2000) gsmRxn%eventNum
                   else
                      write (gsmObj%io%message, 2050) real(gsmRxn%eventNum)
                   endif
                   call gsmObj%io%print(2, 5, gsmObj%io%message)
                endif
             end if
             !>>> }

             ! ===================================================
             ! Tally Projectile (where appropriate)
             ! ===================================================
             if ( mDCMUsed ) then
                ! Count scores from residual projectile nucleus where appropriate
                if ( gsmObj%options%tallyResidualProjectile .and. &
                     & results%projRes%kinEnergy > zro .and. &
                     & results%projRes%numBaryons > 0.9 ) then

                   iz = nint(results%projRes%numProtons)
                   in = nint(results%projRes%numBaryons - results%projRes%numProtons)
                   pevap = gsmObj%genData%molnixenergies%defineenergy (iz, in, 3)

                   ! Momenta (ln), Angular Momenta (amnucl), Energy, and Fission Barrier (bf) Book-keeping
                   call gsmRxn%data%preeq%auxl (results%projRes%numBaryons, &
                        & results%projRes%numProtons, &
                        & results%projRes%angularMom, bf0, ln, erotev, delu)
                   u = results%projRes%kinEnergy + delu
                   ue = u - pevap - erotev ! rotational energy (ev)

                   !>>> CRITICAL REGION (OPENMP) FOR NOW
                   !>>> {
                   if(results%tallySim) then
                      ! Statistics at Various Stages of the Reaction
                      call gsmObj%ststcs (results%projRes%numBaryons, &
                           & results%projRes%numProtons, ue, ln, bf0, 1)

                      if (gsmRxn%output%printHist) then
                         ! Tallying Particle Data
                         call hist_tally(results%projRes%linearMom, results%projRes%angularMom, &
                              & results%projRes%kinEnergy, results%projRes%numBaryons, &
                              & results%projRes%numProtons, results%projExc%numNeutrons, &
                              & results%projExc%numProtons, results%projExc%numHoles )
                      end if
                   end if
                   !>>> }

                   ! Check for fermi break-results%projRes%kinEnergy or other emission
                   if ( gsmObj%genModels%fbu%recommendFermiBreakUp(results%projRes%numBaryons) ) then

                      ! Fermi Break-Up Driving Routine
                      call gsmObj%fermiBreakUpInterface (results%projRes, &
                           & gsmRxn)
                   else
                      ! Preequilibrium/Evaporation decay
                      results%info%wf = zro
                      wam = one

                      ! Calculating Pre-Equilibrium and Equilibrium Particle Emssion
                      call gsmObj%simulateDecay (gsmRxn, &
                           & results%projRes, results%projExc, &
                           & proj%afMultiplier, &
                           & proj%czMultiplier, &
                           & wam )

                      if (wam == -13.d0) then
                         wam = one
                         write(gsmObj%io%message,1100) "Projectile pre-equilibrium decay failure", gsmRxn%eventNum
                         call gsmObj%io%print(3, 3, gsmObj%io%message)

                         go to 30
                      endif
                   endif

                   ! Tally the projectile separately
                   results%projRes%kinEnergy = results%projRes%kinEnergy/1000.d0

                   ! if (gsmRxn%eventNum <= nnnp) call gsmObj%printPartProp (proj, &
                   !      & targ, gsmRxn%outData%ncas, intel)
                   if (results%info%fusion > zro .and. &
                        & gsmRxn%output%accountFission > 0) &
                        & gsmRxn%outData%nfis = gsmRxn%outData%nfis + 1
                   gsmRxn%outData%sfu = gsmRxn%outData%sfu + results%info%wf

                elseif ( gsmObj%options%tallyResidualProjectile ) then

                   ! Proj. residual not physical, check for errors
                   ! Check for physical residual projectile particle
                   if ( results%projRes%kinEnergy < zro .or. results%projRes%numProtons < zro ) then

                      ! Unphysical residual produced, rerun the cascade
                      write(gsmObj%io%message,1100) "An unphysical projectile residual was detected", gsmRxn%eventNum
                      call gsmObj%io%print(3, 3, gsmObj%io%message)

                      if(results%tallySim) then
                         ! Decrement counters
                         gsmRxn%outData%ncas = gsmRxn%outData%ncas - 1
                         gsmRxn%eventNum = gsmRxn%eventNum - 1
                      end if
                     if (associated(gsmRxn%outData)) gsmRxn%outData%intel = &
                          & gsmRxn%outData%intel + 1
                      inelasticEvent = .FALSE.

                      go to 35

                   endif
                endif
             endif

             ! ==========================
             ! Tally Target Residual
             ! ==========================
             if ( gsmObj%genModels%fbu%recommendFermiBreakUp(results%targRes%numBaryons) ) then
                ! Fermi Break-Up Driving Routine
                call gsmObj%fermiBreakUpInterface (results%targRes, &
                     & gsmRxn)
             else
                !   Regular decay of target; ap > 12
                wam = one

                !   Determine semiempirical af/an fits for fission cross sections;
                !   gamma energy varies randomly, so need new value for each iteration:
                !   Corrected by NVMokhov, 10/08/03:
                !   NOTE:
                !     afMultiplier = a_f(CEM)/a_f(RAL)
                !     czMultiplier = C(Z)[CEM]/C(Z)[RAL]
                if (proj%particleFlag == bremsProjFlag) then
                   ! Determine use of CEM or LAQGSM values:
                   if( gsmObj%options%useSDCMAfCzMultipliers ) then
                      call setleveldensity( targ%numBaryons, targ%numProtons, &
                           & proj%kinEnergy, targ%afMultiplier, &
                           & targ%czMultiplier )
                   else
                      call setleveldensity( targ%numBaryons, targ%numProtons, &
                           & proj%kinEnergy, targ%afMultiplier, &
                           & targ%czMultiplier, 1 )
                   endif
                endif

                ! Calculating Pre-Equilibrium and Equilibrium Particle Emssion
                call gsmObj%simulateDecay (gsmRxn, &
                     & results%targRes, results%targExc, &
                     & targ%afMultiplier, &
                     & targ%czMultiplier, &
                     & wam )

                if (wam == -13.d0) then
                   wam = one
                   write(gsmObj%io%message,1100) "Target pre-equilibrium decay failure", gsmRxn%eventNum
                   call gsmObj%io%print(3, 3, gsmObj%io%message)
                   go to 30
                endif
             endif

             ! KKG  10/13/04
             if (proj%particleFlag == bremsProjFlag) then
                call proj%brems%incrementTEqv(proj%kinEnergy)
                call proj%brems%incrementSXAbs( &
                     & photonEG%photocrosssection( &
                     &      1000 * proj%kinEnergy, targ%numBaryons))
             endif

          else

             ! Unallowed A, Z values; run Cascade again, unphysical nucleus produced
             write(gsmObj%io%message,1100) "An unphysical target residual was detected", gsmRxn%eventNum
             call gsmObj%io%print(3, 3, gsmObj%io%message)
             go to 30

          endif

          results%targRes%kinEnergy = results%targRes%kinEnergy/1000.d0   ! converting to GeV
          ! if (gsmRxn%eventNum <= nnnp) call gsmObj%printPartProp (proj, targ, intel)
          if (results%info%fusion > zro .and. gsmRxn%output%accountFission > 0) &
               & gsmRxn%outData%nfis = gsmRxn%outData%nfis + 1
          gsmRxn%outData%sfu = gsmRxn%outData%sfu + results%info%wf

          if(results%tallySim) then
             !>>> CRITICAL REGION
             ! Accumulating Events; evaporation and fission model
             call gsmObj%vlobd( gsmRxn )
             !>>> END CRITICAL REGION
          end if
       endif
    endif

    !>>> THE RN_SEED (FROM RN_next_particle) NEEDS TO BE OPENMP PRIVATE
    ! Successfully simulated a single non-elastic event; go to next RN stride
    if ( inelasticEvent .and. gsmObj%useDefaultRNG ) then
       call RN_next_particle( gsmRxn%eventNum+1, 0_int64, gsmRxn%eventNum )
    end if

    return
! ======================================================================
1000 format (A)
1100 format (A, " (event ", i8, ").  Event will be restarted...")
2000 format (2x, 'nc = ', i8)
2050 format (2x, 'nc = ', es11.3)
! ======================================================================
  end subroutine eventLoop
