
  subroutine cascad ( sDCM, clientProj, clientTarg, results)

! ======================================================================
!
!   Drives the main cascade calculation.
!
!   Called by: CEM03
!
!   Calls: PAULIP PINPN PINPN1 POINTE POINTE1 POTEN TYPINT WIM WOPT
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Edited by AJS,  August, 1997.
!    Edited by AJS,  December, 1997.
!    Modified by AJS,  March, 1999.
!    Edited by S. G. Mashnik, LANL, T-2, December, 1998 to use real
!        binding energies for cascade nucleons.
!    Modified by SGM, 2000-2001, to get CEM2k.
!   "Last" change: 12-AUG-2003 by NVMokhov.
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by A. J. Sierk, LANL T-16, January, 2004.
!    Modified by K. K. Gudima, Feb.-Mar., 2004 to include new refrac
!    Modified by A. J. Sierk, LANL T-16, May, 2004.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================
!
!  Definition of partin:
!                       partin(1); Normalized x coordinate of particle
!                       partin(2); Normalized y coordinate of particle
!                       partin(3); Normalized z coordinate of particle
!                       partin(4); sin(theta), direction of momentum
!                       partin(5); cos(theta), direction of momentum
!                       partin(6); sin(phi), direction of momentum
!                       partin(7); cos(phi), direction of momentum
!                       partin(8); kinetic energy of particle
!                       partin(9); rest mass of particle
!
!  Definition of ipatin:
!                       ipatin(1); charge of particle
!                       ipatin(2); (Equals 1 for gamma; 0 otherwise)
!                       ipatin(3); strangeness of particle
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
!  Definition of spt:
!                       spt(1,k) = sin(theta.k)
!                       spt(2,k) = cos(theta.k)
!                       spt(3,k) = kinetic energy of particle k (GeV)
!                       spt(4,k) = charge of particle k
!                       spt(5,k) = rest mass of particle k
!
!  Definition of parz:
!                       parz(1,k) = particle type #; 1-9 for n-pi+
!                       parz(2,k) = kinetic energy of particle k (GeV)
!                       parz(3,k) = theta of particle k
!                       parz(4,k) = phi of particle k
!                       parz(5,k) = index: <100 for cascade particle,
!  negative for hole;   = 100 for preq., or = 1000 for thermal.
!                       parz(6,k) = electric charge of particle k
!
!  enext has units of GeV
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams,    only: zro, hlf, one, two, ten, pi, twpi, &
         & ato3rd, massPiPM, massPi0
    use standardDCMDataClass, only: sDCMData => StandardDCMData, &
         & properlyConstructed
    use standardDCMData,      only: photonEG

    implicit none
    class(StandardDCM),    intent(inout) :: sDCM
    class(sDCMProjectile), intent(inout) :: clientProj
    class(sDCMData),       intent(inout) :: clientTarg
    class(StandardDCMResults), intent(inout) :: results


    integer(int32) :: i, ia, iipar, ipart, ipa, ipsave, irefl, is1, &
         & itemp, numProg, numZones, mv, nabs, nout, np
    real(real64)   :: aac, b1, cutof1, cutof2, cutof3, empi, &
         & fx, obr, p0, pabs, pote, r, r2, rnd, rtest, sabs, sigabs, sign, &
         & sigp, t0mev, t3, temp, temp1, temp2, temp3, temp4, temp5, &
         & tin1, ts, u, w1, w2, w3, w4, w5, wm, zzc
    real(real64),   dimension(9) :: partin = zro, partne = zro
    integer(int32), dimension(5) :: ipatin = zro, ipatne = zro
    real(real64),   dimension(3) :: v= zro, am0 = zro

! ======================================================================

    ! Calculational values local to this simulation of the sDCM class
    type(sDCMPauliInfo) :: pauliData

    !> \todo Check performance with OpenMP/parallelization when
    !>       (a) local, and (b) when allocated. It _might_ be faster
    !>       given the size of the arrays (~18 kB)
    real(real64),   dimension(9, results%maxProgeny+9), target :: pmemo
    integer(int32), dimension(5, results%maxProgeny+9), target :: imemo

! ======================================================================

    ! Verify all types were properly established:
    ! (check target data type)
    if ( clientTarg%constructionState() /= properlyConstructed ) then
       write(sDCM%io%message, 3000)
       call sDCM%io%print(1, 1, sDCM%io%message)
       write(sDCM%io%message, 3500)
       call sDCM%io%print(1, 1, sDCM%io%message)
       results%simState = 10
       return
    end if

    ! (check that passed in results were properly initialized
    if ( .not.results%initialized ) then
       write(sDCM%io%message, 3100)
       call sDCM%io%print(1, 1, sDCM%io%message)
       write(sDCM%io%message, 3500)
       call sDCM%io%print(1, 1, sDCM%io%message)
       results%simState = 10
       return
    end if

    ! (check that the sDCM class was established)
    if ( .not.sDCM%constructed ) then
       write(sDCM%io%message, 3200)
       call sDCM%io%print(1, 1, sDCM%io%message)
       write(sDCM%io%message, 3500)
       call sDCM%io%print(1, 1, sDCM%io%message)
       results%simState = 10
       return       
    end if

    ! Point results object to local procedure's interacting bank
    results%imemo => imemo
    results%pmemo => pmemo


    numZones = clientTarg%numZones()
    obr = clientTarg%coulombPote( numZones ) / two

10  continue

    ! Reset state and flag information of the sDCM:
    results%simState       = 0_int32
    results%numRestarts(:) = 0_int32

    mv = 0     ! index of last storage cell in imemo/pmemo
    irefl = 0  ! Flags the number of times a particle was reflected
    results%numElastic = 0    ! number of elastic events

    ! Reset all pauli data (data that is passed between cascad and pauli)
    pauliData = sDCMPauliInfo()

    results%excitons%numExcProt = 0  ! number of exciton charge
    results%excitons%numExcNeut = 0  ! number of exciton neutrons
    results%excitons%numExcHoles = 0   ! number of electron-hole pairs
    nabs = 0   ! number of absorbed particles
    empi = zro ! gamma energy correction
    ia = nint(clientTarg%numBaryons()) ! initial a value
!   Incorporation of real binding energies for cascade nucleons
!   Corrected to be the Z and A of the compound nucleus, AJS 10/28/03.
    iipar = clientProj%numBaryons + clientProj%numProtons
    aac = clientTarg%numBaryons() + clientProj%numBaryons
    zzc = clientTarg%numProtons() + clientProj%numProtons
    if (clientProj%numBaryons > 0) &
         & call clientTarg%setSepEnergy( sDCM%bindnuc(iipar, aac, zzc) )

!  Determine random entry point of projectile into nucleus
!  KKG 03/04/04
    if (sDCM%options%boundaryEffects == 0)  then
       call sDCM%pinpn (clientProj, clientTarg, partin, ipatin, nout, am0, p0, t3)
    else
       call sDCM%pinpn1 (clientProj, clientTarg, partin, ipatin, nout, am0, p0, t3)
    endif
!  Projectile baryon no:
    temp1 = dble(clientProj%numBaryons)
!  clientProj%gammaFlag = 1 for gamma projectile:
    temp2 = dble(clientProj%gammaFlag)
!  Projectile charge:
    temp3 = dble(clientProj%numProtons)
    if (temp1 == zro) empi = massPiPM*abs(dble(ipatin(1))) + &
         & massPi0*(one - abs(dble(ipatin(1))))
    results%residual%kinEnergy = t3 + temp1*clientTarg%getSepEnergy() + &
         & (one - temp1)*(one - temp2)*(empi + clientTarg%pionPote())
    results%residual%numBaryons = clientTarg%numBaryons() + temp1
    results%residual%numProtons = clientTarg%numProtons() + temp3
!     ainit = results%residual%numBaryons
!     zinit = results%residual%numProtons
    results%residual%linearMom(1) = zro ! x-momenta of projectile within target nucleus
    results%residual%linearMom(2) = zro ! y-momenta
    results%residual%linearMom(3) = p0 ! z-momenta
    results%numProgeny = 0
!  Angular momentum of nucleus + projectile system:
    results%residual%angularMom(1) = am0(1)
    results%residual%angularMom(2) = am0(2)
    results%residual%angularMom(3) = am0(3)
!   KKG: 10/14/04
!   Absorption of gamma with excitation of GDR if energy < 10MeV.
!   Exciton configuration is 1p1h.
    if (clientProj%gammaFlag.ne.0 .and. clientProj%kinEnergy < 0.010) then
       t0mev = clientProj%kinEnergy*1000.d0
       sabs = photonEG%photoCrossSection (t0mev, clientTarg%numBaryons())
       temp = clientTarg%geomCrossSection()
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message,1000) "163"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       pabs = sabs/temp
       b1 = sDCM%rang()
       if (b1 <= pabs) then
          results%excitons%numExcHoles = 1
          rnd = sDCM%rang()
          temp = clientTarg%numBaryons()
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "173"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          if (rnd <= clientTarg%numProtons()/temp) then
             results%excitons%numExcProt = 1
          else
             results%excitons%numExcNeut = 1
          endif
          go to 40
       else
          results%numElastic = 1
          return
       endif
    endif
!  Start of main loop over cascade particles:

20  temp1 = dble(ipatin(4)) ! cascade particle baryon no.
    temp2 = dble(ipatin(1)) ! cascade particle charge no.
! Incorporation of real binding energies for cascade nucleons:
    if (ipatin(4) > 0) then
       ipart = ipatin(4) + ipatin(1)
       if (   results%residual%numProtons >= one .and. &
            & results%residual%numBaryons >= results%residual%numProtons) &
            & call clientTarg%setSepEnergy &
            & ( sDCM%bindnuc (ipart, results%residual%numBaryons, results%residual%numProtons) )
    endif
!  Nuclear zone no:
    itemp = ipatin(5)
! Nuclear Zone potential
    pote = poten (clientTarg, itemp, ipatin)
    cutof1 = temp1*(pote + temp2*obr) + 0.001d0
    if (partin(8) > cutof1) then
!   Determine interaction point, and characteristics of partner.
!   Add POINTE1; KKG 03/04/04
       if (sDCM%options%boundaryEffects == 0)  then
          call sDCM%pointe (clientTarg, partin, ipatin, partne, ipatne, nabs, nout, &
               & v, u, tin1, sigp, sign, sigabs, t3)
       else
          call sDCM%pointe1 (clientTarg, partin, ipatin, partne, ipatne, nabs, nout, &
               & irefl, v, u, tin1, sigp, sign, sigabs, t3)
       endif
       if (ipatin(4) == 0 .and. nabs == 1) then
!   Absorb stopped pion or gamma and create additional excitons:
          results%excitons%numExcProt   = results%excitons%numExcProt + 1
          results%excitons%numExcNeut   = results%excitons%numExcNeut + 1
          results%excitons%numExcHoles  = results%excitons%numExcHoles  + 2
          go to 40
       endif
!   If cascade particle had a problem in GEOM8, go back, discard and
!   repeat this cascade (Very rare! AJS 3/29/99).
       if (ipatin(1) == 2) then
          results%numRestarts(1) = results%numRestarts(1) + 1
          go to 10
       endif
       if (nabs > 0 .and. ipatin(4) == 1) then
!   If nucleon stays in POINTE for many iterations, assume it is
!   absorbed.
          results%excitons%numExcProt = results%excitons%numExcProt + ipatin(1)
          results%excitons%numExcNeut = results%excitons%numExcNeut + 1 - ipatin(1)
       else
!   Few nucleon interactions
          if (nout <= 0) then
!   Interaction of cascade particle happens:
!
!  New smoother transition method (CEM03) by AJS, April 27, 2004:
!  This modification will be subject to any data being
!  measured between 120 and 230 MeV incident nucleon energy!!
             if (clientProj%kinEnergy > 0.225d0) then
                ipa = 1
                go to 30
             elseif (clientProj%kinEnergy > 0.075d0) then
                fx = hlf*(one + cos((ten*clientProj%kinEnergy - 0.75d0)*pi/1.5d0))
                b1 = sDCM%rang()
                if (b1 > fx) then
                   ipa = 1
                   go to 30
                endif
             endif
!   Original "Proximity" method for determining exit to prequilibrium:
             if (ipatin(4) == 0) then
!   Meson or gamma:
                ipa = 1
             else
!   Baryon:
                cutof3 = pote + 0.051d0
                if (partin(8) <= cutof3) then
                   r2 = partin(1)**2 + partin(2)**2 + partin(3)**2
                   temp = r2
                   if (temp < 0.0d0) then
                      temp = 0.01d0
                      write(sDCM%io%message, 1100) "260"
                      call sDCM%io%print(4, 3, sDCM%io%message)
                   end if
                   r = clientTarg%zoneBoundR( numZones ) * &
                        & sqrt(temp) ! Approximate nucleus radius
!   Asymptotic kinetic energy for particle when outside nucleus:
                   ts = partin(8) - pote
                   is1 = ipatin(1) ! ejectile charge
                   rtest = clientTarg%zoneBoundR( numZones ) - &
                        & r - 0.001d0
                   if (rtest < zro) then
!   r >~ clientTarg%zoneBoundR(n):
                      ipa = 1
                   else
!   If kinetic energy of particle inside the nucleus < 50 MeV + V;
!   average imaginary potential over Fermi distribution of nucleons
!   inside nucleus to determine changeover to preequilibrium decay.
                      wm = zro
                      do i = 1,50
                         w1 = sDCM%wim (clientTarg, 1, partin, ipatin)
                         w2 = sDCM%wim (clientTarg, 0, partin, ipatin)
                         wm = wm + w1 + w2
                      end do
!  0.01*wm is approximate "cascade" imaginary optical potential.
!   w3 is "scattering phenomenology" imaginary optical potential.
                      w3 = sDCM%wopt (clientTarg, ts, r, is1)
                      w4 = 0.01d0*wm - w3
                      temp = w3
                      if (temp < div0Lim .and. temp > -div0Lim) then
                         temp = div0Lim
                         write(sDCM%io%message,1000) "287"
                         call sDCM%io%print(4, 3, sDCM%io%message)
                      end if
                      w5 = abs(w4/temp)
                      if (w5 >= 0.3d0) then
                         ipa = 0
                      else
                         ipa = 1
                      endif
                   endif
                else
!   Kinetic energy > cutof3:
                   ipa = 1
                endif
             endif
!   Calculate the type and results of the interaction:
30           call sDCM%typint (clientTarg, partin, ipatin, partne, ipatne, &
                  & v, u, tin1, sigp, sign, sigabs, mv, np, results)
             if (np <= 0) then
!   Particle is absorbed into nucleus; add an exciton:
                results%excitons%numExcProt = results%excitons%numExcProt + ipatin(1)
                results%excitons%numExcNeut = results%excitons%numExcNeut + 1 - ipatin(1)
                go to 40
             elseif (np > 20) then
!   Rerun cascade if np is strange:
                results%numRestarts(2) = results%numRestarts(2) + 1
                go to 10
             else
                t3 = sqrt(results%pmemo(1,mv+3)**2 + results%pmemo(2,mv+3)**2 + &
                     & results%pmemo(3,mv+3)**2 + results%pmemo(9,mv+3)**2) - &
                     & results%pmemo(9,mv+3)
!   Determine if scattering is allowed by exclusion princlple:
                ipsave = ipatin(2)
                call sDCM%paulip (clientTarg, partin, ipatin, v, &
                     & mv, np, irefl, ipa, pauliData, results)
! FCG corrects error of rejected gamma path, by FCG August 2000:
                if (ipsave.ne.0 .and. ipatin(2).ne.0) t3 = results%residual%kinEnergy
                if (ipa == 2) then
!   Scattering allowed; add an exciton:
                   results%excitons%numExcProt = results%excitons%numExcProt + ipatin(1)
                   results%excitons%numExcNeut = results%excitons%numExcNeut + 1 - ipatin(1)
                   go to 40
                else
                   go to 20
                endif
             endif
          endif
          temp1 = dble(ipatin(4))
          itemp = ipatin(5)
          if (ipatin(1) <= 0) then
             temp2 = zro
          else
             temp2 = dble(ipatin(1))
          endif
          cutof2 = poten (clientTarg, itemp, ipatin) + &
               & (temp2*obr) + 0.001d0
          if (partin(8) <= cutof2) then
!   If energy is below cutof2; absorb particle; create exciton:
             results%excitons%numExcProt = results%excitons%numExcProt + ipatin(1)
             results%excitons%numExcNeut = results%excitons%numExcNeut + 1 - ipatin(1)
          else
             ! Particle no longer interacts (i.e. secondary particle) - update residual and store progeny


             ! Increment number of progeny created
             results%numProgeny = results%numProgeny + 1

             ! Check if progeny array is at capacity
             if ( results%numProgeny > results%maxProgeny ) then

                ! Array of progeny will be exceeded with the next particle; warn and stop sDCM
                write(sDCM%io%message, 2000) results%numProgeny
                call sDCM%io%print(3, 3, sDCM%io%message)
                write(sDCM%io%message, 2010)
                call sDCM%io%print(3, 3, sDCM%io%message)

                ! Set number of progeny to maximum allowed
                results%numProgeny = results%maxProgeny

                ! Flag that no more particles exist in the bank of particles that continue to interact
                mv = 0

             else
                ! Progeny array is not at capacity; can store progeny and update residual

                ! Update residual information
                ! (update A, Z)
                temp1 = dble( ipatin(4) )   ! Baryon number of compound nucleus (i.e. pion)
                temp2 = dble( ipatin(1) )   ! Number of protons in progeny
                results%residual%numBaryons = results%residual%numBaryons - temp1
                results%residual%numProtons = results%residual%numProtons - temp2
                if (   results%residual%numBaryons < 4.d0 .or. &
                     & results%residual%numProtons < 1.d0 .or. &
                     & results%residual%numProtons > results%residual%numBaryons) then
                   ! An unphysical residual was created - flag error and return
                   ! NOTE: This only matters when Z>A; otherwise, could be highly unstable
                   !       compound that will undergo extreme preeq./evap./Fermi Breakup
                   !       decay/deexcitation
                   write(sDCM%io%message, 2100) results%residual%numBaryons, &
                        & results%residual%numProtons
                   call sDCM%io%print(5, 3, sDCM%io%message)
                   results%simState = 1

                   ! Revert A/Z pair of compound to what it was prior to reducing based on progeny
                   results%residual%numBaryons = results%residual%numBaryons + temp1
                   results%residual%numProtons = results%residual%numProtons + temp2
                   results%numProgeny = results%numProgeny - 1
                   mv = 0   ! Claim that bank of further interacting particles is empty (has no meaning now)
                   go to 40

                end if

                ! (update KE)
                temp3 = ipatin(2)
                if (temp1 == zro) then
                   temp5 = dble( ipatin(1) )
                   empi = massPiPM*abs(temp5) + massPi0*(one - abs(temp5))
                endif
                results%residual%kinEnergy = results%residual%kinEnergy - &
                     & (partin(8) + temp1*clientTarg%getSepEnergy() + &
                     & (one - temp1)*(empi + clientTarg%pionPote())*(one - temp3))
                if (results%residual%kinEnergy  <=  0.0001d0) then
                   !  E* <= 100 keV and Z & A unchanged; call it elastic scattering.
                   if (   results%residual%numBaryons == clientTarg%numBaryons() .and. &
                     & results%residual%numProtons == clientTarg%numProtons()) then
                      ! Flag elastic scattering occurred
                      results%numElastic = 1

                      ! Flag that no progeny was created (elastic event instead)
                      results%numProgeny = results%numProgeny - 1
                      return
                   endif
                endif

                ! (update linear momentum)
                !   |p| = m*gamma*v
                temp = partin(8)*(partin(8) + two*partin(9))
                if (temp < 0.0d0) then
                   temp = 0.01d0
                   write(sDCM%io%message, 1100) "360"
                   call sDCM%io%print(4, 3, sDCM%io%message)
                end if
                temp3 = sqrt(temp) ! total momentum
                temp4 = partin(4)*partin(7)   !   p-x/|p|
                temp5 = partin(4)*partin(6)   !   p-y/|p|
                results%residual%linearMom(1) = results%residual%linearMom(1) - &
                     & temp3*temp4 ! x-momenta
                results%residual%linearMom(2) = results%residual%linearMom(2) - &
                     & temp3*temp5 ! y-momenta
                results%residual%linearMom(3) = results%residual%linearMom(3) - &
                  & temp3*partin(5) ! z-momenta

                ! (update angular momentum)
                !  L-x: angularMom(1) = angularMom(1) - p * ( y*(pz/p)- z*(py/p) )
                results%residual%angularMom(1) = results%residual%angularMom(1) + &
                     & temp3*(partin(3)*temp5 - partin(2)*partin(5))
                results%residual%angularMom(2) = results%residual%angularMom(2) + &
                     & temp3*(partin(1)*partin(5) - partin(3)*temp4)
                results%residual%angularMom(3) = results%residual%angularMom(3) + &
                     & temp3*(partin(2)*temp4 - partin(1)*temp5)


                ! Store progeny information for the banked secondary particle
                numProg = results%numProgeny
                results%progenyBnk(numProg)%xCoord = clientTarg%zoneBoundR( numZones ) &
                     & * partin(1)
                results%progenyBnk(numProg)%yCoord = clientTarg%zoneBoundR( numZones ) &
                     & * partin(2)
                results%progenyBnk(numProg)%zCoord = clientTarg%zoneBoundR( numZones ) &
                     & * partin(3)
                results%progenyBnk(numProg)%sinTheta    = partin(4)
                results%progenyBnk(numProg)%cosTheta    = partin(5)
                results%progenyBnk(numProg)%sinPhi      = partin(6)
                results%progenyBnk(numProg)%cosPhi      = partin(7)
                results%progenyBnk(numProg)%kinEnergy   = partin(8)
                results%progenyBnk(numProg)%restMass    = partin(9)
                results%progenyBnk(numProg)%numProtons  = ipatin(1)
                results%progenyBnk(numProg)%numBaryons  = ipatin(4)
                results%progenyBnk(numProg)%strangeness = ipatin(3)
                results%progenyBnk(numProg)%photonFlag  = ipatin(2)
                results%progenyBnk(numProg)%nuclearZone = ipatin(5)
                results%progenyBnk(numProg)%index       = pauliData%ing
                if ( results%progenyBnk(numProg)%numBaryons == 1 .and. &
                     & pauliData%indi == 3 ) then
                   results%progenyBnk(numProg)%index = -results%progenyBnk(numProg)%index
                end if

                pauliData%igs(numProg) = pauliData%ing

             end if
          endif
       endif
    else
!   Kinetic energy is below cutof1; add a particle exciton (absorbed):
       results%excitons%numExcProt = results%excitons%numExcProt + ipatin(1)
       results%excitons%numExcNeut = results%excitons%numExcNeut - ipatin(1) + 1
    endif


    ! Done with last primary particle
40  if (mv > 0) then
       ! Bank is non-zero; obtain banked particle and treat as primary particle
       do i = 1,9
          partin(i) = results%pmemo(i,mv)
       end do
       do i = 1,5
          ipatin(i) = results%imemo(i,mv)
       end do
       pauliData%ing = pauliData%ngen(mv)
       mv = mv - 1 ! decrement number of banked particles
       nabs = 0 ! resetting number of absorbed particles
       nout = 0
       irefl = 0 ! particle not reflected
       go to 20 ! re-running cascade
    else
!   Convert angular momentum to units of h-bar.  Constant =
!   1000./(h-bar*c = 197.32858)
       temp1 = 5.0676896d0*clientTarg%zoneBoundR( numZones )
       results%residual%angularMom(1) = temp1*results%residual%angularMom(1)
       results%residual%angularMom(2) = temp1*results%residual%angularMom(2)
       results%residual%angularMom(3) = temp1*results%residual%angularMom(3)
    endif
    return

! ======================================================================
! (math warnings)
1000 format("Divide by zero error prevented in 'cascad.f90' line(s) ", A)
1100 format("Square root error prevented in 'cascad.f90' line(s) ", A)
! (general warnings)
2000 format(i5, " progeny have been created in the standard ", &
          & "DCM simulation.")
2010 format("   The progeny array bank was exceeded.")
2100 format("A light or invalid residual nucleus (A=", f5.1, ", Z=", f5.1, &
          & ") was created during the Standard DCM simulation.")
3000 format("The target object used by the Standard DCM class was ", &
          & "not properly constructed.")
3100 format("The results object used by the Standard DCM class was ", &
          & "not properly constructed.")
3200 format("The Standard DCM class was NOT constructed.")
3500 format("   Unable to further simulate INC physics.")
! ======================================================================
  end subroutine cascad
