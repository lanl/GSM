
  subroutine peqemt (preeqObj, calcVars, exciton, lm, ipflg, results)

! ======================================================================
!
!    Preequilibrium emission calculation; extracted from old PRECOF
!    routine, 10/08/03.
!
!   Calls: AUXL GAMAGU2 KINEMA MOLNIX ROTOR TKINM3 TRANS8 
!
!    Written by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by K. K. Gudima, June, 2004.  
!    Edited by AJS, January, 2005.
!    Corrected prequilibrium angular distribution, S. G. Mashnik
!    July, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!    Modified by LMK, July 2012 - July 2013 (expanded preeq).
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Modified by LMK, 2014-2015, for expansion of preequilibrium to heavy ions
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use preequilibriumParams, only: zro, one, two, thr, thousand, thrd, twpi, &
         & emnuct, emnucb, electron_mass
    use preequilibriumData, only: t0LogicSwitch, aj, ajthr, zj, emured

! LMK 02/2013
!    use gammaJClass, only : getGammaBeta    (migrated to a passed-in class object)
    use lambda_j, only: gamagu3, kin_energy

    implicit none
    class(Preequilibrium),  intent(inout) :: preeqObj
    type(preequilibriumCalculation), intent(inout) :: calcVars
    type(preeqExcitonData), intent(inout) :: exciton
    integer(int32),         intent(inout) :: lm
    integer(int32),         intent(  out) :: ipflg
    type(preequilibriumResults), intent(inout) :: results

    integer(int32) :: iaa, inj, ip, izj, j, l, n1
    real(real64)   :: ac0, almax, alp, alpx, alpy, alpz, ax, b, b1, b2, &
         & b3, bz, bz1, c, c1, c11, c2, c3, cf, cfr, ct, ct2, ctr, ctrp, &
         & c21, delu, e, eb, emx, emxp, ep1, ep3, f1, f2, fc, fi, fi1, &
         & fi2, fir, g, p2, pm, ppp, sf, sfr, st, st2, str, tl, un, v2
    logical        :: emitok

! ======================================================================

    real(real64), parameter, dimension(2) :: gb1 = [one, one]
    real(real64), dimension(66) :: gj = zro
    real(real64), dimension(3) :: p12 = zro, pl = zro, ps = zro

! ======================================================================

    ! Set gb() for neutron and proton
    calcVars%gb(1) = gb1(1)
    calcVars%gb(2) = gb1(2)


10  if (exciton%numNucleons < one) then
       exciton%numNucleons = exciton%numNucleons + one
       exciton%numHoles = exciton%numHoles + one
       exciton%numTotal = exciton%numTotal + 2
       go to 10
    endif

    calcVars%exn = exciton%numNucleons + exciton%numHoles
    ac0 = calcVars%ac
    ip = nint(exciton%numNucleons)
    iaa = nint(results%residual%numBaryons)

    ! Set up alj(:) array
    call setALJ(calcVars, exciton%numTotal, ip, exciton%numNucleons, calcVars%exn)


! Not needed with new gamma_j model, LMK 04/2015
!   KKG 03/04/04
!   Empirical multipliers for complex particle emission:
!   MIB, KKG 06/21/04
!      gb(3) = gb0(iaa,1)*gbm(1)
!      gb(4) = gb0(iaa,2)*gbm(2)
!      gb(5) = gb0(iaa,2)*gbm(3)
!      gb(6) = gb0(iaa,3)*gbm(4)

! Calculate gamma_beta coefficients (LMK, 04/2015)
    call preeqObj%preeqData%gammaJObj%getGammaBeta( &
         & preeqObj%options%numPreeqType,results%residual%numBaryons, &
         & calcVars%gb )

    g = zro
    do l = 1,preeqObj%options%numPreeqType 
!   gj is the emission rate of particles of type j (n,p,d,t,3He,4He...)
!   into the continuum.
       emitok = exciton%numNucleons >= aj(l) .and. calcVars%exn >= (aj(l)+one) .and. &
            & exciton%numProtons >= zj(l) .and. &
            & (l.ne.1 .or. exciton%numNucleons-exciton%numProtons >= one) .and. &
            & calcVars%rj(l) > zro
       gj(l) = zro
       if (emitok) then
          calcVars%ac = 0.595d0*calcVars%ami(l)
          call gamagu3(gj(l),l,calcVars%uej(l), &
               & calcVars%bj(l), calcVars%vj(l), &
               & calcVars%exn, aj(l), zj(l), calcVars%gb(l), &
               & preeqObj%options%r0Mult, calcVars%ac, &
               & calcVars%zfj(l), calcVars%afj(l), &
               & calcVars%alj(l) )
          gj(l) = calcVars%redpre(l)*gj(l)
          g = g + gj(l)   ! g is total preequilibrium particle emission rate
       endif
    end do

    if (exciton%numTotal <= 0) return
    calcVars%ac = ac0
!   Find rates for changing exciton number with no emission:
    call preeqObj%trans8 (exciton%numNucleons, exciton%numHoles, calcVars%ac, c1, c2, c3, results)
    c = c1 + c2 + c3
    if (c <= zro) then
       ipflg = 1
       return
    endif
    b2 = c + g
    if (c < div0Lim .and. c > -div0Lim) then
       c = div0Lim
       write(preeqObj%io%message,1000) "238, 239"
       call preeqObj%io%print( 4, 3, preeqObj%io%message)
    end if
    c11 = c1/c
    c21 = (c2 + c1)/c
20  b1 = preeqObj%rng()
    if (b2 < div0Lim .and. b2 > -div0Lim) then
       b2 = div0Lim
       write(preeqObj%io%message,1000) "244"
       call preeqObj%io%print( 4, 3, preeqObj%io%message)
    end if
    if (b1 > g/b2) then
!  Transition in exciton number (-2, 0, or +2)
       b3 = preeqObj%rng()
       if (b3 <= c21) then
          bz = preeqObj%rng()
          if (results%residual%numBaryons < div0Lim .and. results%residual%numBaryons > -div0Lim) then
             results%residual%numBaryons = div0Lim
             write(preeqObj%io%message,1000) "253"
             call preeqObj%io%print( 4, 3, preeqObj%io%message)
          end if
          bz1 = results%residual%numProtons/results%residual%numBaryons
          if (b3 <= c11) then
             exciton%numTotal = exciton%numTotal + 2
             exciton%numNucleons = exciton%numNucleons + one
             exciton%numHoles = exciton%numHoles + one
             fc = one
          elseif (b3 > c11) then
             exciton%numTotal = exciton%numTotal - 2
             exciton%numNucleons = exciton%numNucleons - one
             exciton%numHoles = exciton%numHoles - one
             fc = -one
             if (exciton%numProtons <= zro) then
                ipflg = 1
                return
             endif
          endif
          if (bz <= bz1) exciton%numProtons = exciton%numProtons + fc
          ipflg = 1
          return
       else
!  delta-n = 0 option; do not recalculate c's and g's!
          go to 20
       endif
    else
!   Preequilibrium particle is emitted:
       ppp = exciton%numNucleons
!   gj is converted to the sum up to j of all the rates for
!   particles with index less than or equal to j.
       do j = 2,preeqObj%options%numPreeqType           !LMK 07/2012
          gj(j) = gj(j-1) + gj(j)
       end do
!   Randomly select particle to be emitted, according to partial
!   emission rates.
       b = preeqObj%rng()*g
       do j = 1,preeqObj%options%numPreeqType           !LMK 07/2012
          if (b <= gj(j)) then
             lm = j
             go to 30
          endif
       end do
30     continue
       call kin_energy(lm, ep1, calcVars%uej(lm), &
            & calcVars%bj(lm), calcVars%vj(lm), &
            & calcVars%exn, aj(lm), zj(lm), calcVars%gb(lm), &
            & preeqObj%options%r0Mult, calcVars%ac, &
            & calcVars%zfj(lm), calcVars%afj(lm), &
            & calcVars%alj(lm), preeqObj%rng)
       ! ep1 is the KE of emitted particle

       izj = nint(zj(lm))
       inj = nint(aj(lm) - zj(lm))
       if (izj.ne.0) then
          emxp = preeqObj%molEnergy%defineEnergy (izj, inj, 2)
       else
          emxp = 8.071d0
       endif
       ep3 = thousand*emnuct*aj(lm) + emxp - zj(lm) * electron_mass


       ! Update exciton information
       exciton%numNucleons = exciton%numNucleons - aj(lm)
       n1 = nint(aj(lm))
       exciton%numTotal = exciton%numTotal - n1
       exciton%numProtons = exciton%numProtons - zj(lm)


       ! Update residual information
       results%residual%numBaryons  =   calcVars%afj(lm)
       calcVars%athrd   =   calcVars%afjthr(lm)
       results%residual%numProtons  =   calcVars%zfj(lm)
       un   =   results%residual%numBaryons - results%residual%numProtons
       calcVars%iz   =   nint(results%residual%numProtons)
       calcVars%in   =   nint(un)
       calcVars%iz   =   max(1,calcVars%iz)
       calcVars%in   =   max(1,calcVars%in)


       ! Obtain new mass excess and total speed
       emx   =   preeqObj%molEnergy%defineEnergy ( calcVars%iz, calcVars%in, 2 )
       v2   =   results%residual%normSpeed(1)**2 + results%residual%normSpeed(2)**2 + &
            & results%residual%normSpeed(3)**2


       ! Obtain momentum
       pm   =   sqrt(abs(ep1*(ep1 + two*ep3)))


! Randomly samply theta-angle
       b1 = preeqObj%rng()
       !  KKG 06/21/04; CTKALB moved inline by AJS 01/06/05
       !  Use old CEM logic for incident energies above 210 MeV. SGM 7/12/05
       !  Yes, it's a kluge, but time is finite.  AJS
       if (preeqObj%preeqData%kinEnergy() <= t0LogicSwitch) then
          eb = ep1 + calcVars%bj(lm)
          ax = 0.04d0*eb + 1.8d-6*eb**3
          if (ax < div0Lim .and. ax > -div0Lim) then
             ax = div0Lim
             write(preeqObj%io%message,1000) "348"
             call preeqObj%io%print( 4, 3, preeqObj%io%message)
          end if
          ct = -one + log(one + b1*(exp(two*ax) - one))/ax
       else 
          !    Uniform distribution on 0-->pi/2:
          ct = b1                     
          !   Momentum of 250 MeV/c corresponds to about 34 MeV Fermi energy.
          p2 = 250.d0*preeqObj%rng()**thrd
          ct2 = one - two*preeqObj%rng()
          st2 = sqrt(one - ct2*ct2)
          fi2 = twpi*preeqObj%rng()
          !   Factor of 1000. to convert pnx from GeV/c to MeV/c.
          if (ppp< div0Lim .and. ppp > -div0Lim) then
             ppp = div0Lim
             write(preeqObj%io%message,1000) "362"
             call preeqObj%io%print( 4, 3, preeqObj%io%message)
          end if
          f1 = thousand/ppp
          f2 = p2*st2
          p12(1) = f2*cos(fi2) + results%residual%linearMomX*f1
          p12(2) = f2*sin(fi2) + results%residual%linearMomY*f1
          p12(3) = p2*ct2 + results%residual%linearMomZ*f1
       endif
       st = sqrt(abs(one - ct*ct))
       fi = twpi*preeqObj%rng()
       cf = cos(fi)
       sf = sin(fi)
       !   pl is the momentum of the particle with respect to the source:
       pl(1) = pm*st*cf
       pl(2) = pm*st*sf
       pl(3) = pm*ct
       if (v2 > 1.0d-7) then
          !  Nucleus has non-zero velocity; do appropriate kinematic
          !  transformation on preequlibrium particle:
          !   p12 is the momentum of the emitting source in the lab frame:
          ! KKG 06/21/04;  SGM 07/12/05:
          if (preeqObj%preeqData%kinEnergy() <= t0LogicSwitch) then
             ps(1) = pl(1) 
             ps(2) = pl(2) 
             ps(3) = pl(3) 
          else
             call preeqObj%rotor (p12, results%residual%normSpeed, pl, ps)
          endif
          call preeqObj%kinema (ps, results%residual%normSpeed, pl, ct, st, &
               & cf, sf, tl, ep3)
!  After kinema, pl is the lab momentum of the emitted particle:
       else
!  Nuclear velocity essentially zero;
!  Isotropic emission of preequilibrium particle:
          tl = ep1
       endif
       results%residual%linearMomX = results%residual%linearMomX - pl(1)/thousand 
       results%residual%linearMomY = results%residual%linearMomY - pl(2)/thousand 
       results%residual%linearMomZ = results%residual%linearMomZ - pl(3)/thousand 
       if (calcVars%iz > 7 .and. calcVars%in > 7) then
          results%residual%restMass = results%residual%numBaryons * emnucb + &
               & emx/thousand
       else
          results%residual%restMass = results%residual%numBaryons * emnuct + &
               & emx/thousand
       endif
       e = sqrt(results%residual%linearMomX**2 + results%residual%linearMomY**2 + &
            & results%residual%linearMomZ**2 + results%residual%restMass**2)
       results%residual%recEnergy = (e - results%residual%restMass)*thousand 
       if (e < div0Lim .and. e > -div0Lim) then
          e = div0Lim
          write(preeqObj%io%message,1000) "412-414"
          call preeqObj%io%print( 4, 3, preeqObj%io%message)
       end if
       results%residual%normSpeed(1) = results%residual%linearMomX/e
       results%residual%normSpeed(2) = results%residual%linearMomY/e
       results%residual%normSpeed(3) = results%residual%linearMomZ/e
       fi1 = atan2 (sf, cf)
       if (fi1 < zro) fi1 = twpi + fi1
       almax = 0.219327d0 * preeqObj%options%r0Mult *  &
            & ( calcVars%athrd + ajthr(lm) ) * &
            & sqrt(  abs( emured(iaa,n1) * (ep1 - calcVars%vj(lm)) )  )
       alp = almax*sqrt(  abs( preeqObj%rng() )  )
40     continue
       ctr = one - two*preeqObj%rng()
       str = sqrt(abs(one - ctr*ctr))
       fir = twpi*preeqObj%rng()
       cfr = cos(fir)
       sfr = sin(fir)
       ctrp = ctr*ct + str*st*(cfr*cf + sfr*sf)
!   Exclude emission inside of parent nucleus:
       if (ctrp < zro) go to 40
       alpx = alp*(str*sfr*ct - st*sf*ctr)
       alpy = alp*(st*cf*ctr  - str*cfr*ct)
       alpz = alp*str*st*(cfr*sf - cf*sfr) 


       ! Reduce angular momentum
       results%residual%angularMom(1) = results%residual%angularMom(1) - alpx
       results%residual%angularMom(2) = results%residual%angularMom(2) - alpy
       results%residual%angularMom(3) = results%residual%angularMom(3) - alpz

       ! Re-obtain angular momentum (modify is needed), fission barrier, and angular momentum quantum number
       call preeqObj%preeqData%auxl (results%residual%numBaryons, &
            & results%residual%numProtons, results%residual%angularMom, &
            & results%residual%fissBarr, results%residual%angMomFlag, &
            & results%residual%rotEnergy, delu)

       ! Update kinetic energy of residual
       results%residual%kinEnergy = calcVars%uej(lm) + delu - &
            & calcVars%bj(lm) + calcVars%pevapj(lm) - &
            & ep1 + results%residual%rotEnergy

       ! Store emitted fragment information
       results%numProgeny = results%numProgeny + 1
       results%progenyBnk(results%numProgeny)%numBaryons  = aj(lm)         ! Mass Number [Num. nucleons, A]
       results%progenyBnk(results%numProgeny)%numProtons  = zj(lm)         ! Atomic Number [Charge, Z]
       results%progenyBnk(results%numProgeny)%kinEnergy   = tl             ! Kinetic Energy [MeV]
       results%progenyBnk(results%numProgeny)%restMass    = ep3            ! Rest mass [MeV/c**2]
       results%progenyBnk(results%numProgeny)%theta       = atan2 (st, ct) ! Theta of momentum vector
       results%progenyBnk(results%numProgeny)%phi         = fi1            ! Phi   of momentum vector
       results%progenyBnk(results%numProgeny)%origin      = 0.0_real64     ! Fragment origin (=0 from preeq. emission; =1 for Fermi Break-up)

       ! Check for photon emission:
       if ( preeqObj%usePhotonEmission ) then
          call preeqObj%photonEmission( results%progenyBnk(results%numProgeny), results%residual )
       end if


    endif

    return
! ======================================================================
1000 format("Divide by zero error prevented in ", &
          & "'peqemt.f90', line ", A)
! ======================================================================
  end subroutine peqemt
