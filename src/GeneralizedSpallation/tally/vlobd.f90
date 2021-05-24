
  subroutine vlobd (gsmObj, gsmRxn)

! ======================================================================
!
!   Accumulate events in 218 particle-emission output channels
!   and in spectrum and angular distribution arrays.
!
!   All variables referring to coalescence products have a "co" in them.
!   Added separate arrays for spallation residuel particle emission and
!   spectra, fission-fragment evaporation, and prefission spectra.
!   (AJS 09/22/03)
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Modified by A. J. Sierk  May-June, 1996.
!   Modified by SGM at 06/16/2000
!   Modified by KKG in 10/2001 , 2/2003, and 9/2003.
!   "Last change": 13-AUG-2003 by NVMokhov
!   Modified by AJS for "fission", "spallation" and prefission  arrays;
!       9/2003.
!   Modified by AJS to correct the calculation of the mass number of
!       the residual nucleus, 3/5/04.
!  [Modified to tally data for fission events only when idel = 2,
!       AJS, 4/16/04.]  OBSOLETE
!   Modified to tally data for fission events only when mpyld = 2,
!       AJS, 5/13/05.
!   Modified to tally data for different stages of the reaction when
!       misy = 3 by K. K. Gudima, November-December, 2004.
!   Edited by AJS, January, 2005.
!   Removed printout of warning for fission fragments, AJS, June, 2005.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!   Edited by CMJ, XCP-3, 2016-2017, (LAQGSM expansion)
!
! ======================================================================
!
!  Definition of spt:
!                       spt(1,k) = sin(theta.k)
!                       spt(2,k) = cos(theta.k)
!                       spt(3,k) = kinetic energy of particle k (GeV)
!                       spt(4,k) = charge of particle k
!                       spt(5,k) = rest mass of particle k
!  Definition of parz:
!                       parz(1,k) = particle type #; 1-9 for n-pi+;
!                                   A + 999*Z for larger fragments and
!                                   residual nuclei. (1000Z + N)
!                       parz(2,k) = kinetic energy of particle k (GeV)
!                       parz(3,k) = theta of particle k
!                       parz(4,k) = phi of particle k
!                       parz(5,k) = index: <100 for cascade,
!                       negative for hole;
!                       = 100 for preq., or
!                       = 200 for coalescence.
!                       = 1000 for evaporation from spallation residue
!                              (and for residue itself).
!                       = 1500 for Fermi breakup.
!                       = 2000 for evaporation from fission product
!                              (and for fragments themselves).
!                       parz(6,k) = electric charge of particle k
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, int64, real64, error_unit
    use gsm_params, only: zro, one, two, thr, four, fiv, ten, thousand, &
         & radiantodegree
    use evaporationfissiondata, only: ifa, ifz

    implicit none
    class(GSM),        intent(inout) :: gsmObj
    type(GSMReaction), intent(inout), target :: gsmRxn

    integer(int32) :: i, ic, icp, ida, idz, ip, ipar, ipoc, iq, it, itet, &
         & iu, jt, ju, m, nd, nhe3, nhe4, nn, nnc, nne, nnp, nnpof, &
         & nnprf, np, npimin, npiplu, npizer, nt
    real(real64)   :: abeg, am, ama, amaxx, asum, da, dz, eaz, paz, temp, &
         & tet, totas, totzs, tt, vz, ytet, zbeg, zm, zsum
    logical :: doit, fisevent, itt, test

    type(GSMResults),    pointer :: results => NULL()
    type(GSMProjectile), pointer :: proj => NULL()
    type(GSMTarget),     pointer :: targ => NULL()
    type(OutputData),    pointer :: outData => NULL()

! ======================================================================

    real(real64) :: te1, te2, dtt, se, dtdo, d2spec, d2spe
    integer(int32) :: ntt, ntet, nt2, nt3, nti
    common /d2sdto/  te1(200), te2(200), dtt(200), &
         & se(200), dtdo(200,10), d2spec(9,10,200), &
         & d2spe(14,10,200), ntt, ntet, nt2, nt3, nti(4)
    real(real64) :: speco, angco, d2speco
    common /rescoa/  speco(4,200), angco(4,20), d2speco(4,10,200)
    real(real64) :: spef, angf, d2spef
    common /resfis/  spef(6,200), angf(6,20), d2spef(6,10,200)
    real(real64) :: spe, an
    common /resmac/  spe(14,200), an(14,20)
    real(real64) :: specsp, angsp, d2spesp
    common /respal/  specsp(6,200), angsp(6,20), d2spesp(6,10,200)
    real(real64) :: specpf, angpf, d2spepf
    common /resprf/  specpf(6,200), angpf(6,20), d2spepf(6,10,200)
    real(real64) :: spec, ang, chan, dadz
    common /result/  spec(9,200), ang(9,20), chan(218), dadz(351,151)
!   KKG 11/13/04
    real(real64) :: dadz4
    common /result1/ dadz4(4,351,151)

    real(real64) :: totke
    common /tkftot/  totke

! ======================================================================

    ! When using OpenMP, ensure critical region:
    !$OMP critical

    ! Associate quick-access variables
    results => gsmRxn%results
    proj    => results%initialProj
    targ    => results%initialTarg
    outData => gsmRxn%outData

    zbeg = targ%numProtons + fiv
    abeg = targ%numBaryons + fiv
    totas = zro
    totzs = zro
    amaxx = zro
    do i = 1, gsmObj%options%evapOpts%numEvapType
        amaxx = max (amaxx, dble(ifa(i)))
    end do
    fisevent = (results%info%fusion == one)

    if (gsmRxn%output%energySpectra > 0 .or. &
         & gsmRxn%output%angularSpectra > 0 .or. &
         & gsmRxn%output%doubleDiffSpectra > 0 .or. &
         & gsmRxn%output%channelCrossSection > 0 .or. &
         & gsmRxn%output%multiplicities > 0 .or. &
         & gsmRxn%output%nuclideCrossSection > 0) then
       nn = 0
       np = 0
       nd = 0
       nt = 0
       nhe3 = 0
       nhe4 = 0
       npimin = 0
       npizer = 0
       npiplu = 0
       nnc = 0
       nnp = 0
       nne = 0
       nnprf = 0
       nnpof = 0
!  kkg 04/15/2015
       asum=zro
       zsum=zro
!
       do m = 1,results%numProgeny
          am = zro
          zm = zro
          if (results%progenyBnk(m)%restMass < 0.001d0) go to 40
          ipar = nint(results%progenyBnk(m)%typeID) ! ejectile type
          zm = results%progenyBnk(m)%numProtons !ejectile charge

!  Assigning A-value for ejectile; A = 0 for ipar = 7-9 (pions)
          if (ipar < 3) then
             am = one ! for n, p
          elseif (ipar < 4) then
             am = two ! for d
          elseif (ipar < 6) then
             am = thr ! for t and he-3
          elseif (ipar == 6) then
             am = four ! for he-4
          elseif (ipar > 9) then
!   Modified calculation of mass number of residual nucleus:
             am = results%progenyBnk(m)%typeID - 999.d0*zm  !mass not conserved
          endif
!  kkg 04/15/2015
          asum=asum+am
          zsum=zsum+zm
          totas = totas + am
          totzs = totzs + zm
!
          tet = results%progenyBnk(m)%theta*radiantodegree
          tt = results%progenyBnk(m)%kinEnergy*thousand
          ipoc = nint(results%progenyBnk(m)%prodMech)
          doit = ((gsmRxn%output%fisOnly.and.fisevent) .or. &
               & (.not.gsmRxn%output%fisOnly))
          if (doit) then
             if (am > amaxx .and. fisevent) totke = totke + tt
             if (ipar <= 6 .or. &
                  & (ipar > 9 .and. gsmRxn%output%nuclideCrossSection > 0)) then
!   Accumulate nuclide production information:
                da = abeg - am
                dz = zbeg - zm
                ida = nint(da) + 1
                ida = min(nint(abeg), ida)
                ida = min(ida, 350)
                idz = nint(dz) + 1
                idz = min(nint(zbeg) + 1, idz)
                idz = min(idz, 148)
!  04/15/2015
                if (am >= zm) then
                   dadz(ida,idz) = dadz(ida,idz) + one
                   dadz(351,idz) = dadz(351,idz) + one
                   dadz(350,idz) = dadz(350,idz) + tt
                   dadz(349,idz) = dadz(349,idz) + tt**2
                   dadz(ida,150) = dadz(ida,150) + one
                   dadz(ida,151) = dadz(ida,151) + tt
                   dadz(ida,149) = dadz(ida,149) + tt**2
                endif
                if (am >= zm) then
                   ama = results%progenyBnk(m)%restMass
                   eaz = results%progenyBnk(m)%kinEnergy + ama
                   paz = sqrt(abs(eaz**2 - ama**2))
                   temp = eaz
                   if ( temp < div0Lim .and. temp > -div0Lim ) then
                      temp = div0Lim
                      write(gsmObj%io%message,2000) "225"
                      call gsmObj%io%print(4, 3, gsmObj%io%message)
                   end if
!   kkg  04/15/2015
!               vz = paz/temp
                   vz=paz*results%progenyBnk(m)%cosTheta/temp
                   dadz4(1,ida,idz) = dadz4(1,ida,idz) + tt
                   dadz4(1,351,idz) = dadz4(1,351,idz) + tet
                   dadz4(1,350,idz) = dadz4(1,350,idz) + vz
                   dadz4(1,ida,151) = dadz4(1,ida,151) + tet
                   dadz4(1,ida,150) = dadz4(1,ida,150) + vz
                   dadz4(2,ida,idz) = dadz4(2,ida,idz) + tt**2
                   dadz4(2,351,idz) = dadz4(2,351,idz) + tet**2
                   dadz4(2,350,idz) = dadz4(2,350,idz) + vz**2
                   dadz4(2,ida,151) = dadz4(2,ida,151) + tet**2
                   dadz4(2,ida,150) = dadz4(2,ida,150) + vz**2
                   if (tet <= 90.0d0) then
                      dadz4(3,ida,idz) = dadz4(3,ida,idz) + one
                      dadz4(3,351,idz) = dadz4(3,351,idz) + one
                      dadz4(3,350,idz) = dadz4(3,350,idz) + tt
                      dadz4(3,349,idz) = dadz4(3,349,idz) + tt**2
                      dadz4(3,ida,150) = dadz4(3,ida,150) + one
                      dadz4(3,ida,151) = dadz4(3,ida,151) + tt
                      dadz4(3,ida,149) = dadz4(3,ida,149) + tt**2
                   else
                      dadz4(4,ida,idz) = dadz4(4,ida,idz) + one
                      dadz4(4,351,idz) = dadz4(4,351,idz) + one
                      dadz4(4,350,idz) = dadz4(4,350,idz) + tt
                      dadz4(4,349,idz) = dadz4(4,349,idz) + tt**2
                      dadz4(4,ida,150) = dadz4(4,ida,150) + one
                      dadz4(4,ida,151) = dadz4(4,ida,151) + tt
                      dadz4(4,ida,149) = dadz4(4,ida,149) + tt**2
                   endif
                endif
!  End changes 11/13/04.
                if (ipar > 9) go to 30
             endif
          endif
          iu = 0
          it = 0
          itet = 0
          if (gsmRxn%output%angularSpectra > 0 .or. &
               & gsmRxn%output%doubleDiffSpectra > 0) then
             if ( gsmRxn%output%deltaTheta < div0Lim .and. &
                  & gsmRxn%output%deltaTheta > -div0Lim ) then
                gsmRxn%output%deltaTheta = div0Lim
                write(gsmObj%io%message,2000) "268"
                call gsmObj%io%print(4, 3, gsmObj%io%message)
             end if
             ytet = tet/gsmRxn%output%deltaTheta
             itet = int(ytet) + 1
             itet = min (itet, 19)
          endif
          if ((gsmRxn%output%energySpectra > 0 .or. &
               & gsmRxn%output%doubleDiffSpectra > 0) .and. &
               & results%progenyBnk(m)%typeID < ten) then
             do jt = 1,ntt
                if (tt >= te1(jt) .and. tt < te2(jt)) go to 10
             end do
             write (16, 1000) m, (results%progenyBnk(m))
             jt = ntt
10           it = jt
          endif
          itt = it > 0 .and. it <= ntt
          if (gsmRxn%output%doubleDiffSpectra > 0) then
             iu = 0
             do ju = 1,ntet
                if (tet >= gsmRxn%output%angleBins(ju)%lowerBound .and. &
                     & tet < gsmRxn%output%angleBins(ju)%upperBound) then
                   iu = ju
                   go to 20
                endif
             end do
20           continue
          endif
!  Change to accumulate multiplicities ONLY on fission events;
!  AJS  4/16/04; KKG 11/13/04:
          if (doit) then
             if (ipar == 1) then
                nn = nn + 1
                if (ipoc < 100) then
                   ip = 1
                   nnc = nnc + 1
                elseif (ipoc == 100) then
                   ip = 2
                   nnp = nnp + 1
                elseif (ipoc > 100) then
                   ip = 3
                   if (.not.fisevent) then
                      if (ipoc >= 1000) nne = nne + 1
                   else
                      if (ipoc == 1000) nnprf = nnprf + 1
                      if (ipoc == 2000) nnpof = nnpof + 1
                   endif
                endif
             elseif (ipar == 2) then
                np = np + 1
                if (ipoc < 100) then
                   ip = 4
                elseif (ipoc == 100) then
                   ip = 5
                elseif (ipoc > 100) then
                   ip = 6
                endif
             elseif (ipar == 3) then
                nd = nd + 1
                if (ipoc == 200) outData%pcoal(1) = outData%pcoal(1) + one
                if (ipoc < 1000) then
                   if (ipoc == 100) then
                      ip = 7
!               else
!                 ip = 15
                   endif
                else
                   ip = 8
                endif
             elseif (ipar == 4) then
                nt = nt + 1
                if (ipoc == 200) outData%pcoal(2) = outData%pcoal(2) + one
                if (ipoc < 1000) then
                   if (ipoc == 100) then
                      ip = 9
!               else
!                 ip = 16
                   endif
                else
                   ip = 10
                endif
             elseif (ipar == 5) then
                nhe3 = nhe3 + 1
                if (ipoc == 200) outData%pcoal(3) = outData%pcoal(3) + one
                if (ipoc < 1000) then
                   if (ipoc == 100) then
                      ip = 11
!               else
!                 ip = 16
                   endif
                else
                   ip = 12
                endif
             elseif (ipar == 6) then
                nhe4 = nhe4 + 1
                if (ipoc == 200) outData%pcoal(4) = outData%pcoal(4) + one
                if (ipoc < 1000) then
                   if (ipoc == 100) then
                      ip = 13
!               else
!                 ip = 17
                   endif
                else
                   ip = 14
                endif
             elseif (ipar >= 7) then
                if (ipar == 7) then
                   npimin = npimin + 1
                elseif (ipar == 8) then
                   npizer = npizer + 1
                elseif (ipar == 9) then
                   npiplu = npiplu + 1
                endif
             endif
          endif
!  CHANGE TO ACCUMULATE SPECTRA ONLY ON FISSION EVENTS;
!  AJS  4/01/04:
          test = ipar <= 9
          if (gsmRxn%output%fisOnly) test = test.and.fisevent
          if (test) then
             spec(ipar,nt2) = spec(ipar,nt2) + tt
             spec(ipar,nt3) = spec(ipar,nt3) + one
             ang(ipar,20) = ang(ipar,20) + one
             if (itt) spec(ipar,it) = spec(ipar,it) + one
             if (itet > 0) ang(ipar,itet) = ang(ipar,itet) + one
             if (iu > 0) then
                d2spec(ipar,iu,nt2) = d2spec(ipar,iu,nt2) + tt
                d2spec(ipar,iu,nt3) = d2spec(ipar,iu,nt3) + one
                if (itt) d2spec(ipar,iu,it) = d2spec(ipar,iu,it) + one
             endif
          endif
!  CHANGE TO ACCUMULATE SPECTRA ONLY ON FISSION EVENTS;
!  AJS  4/01/04:
!  KKG  11/17/04
          test = ipar <= 6 .and. ip <= 14 .and. ipoc.ne.200
          if (gsmRxn%output%fisOnly) test = test.and.fisevent
          if (test) then
             spe(ip,nt2) = spe(ip,nt2) + tt
             spe(ip,nt3) = spe(ip,nt3) + one
             an(ip,20) = an(ip,20) + one
             if (itt) spe(ip,it) = spe(ip,it) + one
             if (itet > 0) an(ip,itet) = an(ip,itet) + one
             if (iu > 0) then
                d2spe(ip,iu,nt2) = d2spe(ip,iu,nt2) + tt
                d2spe(ip,iu,nt3) = d2spe(ip,iu,nt3) + one
                if (itt) d2spe(ip,iu,it) = d2spe(ip,iu,it) + one
             endif
          endif
          if (ipar <= 6 .and. ipoc == 2000) then
!   Evaporation from fission fragments:
             spef(ipar,nt2) = spef(ipar,nt2) + tt
             spef(ipar,nt3) = spef(ipar,nt3) + one
             angf(ipar,20)  = angf(ipar,20) + one
             if (itt) spef(ipar,it) = spef(ipar,it) + one
             if (itet > 0) angf(ipar,itet) = angf(ipar,itet) + one
             if (iu > 0) then
                d2spef(ipar,iu,nt2) = d2spef(ipar,iu,nt2) + tt
                d2spef(ipar,iu,nt3) = d2spef(ipar,iu,nt3) + one
                if (itt) d2spef(ipar,iu,it) = d2spef(ipar,iu,it) + one
             endif
          endif
          test = ipar <= 6 .and. ipoc >= 1000 .and. .not.fisevent
          if (gsmRxn%output%fisOnly) test = .false.
          if (test) then
!   Evaporation from nonfissioning spallation residues;
!   Also includes Fermi breakup:
             specsp(ipar,nt2) = specsp(ipar,nt2) + tt
             specsp(ipar,nt3) = specsp(ipar,nt3) + one
             angsp(ipar,20)   = angsp(ipar,20)   + one
             if (itt) specsp(ipar,it) = specsp(ipar,it) + one
             if (itet > 0) angsp(ipar,itet) = angsp(ipar,itet) + one
             if (iu > 0) then
                d2spesp(ipar,iu,nt2) = d2spesp(ipar,iu,nt2) + tt
                d2spesp(ipar,iu,nt3) = d2spesp(ipar,iu,nt3) + one
                if (itt) d2spesp(ipar,iu,it) = d2spesp(ipar,iu,it) + one
             endif
          endif
          if (ipar <= 6 .and. fisevent .and. ipoc == 1000) then
!   Prefission evaporation:
             specpf(ipar,nt2) = specpf(ipar,nt2) + tt
             specpf(ipar,nt3) = specpf(ipar,nt3) + one
             angpf(ipar,20)   = angpf(ipar,20)   + one
             if (itt) specpf(ipar,it) = specpf(ipar,it) + one
             if (itet > 0) angpf(ipar,itet) = angpf(ipar,itet) + one
             if (iu > 0) then
                d2spepf(ipar,iu,nt2) = d2spepf(ipar,iu,nt2) + tt
                d2spepf(ipar,iu,nt3) = d2spepf(ipar,iu,nt3) + one
                if (itt) d2spepf(ipar,iu,it) = d2spepf(ipar,iu,it) + one
             endif
          endif
!  Coalescence:
!  CHANGE TO ACCUMULATE SPECTRA ONLY ON FISSION EVENTS;
!  AJS  4/01/04:
          test = ipar >= 3 .and. ipar <= 6 .and. ipoc == 200
          if (gsmRxn%output%fisOnly) test = test.and.fisevent
          if (test) then
             iq = ipar - 2
             speco(iq,nt2) = speco(iq,nt2) + tt
             speco(iq,nt3) = speco(iq,nt3) + one
             angco(iq,20)  = angco(iq,20) + one
             if (itt) speco(iq,it) = speco(iq,it) + one
             if (itet > 0) angco(iq,itet) = angco(iq,itet) + one
             if (iu > 0) then
                d2speco(iq,iu,nt2) = d2speco(iq,iu,nt2) + tt
                d2speco(iq,iu,nt3) = d2speco(iq,iu,nt3) + one
                if (itt) d2speco(iq,iu,it) = d2speco(iq,iu,it) + one
             endif
          endif
30        continue
       end do
!  End of m DO loop ^
40     continue
!  Check that A and Z are conserved following suggestion of
!  G. McKinney and K. Gudima; print diagnostic on "info" file
!  and to terminal.  AJS (02/07/06):
       if (abs(totas - targ%numBaryons - dble(proj%numBaryons)) > 0.1d0 .or. &
            & abs(totzs - targ%numProtons - dble(proj%numProtons)) > 0.1d0) then
!  kkg  04/15/2015
          if ( gsmRxn%incFlag == sDCMFlagged ) then
             write ( gsmObj%io%message, 1100) gsmRxn%eventAttempt, totzs, totas
             call gsmObj%io%print(2, 2, gsmObj%io%message)
             write (16, 1100) gsmRxn%eventAttempt, totzs, totas
!           stop
          endif
       endif
!
       if (gsmRxn%output%channelCrossSection <= 0) go to 50
       if (nt > 0) chan(214) = chan(214) + dble(nt)
       if (nhe3 > 0) chan(215) = chan(215) + dble(nhe3)
       if (npimin > 0) chan(216) = chan(216) + dble(npimin)
       if (npizer > 0) chan(217) = chan(217) + dble(npizer)
       if (npiplu > 0) chan(218) = chan(218) + dble(npiplu)
       ic = 0
       if (nn == 1) then
!    This is (1n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 1
!    This is (np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 26
!    This is (n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 45
!    This is (n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 63
!    This is (n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 80
!    This is (n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 96
!    This is (nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 112
!    This is (npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 132
!    This is (n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 147
!    This is (nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 162
!    This is (npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 179
!    This is (n2He4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 2) ic = 176
       elseif (nn == 2) then
!    This is (2n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 2
!    This is (2np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 27
!    This is (2n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 46
!    This is (2n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 64
!    This is (2n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 81
!    This is (2n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 97
!    This is (2nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 113
!    This is (2npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 133
!    This is (2n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 148
!    This is (2nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 163
!    This is (2npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 180
!    This is (2n2He4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 2) ic = 177
       elseif (nn == 3) then
!    This is (3n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 3
!    This is (3np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 28
!    This is (3n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 47
!    This is (3n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 65
!    This is (3n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 82
!    This is (3n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 98
!    This is (3nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 114
!    This is (3npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 134
!    This is (3n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 149
!    This is (3nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 164
!    This is (3npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 181
       elseif (nn == 4) then
!    This is (4n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 4
!    This is (4np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 29
!    This is (4n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 48
!    This is (4n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 66
!    This is (4n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 83
!    This is (4n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 99
!    This is (4nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 115
!    This is (4npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 135
!    This is (4n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 150
!    This is (4nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 165
!    This is (4npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 182
       elseif (nn == 5) then
!    This is (5n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 5
!    This is (5np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 30
!    This is (5n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 49
!    This is (5n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 67
!    This is (5n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 84
!    This is (5n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 100
!    This is (5nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 116
!    This is (5npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 136
!    This is (5n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 151
!    This is (5nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 166
!    This is (5npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 183
       elseif (nn == 6) then
!    This is (6n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 6
!    This is (6np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 31
!    This is (6n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 50
!    This is (6n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 68
!    This is (6n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 85
!    This is (6n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 101
!    This is (6nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 117
!    This is (6npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 137
!    This is (6n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 152
!    This is (6nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 167
!    This is (6npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 184
       elseif (nn == 7) then
!    This is (7n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 7
!    This is (7np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 32
!    This is (7n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 51
!    This is (7n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 69
!    This is (7n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 86
!    This is (7n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 102
!    This is (7nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 118
!    This is (7npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 138
!    This is (7n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 153
!    This is (7nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 168
!    This is (7npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 185
       elseif (nn == 8) then
!    This is (8n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 8
!    This is (8np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 33
!    This is (8n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 52
!    This is (8n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 70
!    This is (8n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 87
!    This is (8n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 103
!    This is (8nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 119
!    This is (8npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 139
!    This is (8n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 154
!    This is (8nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 169
!    This is (8npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 186
       elseif (nn == 9) then
!    This is (9n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 9
!    This is (9np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 34
!    This is (9n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 53
!    This is (9n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 71
!    This is (9n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 88
!    This is (9n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 104
!    This is (9nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 120
!    This is (9npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 140
!    This is (9n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 155
!    This is (9nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 170
!    This is (9npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 187
       elseif (nn == 10) then
!    This is (10n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 10
!    This is (10np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 35
!    This is (10n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 54
!    This is (10n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 72
!    This is (10n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 89
!    This is (10n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 105
!    This is (10nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 121
!    This is (10npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 141
!    This is (10n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 156
!    This is (10nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 171
!    This is (10npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 188
       elseif (nn == 11) then
!    This is (11n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 11
!    This is (11np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 36
!    This is (11n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 55
!    This is (11n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 73
!    This is (11n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 90
!    This is (11n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 106
!    This is (11nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 122
!    This is (11npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 142
!    This is (11n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 157
!    This is (11nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 172
!    This is (11npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 189
       elseif (nn == 12) then
!    This is (12n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 12
!    This is (12np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 37
!    This is (12n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 56
!    This is (12n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 74
!    This is (12n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 91
!    This is (12n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 107
!    This is (12nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 123
!    This is (12npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 143
!    This is (12n2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 158
!    This is (12nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 173
!    This is (12npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 190
       elseif (nn == 13) then
!    This is (13n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 13
!    This is (13np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 38
!    This is (13n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 57
!    This is (13n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 75
!    This is (13n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 92
!    This is (13n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 108
!    This is (13nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 124
!    This is (13npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 144
!    This is (13nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 174
!    This is (13npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 191
       elseif (nn == 14) then
!    This is (14n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 14
!    This is (14np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 39
!    This is (14n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 58
!    This is (14n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 76
!    This is (14n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 93
!    This is (14n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 109
!    This is (14nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 125
!    This is (14npd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 145
!    This is (14nHe4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 175
!    This is (14npHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 192
       elseif (nn == 15) then
!    This is (15n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 15
!    This is (15np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 40
!    This is (15n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 59
!    This is (15n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 77
!    This is (15n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 94
!    This is (15n5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 110
!    This is (15nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 126
       elseif (nn == 16) then
!    This is (16n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 16
!    This is (16np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 41
!    This is (16n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 60
!    This is (16n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 78
!    This is (16n4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 95
!    This is (16nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 127
       elseif (nn == 17) then
!    This is (17n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 17
!    This is (17np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 42
!    This is (17n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 61
!    This is (17n3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 79
!    This is (17nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 128
       elseif (nn == 18) then
!    This is (18n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 18
!    This is (18np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 43
!    This is (18n2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 62
!    This is (18nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 129
       elseif (nn == 19) then
!    This is (19n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 19
!    This is (19np) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 44
!    This is (19nd) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 130
       elseif (nn == 20) then
!    This is (20n) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 20
       elseif (nn == 0) then
!    This is (1p) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 21
!    This is (2p) channel:
          if (np == 2 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 22
!    This is (3p) channel:
          if (np == 3 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 23
!    This is (4p) channel:
          if (np == 4 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 24
!    This is (5p) channel:
          if (np == 5 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 25
!    This is (d) channel:
          if (np == 0 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 111
!    This is (pd) channel:
          if (np == 1 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 131
!    This is (2pd) channel:
          if (np == 2 .and. nd == 1 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 0) ic = 146
!    This is (He4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 159
!    This is (2He4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 2) ic = 160
!    This is (3He4) channel:
          if (np == 0 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 3) ic = 161
!    This is (pHe4) channel:
          if (np == 1 .and. nd == 0 .and. nt == 0 .and. nhe3 == 0 &
               & .and. nhe4 == 1) ic = 178
       endif
       icp = 0
       if (nt == 0 .and. nhe3 == 0) then
          if (nd == 0 .and. nhe4 == 0) then
             if (nn == 0) then
!    This is (Xp) channel:
                if (np > 0) icp = 194
             elseif (nn > 0) then
!    This is (Xn) channel:
                if (np == 0) icp = 193
!    This is (Xn+Yp) channel:
                if (np > 0) icp = 195
             endif
          elseif (np == 0 .and. nn > 0) then
!    This is (Xn+d) channel:
             if (nd == 1 .and. nhe4 == 0) icp = 196
!    This is (Xn+2d) channel:
             if (nd == 2 .and. nhe4 == 0) icp = 197
!    This is (Xn+He4) channel:
             if (nd == 0 .and. nhe4 == 1) icp = 208
!    This is (Xn+2He4) channel:
             if (nd == 0 .and. nhe4 == 2) icp = 209
          elseif (nhe4 == 0 .and. nn > 0) then
!    This is (Xn+p+d) channel:
             if (nd == 1 .and. np == 1) icp = 198
!    This is (Xn+2p+d) channel:
             if (nd == 1 .and. np == 2) icp = 199
!    This is (Xn+p+2d) channel:
             if (nd == 2 .and. np == 1) icp = 200
!    This is (Xn+2p+2d) channel:
             if (nd == 2 .and. np == 2) icp = 201
          elseif (nd == 0 .and. nn > 0) then
!    This is (Xn+p+He4) channel:
             if (np == 1 .and. nhe4 == 1) icp = 210
!    This is (Xn+2p+He4) channel:
             if (np == 2 .and. nhe4 == 1) icp = 211
!    This is (Xn+3p+He4) channel:
             if (np == 3 .and. nhe4 == 1) icp = 212
          elseif (nd == 1 .and. nn > 0) then
!    This is (Xn+d+He4) channel:
             if (nhe4 == 1 .and. np == 0) icp = 213
          endif
       elseif (nt > 0 .and. nhe3 == 0) then
!    This is (Xn+t) channel:
          if (np == 0 .and. nn > 0 .and. nd == 0 .and. nhe4 == 0 &
               & .and. nt == 1) icp = 202
!    This is (Xn+2t) channel:
          if (np == 0 .and. nn > 0 .and. nd == 0 .and. nhe4 == 0 &
               & .and. nt == 2) icp = 203
!    This is (Xn+p+t) channel:
          if (np == 1 .and. nn > 0 .and. nd == 0 .and. nhe4 == 0 &
               & .and. nt == 1) icp = 204
!    This is (Xn+2p+t) channel:
          if (np == 2 .and. nn > 0 .and. nd == 0 .and. nhe4 == 0 &
               & .and. nt == 1) icp = 205
!    This is (Xn+d+t) channel:
          if (np == 0 .and. nn > 0 .and. nd == 1 .and. nhe4 == 0 &
               & .and. nt == 1) icp = 206
!    This is (Xn+p+d+t) channel:
          if (np == 1 .and. nn > 0 .and. nd == 1 .and. nhe4 == 0 &
               & .and. nt == 1) icp = 207
       endif
       if (ic.ne.0) chan(ic) = chan(ic) + one
       if (icp > 0) chan(icp) = chan(icp) + one
    endif
50  continue
    if (gsmRxn%output%nuclideCrossSection == 3) then
!   KKG  12/02/04   Accumulate residual nuclei information:
       call gsmObj%resdist (results%info, zbeg, abeg, fisevent)
!      Accumulate distribution of fission-fragment opening angles:
       call gsmObj%opandist (results%info, fisevent, nn)
!      Accumulate distribution of neutron multiplicities:
       call gsmObj%disnmul (fisevent, nn, nnc, nnp, nne, nnprf, nnpof)
    endif

    ! Tally for PISA now
    if(gsmRxn%output%printPISA) call gsmObj%pisaSpectra(gsmRxn%output, results)

    !$OMP end critical

    return

! ======================================================================
1000 format (1x,'Tallying (vlobd): particle # ',i4,' has an energy outside ', &
         & 'the range for the spectrum; progeny details: ='/, &
         & 2f5.1, /, & ! A, Z
         & 2f10.6, /, & ! E_k, E_0
         & 2f6.2, 2f6.4, /, & ! phi, theta, sin(Theta), cos(Theta)
         & 2f8.4 )
1100 format (1x,'Tallying (vlobd): for event # ',i8,', nonconservation of ', &
         & 'z = ',f5.1,' and/or a = ',f5.1,' !')
2000 format("Divide by zero error in 'vlobd.f90', line(s) ", A)
! ======================================================================
  end subroutine vlobd
