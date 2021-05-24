
  subroutine typeout (gsmObj, proj, targ, output, outData)

! ======================================================================
!
!     Generates printout of results on unit 31 (cem2kgemfit.res).
!     Simplified by putting various loops into a separate subroutine
!     PRTDIST, AJS, 10/24/03. Also removed the calculation and printing
!     of the particle multiplicity table to PRTMULT, 12/15/03.
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Modified by A. J. Sierk, April-May, 1996.
!   Modified by AJS, Dec 1997; March, 1998.
!   Modified by AJS, March, 1999.
!   Modified by SGM at 06/16/2000
!   Modified by SGM at 05/08/2001
!   Modified by SGM on 12/04/01
!   "Last change": 13-AUG-2003 by NVMokhov
!   Modfied by AJS; September-October, 2003.
!   Modfied by AJS; December, 2003.
!   Modfied by KKG; September-November, 2004.
!   Edited by AJS, January, 2005.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!   Edited by CMJ, XCP-3, June 2017 (LAQGSM expansion)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64, int64
    use OutputDataMod, only: OutputData
    use gsm_params, only: zro, one, two, fiv, thousand, degreeToRadian, &
         & twpi
    use standarddcmdata, only: photonEG

    implicit none
    class(GSM),           intent(inout) :: gsmObj
    class(GSMProjectile), intent(in   ) :: proj
    class(GSMTarget),     intent(in   ) :: targ
    class(GSMOutput),     intent(in   ) :: output
    class(OutputData),    intent(inout) :: outData

    integer(int32)    :: i, ia, iaa, iamax, ic, ii, iii, ij, ik, in, &
         & ip, ipar, ipt, ipz, iz, iz1, iz2, iz3, izmax, j, k, lpr, &
         & lprm, m, n
    integer(int64)    :: npt
    real(real64)      :: aav, abeg, aeqav, apqav, arav, artav, atk, atks, &
         & bfav, datk, datks, def, defmu, defp, defpmu, dfis, dsigf, dy1, &
         & dy2, dy3, dys1, dys2, dys3, dyy, dyys, eeav, ef, efmu, efp, &
         & efpmu, elav, eleav, elpav, elrav, elrtav, epav, eqqv, erkav, &
         & estarav, faca, fact, fe, fision, fl, fn, fneq, fnfis, fnn, fnpreq, &
         & fnres, fnret, ppt, s1, sfu1, siga, sigaeq, sigapq, sigar, &
         & sigart, sigbf, sigee, sigel, sigep, sigeq, sigerk, sigest, sigfis, &
         & sigfw, &
         & sigg, sigin, siginn, sigl, sigla, sigle, siglp, siglr, siglrt, &
         & sigz, sigzeq, sigzpq, sigzr, sigzrt, ss, t0mev, temp, test, &
         & tetmax, tetmin, tkfbar, y1, y2, y3, ys1, ys2, ys3, yy, yys, zav, &
         & zbeg, zeqav, zpqav, zrav, zrtav
    logical           :: printit, prtz1, prtz2, prtz3
    character(LEN=10) :: dummc = ""
    character(LEN=24) :: ddate = ""

    integer(int32),    dimension(151) :: izr = 0_int32
    integer(int32),    dimension(351) :: iar = 0_int32
    real(real64),      dimension( 20) ::  om = zro, dom = zro, &
         & const = zro, doo = zro
    real(real64),      dimension(151) ::  z = zro, dasum = zro
    real(real64),      dimension(351) :: ar = zro
    logical,           dimension(151) :: prtz = .FALSE.
    character(LEN=10), dimension(218) :: gch = ""

! ======================================================================

    character(LEN=10), parameter, dimension(9) :: par2 = &
         & [' neutrons ', ' protons  ', 'deuterons ', ' tritons  ', &
         &  ' Helium-3 ', '  alphas  ', 'neg. pions', 'neut pions', &
         &  'pos. pions' ]

    character(LEN=10), parameter, dimension(218) :: ch = [ &
         ! 1-72
         & '(N,n)     ', '(N,2n)    ', '(N,3n)    ', '(N,4n)    ', &
         & '(N,5n)    ', '(N,6n)    ', '(N,7n)    ', '(N,8n)    ', &
         & '(N,9n)    ', '(N,10n)   ', '(N,11n)   ', '(N,12n)   ', &
         & '(N,13n)   ', '(N,14n)   ', '(N,15n)   ', '(N,16n)   ', &
         & '(N,17n)   ', '(N,18n)   ', '(N,19n)   ', '(N,20n)   ', &
         & '(N,p)     ', '(N,2p)    ', '(N,3p)    ', '(N,4p)    ', &
         & '(N,5p)    ', '(N,n+p)   ', '(N,2n+p)  ', '(N,3n+p)  ', &
         & '(N,4n+p)  ', '(N,5n+p)  ', '(N,6n+p)  ', '(N,7n+p)  ', &
         & '(N,8n+p)  ', '(N,9n+p)  ', '(N,10n+p) ', '(N,11n+p) ', &
         & '(N,12n+p) ', '(N,13n+p) ', '(N,14n+p) ', '(N,15n+p) ', &
         & '(N,16n+p) ', '(N,17n+p) ', '(N,18n+p) ', '(N,19n+p) ', &
         & '(N,n+2p)  ', '(N,2n+2p) ', '(N,3n+2p) ', '(N,4n+2p) ', &
         & '(N,5n+2p) ', '(N,6n+2p) ', '(N,7n+2p) ', '(N,8n+2p) ', &
         & '(N,9n+2p) ', '(N,10n+2p)', '(N,11n+2p)', '(N,12n+2p)', &
         & '(N,13n+2p)', '(N,14n+2p)', '(N,15n+2p)', '(N,16n+2p)', &
         & '(N,17n+2p)', '(N,18n+2p)', '(N,n+3p)  ', '(N,2n+3p) ', &
         & '(N,3n+3p) ', '(N,4n+3p) ', '(N,5n+3p) ', '(N,6n+3p) ', &
         & '(N,7n+3p) ', '(N,8n+3p) ', '(N,9n+3p) ', '(N,10n+3p)', & ! 1-72
         ! 73-144
         & '(N,11n+3p)', '(N,12n+3p)', '(N,13n+3p)', '(N,14n+3p)', &
         & '(N,15n+3p)', '(N,16n+3p)', '(N,17n+3p)', '(N,n+4p)  ', &
         & '(N,2n+4p) ', '(N,3n+4p) ', '(N,4n+4p) ', '(N,5n+4p) ', &
         & '(N,6n+4p) ', '(N,7n+4p) ', '(N,8n+4p) ', '(N,9n+4p) ', &
         & '(N,10n+4p)', '(N,11n+4p)', '(N,12n+4p)', '(N,13n+4p)', &
         & '(N,14n+4p)', '(N,15n+4p)', '(N,16n+4p)', '(N,n+5p)  ', &
         & '(N,2n+5p) ', '(N,3n+5p) ', '(N,4n+5p) ', '(N,5n+5p) ', &
         & '(N,6n+5p) ', '(N,7n+5p) ', '(N,8n+5p) ', '(N,9n+5p) ', &
         & '(N,10n+5p)', '(N,11n+5p)', '(N,12n+5p)', '(N,13n+5p)', &
         & '(N,14n+5p)', '(N,15n+5p)', '(N,d)     ', '(N,n+d)   ', &
         & '(N,2n+d)  ', '(N,3n+d)  ', '(N,4n+d)  ', '(N,5n+d)  ', &
         & '(N,6n+d)  ', '(N,7n+d)  ', '(N,8n+d)  ', '(N,9n+d)  ', &
         & '(N,10n+d) ', '(N,11n+d) ', '(N,12n+d) ', '(N,13n+d) ', &
         & '(N,14n+d) ', '(N,15n+d) ', '(N,16n+d) ', '(N,17n+d) ', &
         & '(N,18n+d) ', '(N,19n+d) ', '(N,p+d)   ', '(N,n+p+d) ', &
         & '(N,2n+p+d)', '(N,3n+p+d)', '(N,4n+p+d)', '(N,5n+p+d)', &
         & '(N,6n+p+d)', '(N,7n+p+d)', '(N,8n+p+d)', '(N,9n+p+d)', &
         & '(N,10n+pd)', '(N,11n+pd)', '(N,12n+pd)', '(N,13n+pd)', & ! 73-144
         ! 145-218
         & '(N,14n+pd)', '(N,2p+d)  ', '(N,n+2p+d)', '(N,2n+2pd)', &
         & '(N,3n+2pd)', '(N,4n+2pd)', '(N,5n+2pd)', '(N,6n+2pd)', &
         & '(N,7n+2pd)', '(N,8n+2pd)', '(N,9n+2pd)', '(N,10n2pd)', &
         & '(N,11n2pd)', '(N,12n2pd)', '(N,He4)   ', '(N,2He4)  ', &
         & '(N,3He4)  ', '(N,n+He4) ', '(N,2n+He4)', '(N,3n+He4)', &
         & '(N,4n+He4)', '(N,5n+He4)', '(N,6n+He4)', '(N,7n+He4)', &
         & '(N,8n+He4)', '(N,9n+He4)', '(N,10nHe4)', '(N,11nHe4)', &
         & '(N,12nHe4)', '(N,13nHe4)', '(N,14nHe4)', '(N,n+2He4)', &
         & '(N,2n2He4)', '(N,p+He4) ', '(N,np+He4)', '(N,2npHe4)', &
         & '(N,3npHe4)', '(N,4npHe4)', '(N,5npHe4)', '(N,6npHe4)', &
         & '(N,7npHe4)', '(N,8npHe4)', '(N,9npHe4)', '(N,10npHe)', &
         & '(N,11npHe)', '(N,12npHe)', '(N,13npHe)', '(N,14npHe)', &
         & '(N,xn)    ', '(N,xp)    ', '(N,xn+yp) ', '(N,xn+d)  ', &
         & '(N,xn+2d) ', '(N,xn+p+d)', '(N,xn+2pd)', '(N,xn+p2d)', &
         & '(N,xn2p2d)', '(N,xn+t)  ', '(N,xn+2t) ', '(N,xn+p+t)', &
         & '(N,xn+2pt)', '(N,xn+d+t)', '(N,xn+pdt)', '(N,xn+He4)', &
         & '(N,xn2He4)', '(N,xnpHe4)', '(N,xn2pHe)', '(N,xn3pHe)', &
         & '(N,xndHe4)', '(N,t+x)   ', '(N,He3+x) ', '(N,pi- +x)', &
         & '(N,pi0 +x)', '(N,pi+ +x)'] ! 145-218

! ======================================================================

    integer(int32) :: ntt, ntet, nt2, nt3, nti
    real(real64)   :: te1, te2, dtt, se, dtdo, d2spec, d2spe
    common /d2sdto/ te1(200), te2(200), dtt(200), &
         & se(200), dtdo(200,10), d2spec(9,10,200), &
         & d2spe(14,10,200), ntt, ntet, nt2, nt3, nti(4)
    integer(int64) :: neq
    real(real64)   :: eestt, eestsq, aeqtot, aeqsq, zeqtot, zeqsq, &
         & aemin, aemax, eletot, elesq, elemin, elemax, &
         & eestrmn, eestrmx, zemin, zemax
    common /eqidat/  eestt, eestsq, aeqtot, aeqsq, zeqtot, zeqsq, &
         & aemin, aemax, eletot, elesq, elemin, elemax, &
         & eestrmn, eestrmx, zemin, zemax, neq
    real(real64)   :: estart, estarsq, estrmn, estrmx, atot, atsq, &
         & ammin, ammax, eltot, elsq, elmmin, elmmax, ztot, &
         & ztsq, zmmin, zmmax, bftot, bfsq, bfmin, bfmax
    common /fisda2/  estart, estarsq, estrmn, estrmx, atot, atsq, &
         & ammin, ammax, eltot, elsq, elmmin, elmmax, ztot, &
         & ztsq, zmmin, zmmax, bftot, bfsq, bfmin, bfmax
    integer(int32) :: no
    integer(int64) :: npreq
    real(real64)   :: epstt, epstsq, apqtot, apqsq, zpqtot, zpqsq, &
         & apmin, apmax, elptot, elpsq, elpmin, elpmax, &
         & epstrmn, epstrmx, zpmin, zpmax
    common /predat/  epstt, epstsq, apqtot, apqsq, zpqtot, zpqsq, &
         & apmin, apmax, elptot, elpsq, elpmin, elpmax, &
         & epstrmn, epstrmx, zpmin, zpmax, npreq
    integer(int64) :: nres
    real(real64)   :: erkt, erksq, artot, arsq, zrtot, zrsq, &
         & armin, armax, elrtot, elrsq, elrmin, elrmax, &
         & erkmn, erkmx, zrmin, zrmax
    common /resdat/  erkt, erksq, artot, arsq, zrtot, zrsq, &
         & armin, armax, elrtot, elrsq, elrmin, elrmax, &
         & erkmn, erkmx, zrmin, zrmax, nres
    real(real64)   :: spec, ang, chan, dadz
    common /result/  spec(9,200), ang(9,20), chan(218), dadz(351,151)
    integer(int64) :: nret
    real(real64)   :: arttot, artsq, artmin, artmax, zrttot, zrtsq, &
         & zrtmin, zrtmax, elrttot, elrtsq, elrtmin, &
         & elrtmax
    common /retdat/  arttot, artsq, artmin, artmax, zrttot, zrtsq, &
         & zrtmin, zrtmax, elrttot, elrtsq, elrtmin, &
         & elrtmax, nret
    integer(int32) :: istp
    common /stopr/   istp
    real(real64)   :: totke
    common /tkftot/  totke

! ======================================================================

    ! Require only one thread to enter output printing at a time
    !$OMP critical

    fl = dble(outData%limcc)
    fn = dble(outData%ncas)
    fe = dble(outData%intel)
    fneq = dble(neq)
    fnpreq = dble(npreq)
    fnres = dble(nres)
    fnret = dble(nret)
    fnfis = dble(outData%nfis)
    temp = fl
    if ( temp < div0Lim .and. temp > -div0Lim ) then
       temp = div0Lim
       write(gsmObj%io%message, 5000) "229"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    sigin = outData%sigom * fn / temp
    t0mev = proj%kinEnergy * thousand
!  KKG 09/17/04
    sigla = sigin   ! use the calculated MC cross section by default
    if (proj%numBaryons >= 1) then
       ! 'iii' id's a n, p, d, t, he3, or he4 particle (1-6, respectively)
       iii = proj%numProtons + proj%numBaryons
       if ( iii <= 6 .and. iii > 0 ) then
          sigla = gsmObj%sinvla (iii, targ%numBaryons, targ%numProtons, t0mev)
       endif
    endif
! Photon interaction:
!   KKG  10/13/04
    if ( fn < div0Lim .and. fn > -div0Lim ) then
       fn = div0Lim
       write(gsmObj%io%message, 5000) "217,218, 232, 234"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    if (proj%particleFlag == photonProjFlag .or. &
        & proj%particleFlag == bremsProjFlag) then
       if (proj%particleFlag == bremsProjFlag) then
!         eqqv = proj%brems%tEqv()/(fn*(proj%brems%tMax() - proj%brems%tMin()))
!          if ( proj%brems%tMax() < div0Lim .and. proj%brems%tMax() > -div0Lim ) then
!             call proj%brems%setTMax(div0Lim)
!             write(gsmObj%io%message, 5000) "217"
!             call gsmObj%io%print(4, 3, gsmObj%io%message)
!          end if
          eqqv = proj%brems%tEqv()/(fn*proj%brems%tMax())
          sigla = proj%brems%sxabs()/fn
          if ( eqqv < div0Lim .and. eqqv > -div0Lim ) then
             eqqv = div0Lim
             write(gsmObj%io%message, 5000) "220"
             call gsmObj%io%print(4, 3, gsmObj%io%message)
          end if
          sigeq = sigla/eqqv
          write (31, 999) eqqv, sigeq, sigla
          sigla = sigeq
       else
          sigla = photonEG%photocrosssection (t0mev, targ%numBaryons)
       endif
       test = 0.01d0
    else
       test = one
    endif
    write (31, 1000) sigla, sigin
    sigin = sigla
    outData%sigom = sigla*fl/fn
    sigel = outData%sigom*fe/temp
    siginn = sigin/fn
    dfis = zro
    dsigf = zro
    estarav = zro
    sigest = zro
    zav = zro
    aav = zro
    elav = zro
    sigz = zro
    siga = zro
    sigl = zro
    if (outData%sfu.ne.zro) then
       sfu1 = outData%sfu/fn
       sigfw = sfu1*sigin
       if (outData%nfis.ne.0) then
          fision = fnfis/fn
          dfis = sqrt(abs(fnfis))/fn
          sigfis = fision*sigin
          tkfbar = totke/fnfis
          dsigf = dfis*sigin
          estarav = estart/fnfis
          sigest = sqrt(abs(estarsq/fnfis - estarav**2))
          bfav = bftot/fnfis
          sigbf = sqrt(abs(bfsq/fnfis - bfav**2))
          zav = ztot/fnfis
          sigz = sqrt(abs(ztsq/fnfis - zav**2))
          aav = atot/fnfis
          siga = sqrt(abs(atsq/fnfis - aav**2))
          elav = eltot/fnfis
          sigl = sqrt(abs(elsq/fnfis - elav**2))
!         if (iz0 > 0) then
!           a0av = azro/dble(iz0)
!           el0av = elzro/dble(iz0)
!         endif
       endif
    endif
    if ( fnpreq < div0Lim .and. fnpreq > -div0Lim ) then
       fnpreq = div0Lim
       write(gsmObj%io%message, 5000) "272-279"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    epav = epstt/fnpreq
    zpqav = zpqtot/fnpreq
    apqav = apqtot/fnpreq
    elpav = elptot/fnpreq
    sigep = sqrt(abs(epstsq/fnpreq - epav**2))
    sigzpq = sqrt(abs(zpqsq/fnpreq - zpqav**2))
    sigapq = sqrt(abs(apqsq/fnpreq - apqav**2))
    siglp = sqrt(abs(elpsq/fnpreq - elpav**2))
    if (fneq > zro) then
       eeav = eestt/fneq
       zeqav = zeqtot/fneq
       aeqav = aeqtot/fneq
       eleav = eletot/fneq
       sigee = sqrt(abs(eestsq/fneq - eeav**2))
       sigzeq = sqrt(abs(zeqsq/fneq - zeqav**2))
       sigaeq = sqrt(abs(aeqsq/fneq - aeqav**2))
       sigle = sqrt(abs(elesq/fneq - eleav**2))
    endif
    if (fnres > zro) then
       erkav = erkt/fnres
       zrav = zrtot/fnres
       arav = artot/fnres
       elrav = elrtot/fnres
       sigerk = sqrt(abs(erksq/fnres - erkav**2))
       sigzr = sqrt(abs(zrsq/fnres - zrav**2))
       sigar = sqrt(abs(arsq/fnres - arav**2))
       siglr = sqrt(abs(elrsq/fnres - elrav**2))
    endif
    if (fnret > zro) then
       zrtav = zrttot/fnret
       artav = arttot/fnret
       elrtav = elrttot/fnret
       sigzrt = sqrt(abs(zrtsq/fnret - zrtav**2))
       sigart = sqrt(abs(artsq/fnret - artav**2))
       siglrt = sqrt(abs(elrtsq/fnret - elrtav**2))
    endif

    call fdate (ddate)
    write ( *, 1900) ddate
    write (31, 1900) ddate
    if (proj%particleFlag == nucleusProjFlag .or. &
        & proj%particleFlag == pionProjFlag) then
       if ( t0mev < 1000) then
          write (31, 1700) t0mev, proj%numProtons, proj%numBaryons, &
              & targ%numProtons, targ%numBaryons, outData%ncas, outData%intel, &
              & sigin, sigel
       else
          write (31, 1725) t0mev/1000, proj%numProtons, proj%numBaryons, &
              & targ%numProtons, targ%numBaryons, outData%ncas, &
              & outData%intel, sigin, sigel
       endif
    else
       if ( t0mev < 1000) then
          write (31, 1750) t0mev, proj%numProtons, proj%numBaryons, &
              & targ%numProtons, targ%numBaryons, outData%ncas, outData%intel, &
              & sigin, sigel
       else
          write (31, 1775) t0mev/1000, proj%numProtons, proj%numBaryons, &
              & targ%numProtons, targ%numBaryons, outData%ncas, &
              & outData%intel, sigin, sigel
       endif
    endif
      ! Displays reaction to user at end of calculation
    if ( t0mev < 1000 ) then
       ! Displays t0mev in units of MeV
       write ( *, 1800) t0mev, proj%numProtons, proj%numBaryons, &
           &  targ%numProtons, targ%numBaryons
    elseif ( t0mev < 1000000 ) then
       ! Displays t0mev in units of GeV
       write ( *, 1825) (t0mev/1000), proj%numProtons, proj%numBaryons, &
           & targ%numProtons, targ%numBaryons
    elseif ( t0mev < 1000000000 ) then
         ! Displays t0mev in units of TeV
       write ( *, 1850) (t0mev/1000000), proj%numProtons, proj%numBaryons, &
           & targ%numProtons, targ%numBaryons
    else
         ! Displays in units of MeV in scientific notation
       write ( *, 1875) t0mev, proj%numProtons, proj%numBaryons, &
           & targ%numProtons, targ%numBaryons
    endif

    if (npreq > 0) &
         & write (31, 2200) npreq, epav, sigep, epstrmn, epstrmx, &
         & zpqav, sigzpq, zpmin, zpmax, apqav, sigapq, &
         & apmin, apmax, elpav, siglp, elpmin, elpmax
    if (istp > zro) write (31, 2300) istp
    if (outData%ifermi > zro) write (31, 4800) outData%ifermi
    if (nret > 0) &
         & write (31, 2400) nret, zrtav, sigzrt, zrtmin, zrtmax, artav, &
         & sigart, artmin, artmax, elrtav, siglrt, &
         & elrtmin, elrtmax
    if (gsmObj%options%preeqOpts%excludePreeq == 0 .and. neq > 0) then
       write (31, 3100) neq, eeav, sigee, eestrmn, eestrmx, &
            & zeqav, sigzeq, zemin, zemax, aeqav, sigaeq, &
            & aemin, aemax, eleav, sigle, elemin, elemax
    endif
    if (nres > 0) &
         & write (31, 2600) nres, erkav, sigerk, erkmn, erkmx, &
         & zrav, sigzr, zrmin, zrmax, arav, sigar, &
         & armin, armax, elrav, siglr, elrmin, elrmax
    if (outData%sfu.ne.zro .and. outData%nfis == 0) write (31, 2100) sfu1, sigfw
    if (outData%nfis > 0) then
       write (31, 2500) outData%nfis, estarav, sigest, estrmn, estrmx, &
            & zav, sigz, zmmin, zmmax, aav, siga, &
            & ammin, ammax, elav, sigl, elmmin, elmmax, &
            & bfav, sigbf, bfmin, bfmax
       write (31, 4900) tkfbar
       write (31, 2000) fision, dfis, sigfis, dsigf
       if (outData%sfu.ne.zro) write (31, 2100) sfu1, sigfw
!       if (infis == 1) write (31, 2700) infis
!       if (infis > 1) write (31, 2800) infis
!       if (iz0 == 1) then
!         write (31, 3000) iz0, a0av, el0av
!       elseif (iz0 > 1) then
!         write (31, 2900) iz0, a0av, el0av
!       endif
    endif
!  Print particle-multiplicity table:
    if (output%multiplicities > 0 .or. &
         & output%energySpectra > 0 .or. &
         & output%doubleDiffSpectra > 0. .or. &
         & output%nuclideCrossSection > 0 .or. &
         & output%angularSpectra > 0) then
       sigg = sigin
       fnn = fn
       if (output%fisOnly) then
!  Multiplicity for fissioning nuclei only:
          sigg = sigfis
          fnn = fnfis
       endif
       if (output%multiplicities > 0 .or. &
            & output%energySpectra > 0 .or. &
            & output%doubleDiffSpectra > 0. &
            & .or. output%angularSpectra > 0) then
          do ipar = 1,9
             call gsmObj%prtmult (output, outData, ipar, fnn, sigg, nt2, nt3)
          end do
       endif
    endif
!
!   Print out cross sections for various particle-emission channels:
!
    if (output%channelCrossSection > 0) then
       fact = siginn
       if (proj%particleFlag == nucleusProjFlag .or. &
           & proj%particleFlag == pionProjFlag) then
          write (31, 3200)
       else
          write (31, 3250)
       endif
       m = 1
10     if (m <= 192) then
          if (fact*chan(m) >= test) then
             ef = chan(m)*fact
             def = sqrt(abs(chan(m)))*fact
             if (m < 192) then
                do in = m+1,192
                   if (fact*chan(in) >= test) then
                      n = in
                      efp = chan(n)*fact
                      defp = sqrt(abs(chan(n)))*fact
                      if (proj%particleFlag == nucleusProjFlag .or. &
                          & proj%particleFlag == pionProjFlag) then
                         write (31, 3400) ch(m), ef, def, ch(n), efp, defp
                      else
                         dummc = ch(m)
                         dummc(2:2) = 'g'
                         gch(m) = dummc
                         dummc = ch(n)
                         dummc(2:2) = 'g'
                         gch(n) = dummc
                         efmu = thousand*ef
                         defmu = thousand*def
                         efpmu = thousand*efp
                         defpmu = thousand*defp
                         write (31, 3450) gch(m), efmu, defmu, gch(n), efpmu, &
                              & defpmu
                      endif
                      m = n + 1
                      go to 10
                   endif
                end do
             endif
             write (31, 3500) ch(m), ef, def
          else
             m = m + 1
             go to 10
          endif
       endif
       write (31, 3300)
       do ii = 1, 13
          ij = 192 + 2*ii - 1
          ik = 192 + 2*ii
          ef = chan(ij)*fact
          def = sqrt(abs(chan(ij)))*fact
          efp = chan(ik)*fact
          defp = sqrt(abs(chan(ik)))*fact
          if (ef >= test .or. efp >= test) then
             if (proj%particleFlag == nucleusProjFlag .or. &
                 & proj%particleFlag == pionProjFlag) then
                write (31, 3400) ch(ij), ef, def, ch(ik), efp, defp
             else
                dummc = ch(ij)
                dummc(2:2) = 'g'
                gch(ij) = dummc
                dummc = ch(ik)
                dummc(2:2) = 'g'
                gch(ik) = dummc
                efmu = thousand*ef
                defmu = thousand*def
                efpmu = thousand*efp
                defpmu = thousand*defp
                write (31, 3450) gch(ij), efmu, defmu, gch(ik), efpmu, &
                     & defpmu
             endif
          endif
       end do
    endif
!   End of channelCrossSection loop ^
!
!   Print out isotope formation cross sections:
!
    if (output%nuclideCrossSection > 0) then
       abeg = targ%numBaryons + fiv
       zbeg = targ%numProtons + fiv
       izmax = nint(zbeg)
       izmax = min(izmax, 148)
       iamax = nint(abeg) + proj%numBaryons
       iamax = min(iamax, 348)
       do iz = 1,izmax + 1
          z(iz) = max (-one, zbeg - dble(iz-1))
          dasum(iz) = zro
          do ia = 1,iamax
             ar(ia) = max(abeg - dble(ia-1), zro)
             dasum(iz) = dasum(iz) + dadz(ia,iz)
          end do
          prtz(iz) = z(iz) >= zro .and. dasum(iz) > zro
       end do
       write (31, 3700)
       lprm = (izmax + 3)/3
       lprm = max(lprm, 1)
       do lpr = 1,lprm
          iz1 = 3*lpr - 2
          iz2 = 3*lpr - 1
          iz3 = 3*lpr
          prtz1 = prtz(iz3)
          prtz2 = .not.prtz(iz3) .and. prtz(iz2)
          prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. prtz(iz1)
          if (prtz1) then
             write (31, 3800) z(iz1), z(iz2), z(iz3)
          elseif (prtz2) then
             write (31, 3900) z(iz1), z(iz2)
          elseif (prtz3) then
             write (31, 4000) z(iz1)
          endif
          if ( fn < div0Lim .and. fn > -div0Lim ) then
             fn = div0Lim
             write(gsmObj%io%message, 5000) "487"
             call gsmObj%io%print(4, 3, gsmObj%io%message)
          end if
          faca = sigin/fn
          ppt = zro
          do k = 1,iamax
             y1 = zro
             y2 = zro
             y3 = zro
             y1 = dadz(k,iz1)*faca
             dy1 = sqrt(abs(dadz(k,iz1)))*faca
             y2 = dadz(k,iz2)*faca
             dy2 = sqrt(abs(dadz(k,iz2)))*faca
             y3 = dadz(k,iz3)*faca
             dy3 = sqrt(abs(dadz(k,iz3)))*faca
             printit = y1 > zro .or. y2 > zro .or. y3 > zro
             if (printit) then
!  KKG 11/15/04
                if (ar(k) > zro .and. ar(k) >= z(iz3)) then
                   iaa = nint(ar(k))
                   prtz1 = prtz(iz3)
                   prtz2 = .not.prtz(iz3) .and. prtz(iz2)
                   prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. &
                        & prtz(iz1)
                   if (prtz1) then
                      ppt = ppt + one
                      write (31, 4100) iaa, y1, dy1, y2, dy2, y3, dy3
                   elseif (prtz2) then
                      ppt = ppt + one
                      write (31, 4200) iaa, y1, dy1, y2, dy2
                   elseif (prtz3) then
                      ppt = ppt + one
                      write (31, 4300) iaa, y1, dy1
                   endif
                endif
             endif
          end do
          if (ppt > zro) then
             ipt = nint(ppt)
             ys1  = dasum(iz1)*faca
             dys1 = sqrt(abs(dasum(iz1)))*faca
             ys2  = dasum(iz2)*faca
             dys2 = sqrt(abs(dasum(iz2)))*faca
             ys3  = dasum(iz3)*faca
             dys3 = sqrt(abs(dasum(iz3)))*faca
             if (ys3 > zro) then
                write (31, 4400) ipt, ys1, dys1, ys2, dys2, ys3, dys3
             elseif (ys2 > zro) then
                write (31, 4500) ipt, ys1, dys1, ys2, dys2
             elseif (ys1 > zro) then
                write (31, 4600) ipt, ys1, dys1
             endif
          endif
       end do
       write (31, 4700)
       ip = nint(abeg)
       ipz = nint(zbeg)
       do i=1,ip
          iar(i) = ip - i + 1
       end do
       do i=1,ipz
          izr(i) = ipz - i + 1
       end do
       ss = faca
       write (31, 1100)
       yys   = zro
       atks  = zro
       datks = zro
       npt = 0
       do i = 1,ip
          yy = dadz(i,150)*ss
          dyy = sqrt(abs(dadz(i,150)))*ss
          yys = yys + dadz(i,150)
          if (dadz(i,150) > zro) then
             atk   = dadz(i,151)/dadz(i,150)
             atks  = atks + dadz(i,151)
             datk  = sqrt(abs(dadz(i,149)/dadz(i,150) - atk**2))
             datks = datks + dadz(i,149)
          else
             atk  = zro
             datk = zro
          endif
          printit = yy > zro
          if (printit) then
             npt = npt + 1
             write (31, 1200) iar(i), yy, dyy, atk, datk
          endif
       end do
       dyys = sqrt(abs(yys))*ss
       if (yys > zro) then
          atks  = atks/yys
          datks = sqrt(abs(datks/yys - atks**2))
       endif
       yys = yys*ss
       write (31, 1500) npt, yys, dyys, atks, datks
       write (31, 1300)
       yys = zro
       atks  = zro
       datks = zro
       npt = 0
       do i = 1,ipz+1
          yy = dadz(351,i)*ss
          dyy = sqrt(abs(dadz(351,i)))*ss
          yys = yys + dadz(351,i)
          if (dadz(351,i) > zro) then
             atk   = dadz(350,i)/dadz(351,i)
             atks  = atks + dadz(350,i)
             datk  = sqrt(abs(dadz(349,i)/dadz(351,i) - atk**2))
             datks = datks + dadz(349,i)
          else
             atk  = zro
             datk = zro
          endif
          printit = yy > zro
          if (printit) then
             npt = npt + 1
             write (31, 1400) izr(i), yy, dyy, atk, datk
          endif
       end do
       dyys = sqrt(abs(yys))*ss
       if (yys > zro) then
          atks  = atks/yys
          datks = sqrt(abs(datks/yys - atks**2))
       endif
       yys = yys*ss
       write (31, 1600) npt, yys, dyys, atks, datks
    endif
!   End of nuclideCrossSection loop ^
!   Print forward/backward isotope production correlations:
    if (output%nuclideCrossSection >= 2) &
        & call gsmObj%prtdadz (targ%numBaryons, targ%numProtons, proj%numBaryons, sigin, fn)
    if (output%nuclideCrossSection == 3) &
        & call gsmObj%prrdis (targ%numBaryons, targ%numProtons, proj%numBaryons, fn)
    if (output%nuclideCrossSection == 3 .and. outData%nfis.ne.0) &
        & call gsmObj%propan (fnfis)
    if (output%nuclideCrossSection == 3) &
        & call gsmObj%pdisnm (fn)
!
!   Print out angular distributions and/or spectra:
!
    if (output%energySpectra > 0 .or. &
         & output%doubleDiffSpectra > 0 .or. &
         & output%angularSpectra > 0) then
       if (output%angularSpectra > 0) then
          do m = 1,20
             om(m) = output%deltaTheta/two + dble(m-1)*output%deltaTheta
             tetmin = dble(m-1)*output%deltaTheta*degreeToRadian
             tetmax = dble(m)*output%deltaTheta*degreeToRadian
             dom(m) = (cos(tetmin) - cos(tetmax))*twpi
             temp = fl*dom(m)
             if ( temp < div0Lim .and. temp > -div0Lim ) then
                temp = div0Lim
                write(gsmObj%io%message, 5000) "629"
                call gsmObj%io%print(4, 3, gsmObj%io%message)
             end if
             const(m) = outData%sigom/(temp)
          end do
       endif
       if (output%energySpectra > 0 .or. output%doubleDiffSpectra > 0) then
          do m = 1, ntt
             temp = fl*dtt(m)
             if ( temp < div0Lim .and. temp > -div0Lim ) then
                temp = div0Lim
                write(gsmObj%io%message, 5000) "636"
                call gsmObj%io%print(4, 3, gsmObj%io%message)
             end if
             se(m) = outData%sigom/(temp)
          end do
          if (output%doubleDiffSpectra > 0) then
             do j = 1,ntt
                do i = 1,ntet
                   temp = twpi*(cos(output%angleBins(i)%lowerBound*degreeToRadian) - &
                        & cos(output%angleBins(i)%upperBound*degreeToRadian))
                   if ( temp < div0Lim .and. temp > -div0Lim ) then
                      temp = div0Lim
                      write(gsmObj%io%message, 5000) "643"
                      call gsmObj%io%print(4, 3, gsmObj%io%message)
                   end if
                   dtdo(j,i) = se(j)/(temp)
                end do
             end do
          endif
       endif
       do i = 1,ntet
          temp = fl*twpi*(cos(output%angleBins(i)%lowerBound*degreeToRadian) - &
               & cos(output%angleBins(i)%upperBound*degreeToRadian))
          if ( temp < div0Lim .and. temp > -div0Lim ) then
             temp = div0Lim
             write(gsmObj%io%message, 5000) "651"
             call gsmObj%io%print(4, 3, gsmObj%io%message)
          end if
          doo(i) = outData%sigom/(temp)
       end do
       if (output%energySpectra > 0 .or. output%doubleDiffSpectra > 0) then
90        s1 = spec(1,ntt) + spec(2,ntt) + spec(3,ntt) + spec(4,ntt) + &
               & spec(5,ntt) + spec(6,ntt) + spec(7,ntt) + spec(8,ntt) + &
               & spec(9,ntt)
          if (s1 == zro) then
             ntt = ntt - 1
             go to 90
          endif
       endif
       ejectileTally: do ic = output%minEjectileRange, output%maxejectileRange
          if(.not.output%spectraEjectiles(ic)) cycle ejectileTally
!         goahead = spec(ic,nt3) > zro .or. d2spec(ic,
          if (spec(ic,nt3) > zro .or. ang(ic,20) > zro) then
             if (spec(ic,nt3) > zro) write (31, 3600) par2(ic)
             call gsmObj%prtdist(output, ic, om, const, doo, siginn)
          endif
       end do ejectileTally
    endif
! LMK 06/2012
    if(output%printPISA)  call gsmObj%pisaprint(output, sigin, outData%ncas)

    !$OMP end critical
!
    return

! ======================================================================
999 format (/2x,'Number of equivalent gamma quanta = ',1pe13.6/ &
         & 2x,'Inelasic cross section per eqqv   = ',e13.6/ &
         & 2x,'Averaged absorption cross section = ',e13.6/ &
         & 2x,'Results are normalized to eqqv.' )
1000 format (/3x,'Inelastic cross section used here = ',f7.2,' mb'/ &
          & 1x,'Monte Carlo inelastic cross section = ',f7.2,' mb'/)
1100 format (/'Mass yield [mb] and the mean and variance of the ', &
          & 'kinetic energy [MeV]'/1x,'of residual nuclei:')
1200 format (1x,'A =',i4,1x,1pe9.3,' +/- ',e8.2,2x,e9.3,' +/- ',e8.2)
1300 format (/'Charge yield [mb] and the mean and variance of the ', &
          & ' kinetic energy [MeV]'/1x,'of residual nuclei:')
1400 format (1x,'Z =',i3,1x,1pe9.3,' +/- ',e8.2,2x,e9.3,' +/- ',e8.2)
1500 format (1x,'S =',i4,1x,1pe9.3,' +/- ',e8.2,2x,e9.3,' +/- ',e8.2)
1600 format (1x,'S =',i3,1x,1pe9.3,' +/- ',e8.2,2x,e9.3,' +/- ',e8.2)
1700 format (/15x,f6.1,' MeV (Z = ',i2,', A = ',i1,') + (Z = ',f3.0, &
          & ', A = ',f4.0,')' &
          & //1x,'Number of inelastic interactions = ',i8,',' &
          & /1x,'Number of elastic interactions   = ',i8,',' &
          & //1x,'Reaction cross section = ',f7.2,' mb,', &
          & 1x,'Elastic cross section = ',f7.2,' mb.')
1725 format (/15x,f6.1,' GeV (Z = ',i3,', A = ',i3,') + (Z = ',f3.0, &
          & ', A = ',f4.0,')' &
          & //1x,'Number of inelastic interactions = ',i8,',' &
          & /1x,'Number of elastic interactions   = ',i8,',' &
          & //1x,'Reaction cross section = ',f7.2,' mb,', &
          & 1x,'Elastic cross section = ',f7.2,' mb.')
1750 format (/15x,f6.1,' MeV (Z = ',i2,', A = ',i1,') + (Z = ',f3.0, &
          & ', A = ',f4.0,')' &
          & //1x,'Number of inelastic interactions = ',i12,',' &
          & /1x,'Number of elastic interactions   = ',i12,',' &
          & //1x,'Reaction cross section = ',f7.2,' mb,', &
          & 1x,'Elastic cross section = ',f7.2,' mb.')
1775 format (/15x,f6.1,' GeV (Z = ',i3,', A = ',i3,') + (Z = ',f3.0, &
          & ', A = ',f4.0,')' &
          & //1x,'Number of inelastic interactions = ',i12,',' &
          & /1x,'Number of elastic interactions   = ',i12,',' &
          & //1x,'Reaction cross section = ',f7.2,' mb,', &
          & 1x,'Elastic cross section = ',f7.2,' mb.')
1800 format (/15x,f6.1,' MeV (Z = ',i3,', A = ',i3,') + (Z = ',f3.0, &
          & ', A = ',f4.0,')')
1825 format (/15x,f6.1,' GeV (Z = ',i3,', A = ',i3,') + (Z = ',f3.0, &
          & ', A = ',f4.0,')')
1850 format (/15x,f6.1,' TeV (Z = ',i3,', A = ',i3,') + (Z = ',f3.0, &
          & ', A = ',f4.0,')')
1875 format (/15x,es9.2,' MeV (Z = ',i3,', A = ',i3,') + (Z = ',f3.0, &
          & ', A = ',f4.0,')')
1900 format (5x,a24)
2000 format (/15x,'Direct Monte Carlo Simulation Method:'/ &
          & 1x,'Fissility = ',f7.4,' +/- ',f7.4,','/1x,'Fission ', &
          & 'cross section = ',1pe11.5,' +/- ',e8.2,' mb.')
2100 format (/15x,'Statistical Weight Functions Method:'/ &
          & 1x,'Fissility = ',f7.4,','/1x,'Fission cross section = ', &
          & 1pe11.5,' mb.')
2200 format (/2x,'The mean excitation energy, charge, mass, and ', &
          & 'angular momentum'/2x,'of the',i9,' nuclei after the ', &
          & /2x,'cascade and before preequilibrium decay are:'/2x, &
          & 'E*av =', f6.1,' +/-',f5.1, &
          & ' MeV;  E*min =',f6.1,';  E*max =',f6.1, &
          & /2x,' Zav = ',f5.1,' +/- ',f4.1,';       Zmin = ',f4.0, &
          & ';    Zmax = ',f4.0/2x,' Aav = ',f5.1,' +/- ',f4.1,';',7x, &
          & 'Amin = ',f4.0,';    Amax = ',f4.0 &
          & /2x,' Lav =  ',f4.1,' +/- ',f4.1,' h-bar; Lmin = ',f4.0, &
          & ';    Lmax = ',f4.0)
2300 format (/2x,'The code entered PRECOF with Z < 3 or A < 5 ', &
          & i7,' times.')
2400 format (/2x,'The mean charge, mass, and angular momentum '/2x, &
          & 'of the ',i7,' residual nuclei with less than ', &
          & /2x,'3 MeV of excitation energy after the cascade are:' &
          & /2x,' Zav = ',f5.1,' +/-',f5.1,';       Zmin = ',f4.0, &
          & ';   Zmax = ',f4.0/2x,' Aav = ',f5.1,' +/-',f5.1,';',7x, &
          & 'Amin = ',f4.0,';   Amax = ',f4.0 &
          & /2x,' Lav =  ',f4.1,' +/-',f5.1,' h-bar; Lmin = ',f4.0, &
          & ';   Lmax = ',f4.0)
2500 format (/2x,'The mean excitation energy, charge, mass, ', &
          & 'angular momentum, and '/2x,'fission barrier height of ', &
          & 'the ',i7,' fissioning nuclei are:'/2x,'E*av =',f6.1, &
          & ' +/-',f5.1,' MeV; E*min =',f6.1,'; E*max = ',f6.1, &
          & /2x,' Zav = ',f5.1,' +/-',f4.1,';        Zmin = ',f4.0, &
          & ';  Zmax = ',f4.0/2x,' Aav = ',f5.1,' +/-',f4.1,';',8x, &
          & 'Amin = ',f4.0,';  Amax = ',f4.0 &
          & /2x,' Lav =  ',f4.1,' +/-',f4.1,' h-bar;  Lmin = ',f4.0, &
          & ';  Lmax = ',f4.0 &
          & /2x,'Bfav =  ',f4.1,' +/-',f4.1,' MeV;   Bfmin = ',f4.1, &
          & '; Bfmax = ',f4.1)
2600 format (/2x,'The mean kinetic energy, charge, mass, and ', &
          & 'angular momentum '/2x,'of ', &
          & 'the ',i7,' residual nuclei are:'/2x,'Ekav = ',f5.1, &
          & ' +/-',f5.1,' MeV;  Ekmin = ',f4.1,';  Ekmax =',f6.1, &
          & /2x,' Zav = ',f5.1,' +/-',f5.1,';       Zmin = ',f4.0, &
          & ';   Zmax = ',f4.0/2x,' Aav = ',f5.1,' +/-',f5.1,';',7x, &
          & 'Amin = ',f4.0,';   Amax = ',f4.0 &
          & /2x,' Lav =  ',f4.1,' +/-',f5.1,' h-bar; Lmin = ',f4.0, &
          & ';   Lmax = ',f4.0)
!2700 format (/1x,i7,' nucleus with delta-Z > 29, delta-A > 59, or ', &
!          & 'L > 99 fissioned.')
! 2800 format (/1x,i7,' nuclei with delta-Z > 29, delta-A > 59, or ', &
!           & 'L > 99 fissioned.')
!2900 format (/1x,i5,' nuclei with Z > Zinit fissioned; their average ', &
!          & 'A is ',f5.1/2x,'their average L is ',f4.1,'.')
!3000 format (/1x,i5,' nucleus with Z > Zinit fissioned; with ', &
!          & 'A = ',f5.1,' and L = ',f4.1,'.')
3100 format (/2x,'The mean excitation energy, charge, mass, and ', &
          & 'angular momentum'/2x,'of the ',i7,' nuclei after ', &
          & 'preequilibrium'/2x,'decay and before the start of ', &
          & 'statistical decay are:'/2x,'E*av =',f6.1,' +/- ',f5.1, &
          & ' MeV; E*min = ',f5.1,';  E*max =',f6.1, &
          & /2x,' Zav = ',f5.1,' +/- ',f4.1,';       Zmin = ',f4.0, &
          & ';   Zmax = ',f4.0/2x,' Aav = ',f5.1,' +/- ',f4.1,';',7x, &
          & 'Amin = ',f4.0,';   Amax = ',f4.0 &
          & /2x,' Lav =  ',f4.1,' +/- ',f4.1,' h-bar; Lmin = ',f4.0, &
          & ';   Lmax = ',f4.0)
3200 format (/1x, 'Yields of different channels (with > 1 mb):'/)
3250 format (/1x, 'Yields of different channels (with > 10 micro-b):'/)
3300 format (/1x, 'Yields of several summed channels:'/)
3400 format (1x,a,': ',f7.2,' +/- ',f5.2,' mb,',4x,a,': ',f7.2, &
          & ' +/- ',f5.2,' mb,')
3450 format (1x,a,': ',f7.1,' +/- ',f5.1,' mub,',4x,a,': ',f7.1, &
          & ' +/- ',f5.1,' mub,')
3500 format (1x,a,': ',f7.2,' +/- ',f5.2,' mb.')
3600 format (/1x,'********************************** ',a,' *********', &
          & '************************')
3700 format (/1x,'*************** Nuclide yields [mb]  (zero values ', &
          & 'suppressed) *****************')
3800 format (/17x,'Z = ',f4.0,15x,'Z = ',f4.0,15x,'Z = ',f4.0)
3900 format (/17x,'Z = ',f4.0,15x,'Z = ',f4.0)
4000 format (/17x,'Z = ',f4.0)
4100 format (1x,'A =',i4,3(1x,1pe9.3,' +/- ',e8.2)) !                                                                 
4200 format (1x,'A =',i4,2(1x,1pe9.3,' +/- ',e8.2))
4300 format (1x,'A =',i4,1(1x,1pe9.3,' +/- ',e8.2))
4400 format (1x,'S =',i4,3(1x,1pe9.3,' +/- ',e8.2))
4500 format (1x,'S =',i4,2(1x,1pe9.3,' +/- ',e8.2))
4600 format (1x,'S =',i4,1(1x,1pe9.3,' +/- ',e8.2))
4700 format (/1x,'End of nuclide yields.')
4800 format (/1x,'The program called Fermi breakup ',i8,' times.')
4900 format (/1x,'The mean total fission product kinetic energy ', &
          & 'after neutron emission is ',f6.2,' MeV.')
5000 format("Divide by zero error prevented in 'typeout.f90', line(s) ", A)
! ======================================================================
  end subroutine typeout
