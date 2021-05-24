
  subroutine initial (gsmObj, proj, output)

! ======================================================================
!
!   CEM95 written by S. G. Mashnik
!
!   Initial contains the initialization of various constants and
!   variables originally set in the main program of cem95.  This is in
!   order to shorten the program, and to remove functionality to a
!   single-purpose routine.
!
! ======================================================================
!
!    Written by A. J. Sierk,  LANL  T-2  May, 1996.
!    Modified by AJS, January-March, 1999.
!    Modified by SGM at 06/16/2000
!   "Last" change: 12-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by K. K. Gudima, October-November, 2004.
!    Edited by AJS, January, 2005.
!    Modified by KKG, February, 2005.
!    Modified : 16-Dec-2005 by AJS, following suggestion of REP.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, Sept. 2016 (Added LAQGSM variables)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, int64, real64
    use gsm_params, only: zro, one, two, thr, thousand

    implicit none
    class(GSM),           intent(inout) :: gsmObj
    class(GSMProjectile), intent(in   ) :: proj
    class(GSMOutput),     intent(inout) :: output

    integer(int32) :: i, j, j1, j2, j3, k, l
    real(real64)   :: temp, tkmax
    logical        :: log1, log2
    logical, save  :: first = .TRUE.

! ======================================================================

    real(real64) :: te1, te2, dtt, se, dtdo, d2spec, d2spe
    integer(int32) :: ntt, ntet, nt2, nt3, nti
    common /d2sdto/  te1(200), te2(200), dtt(200), &
         & se(200), dtdo(200,10), d2spec(9,10,200), &
         & d2spe(14,10,200), ntt, ntet, nt2, nt3, nti(4)
    real(real64) :: eestt, eestsq, aeqtot, aeqsq, zeqtot, zeqsq, &
         & aemin, aemax, eletot, elesq, elemin, elemax, &
         & eestrmn, eestrmx, zemin, zemax
    integer(int64) :: neq
    common /eqidat/  eestt, eestsq, aeqtot, aeqsq, zeqtot, zeqsq, &
         & aemin, aemax, eletot, elesq, elemin, elemax, &
         & eestrmn, eestrmx, zemin, zemax, neq
    real(real64) :: estart, estarsq, estrmn, estrmx, atot, atsq, &
         & ammin, ammax, eltot, elsq, elmmin, elmmax, ztot, &
         & ztsq, zmmin, zmmax, bftot, bfsq, bfmin, bfmax
    common /fisda2/  estart, estarsq, estrmn, estrmx, atot, atsq, &
         & ammin, ammax, eltot, elsq, elmmin, elmmax, ztot, &
         & ztsq, zmmin, zmmax, bftot, bfsq, bfmin, bfmax
    real(real64) :: pmul, dpmul, y, dy, emean, pmulc, dpmulc, yc, dyc, &
         & emeanc, pmulpf, dpmulpf, ypf, dypf, emeanpf, pmulsp, dpmulsp, &
         & ysp, dysp, emeansp, pmulf, dpmulf, yf, dyf, emeanf, pmcoa, &
         & dpmcoa, ycoal, dycoal, emeanco
    common /multip/ pmul(9),   dpmul(9),   y(9),   dy(9),   emean(9),  &
         & pmulc(14), dpmulc(14), yc(14), dyc(14), emeanc(14), &
         & pmulpf(6), dpmulpf(6), ypf(6), dypf(6), emeanpf(6), &
         & pmulsp(6), dpmulsp(6), ysp(6), dysp(6), emeansp(6), &
         & pmulf(6),  dpmulf(6),   yf(6),  dyf(6),  emeanf(6), &
         & pmcoa(4), dpmcoa(4), ycoal(4), dycoal(4), emeanco(4)
    integer(int32) :: nevtype
    common /nevtyp/  nevtype
    real(real64) :: epstt, epstsq, apqtot, apqsq, zpqtot, zpqsq, &
         & apmin, apmax, elptot, elpsq, elpmin, elpmax, &
         & epstrmn, epstrmx, zpmin, zpmax
    integer(int64) ::npreq
    common /predat/  epstt, epstsq, apqtot, apqsq, zpqtot, zpqsq, &
         & apmin, apmax, elptot, elpsq, elpmin, elpmax, &
         & epstrmn, epstrmx, zpmin, zpmax, npreq
    real(real64) :: speco, anco, d2speco
    common /rescoa/  speco(4,200), anco(4,20), d2speco(4,10,200)
    real(real64) :: erkt, erksq, artot, arsq, zrtot, zrsq, &
         & armin, armax, elrtot, elrsq, elrmin, elrmax, &
         & erkmn, erkmx, zrmin, zrmax
    integer(int64) :: nres
    common /resdat/  erkt, erksq, artot, arsq, zrtot, zrsq, &
         & armin, armax, elrtot, elrsq, elrmin, elrmax, &
         & erkmn, erkmx, zrmin, zrmax, nres
    real(real64) :: spef, anf, d2spef
    common /resfis/  spef(6,200), anf(6,20), d2spef(6,10,200)
    real(real64) :: spe, an
    common /resmac/  spe(14,200), an(14,20)
    real(real64) :: specsp, angsp, d2spesp
    common /respal/  specsp(6,200), angsp(6,20), d2spesp(6,10,200)
    real(real64) :: specpf, angpf, d2spepf
    common /resprf/  specpf(6,200), angpf(6,20), d2spepf(6,10,200)
    real(real64) :: spec, ang, chan, dadz
    common /result/  spec(9,200), ang(9,20), chan(218), dadz(351,151)
    real(real64) :: arttot, artsq, artmin, artmax, zrttot, zrtsq, zrtmin, &
         & zrtmax, elrttot, elrtsq, elrtmin, elrtmax
    integer(int64) :: nret
    common /retdat/  arttot, artsq, artmin, artmax, zrttot, zrtsq,  &
         & zrtmin, zrtmax, elrttot, elrtsq, elrtmin,  &
         & elrtmax, nret
    integer(int32) :: istp
    common /stopr/   istp
    real(real64) :: totke
    common /tkftot/  totke
    !   KKG 11/13/04:
    real(real64) :: dadz4
    common /result1/ dadz4(4,351,151)
    real(real64) :: rdis, dex, dpm
    common /resdis/  rdis(5,5,250), dex, dpm
    real(real64) :: opan, dth12
    common /fisopa/  opan(7,185), dth12
    real(real64) :: disnm
    common /disnmu/  disnm(6,155)

! ======================================================================

    ! Only allow one thread to set initial output data:
    !$OMP critical

    do i = 1,9
       pmul(i) = zro
       dpmul(i) = zro
       y(i) = zro
       dy(i) = zro
       emean(i) = zro
    end do

    do k = 1,200
       do j = 1,9
          spec(j,k) = zro
          if (j <= 6) then
             spef(j,k) = zro
             specpf(j,k) = zro
             specsp(j,k) = zro
          endif
          if (j <= 4) speco(j,k) = zro
          do i = 1,10
             d2spec(j,i,k) = zro
             if (j <= 4) d2speco(j,i,k) = zro
             if (j <= 6) then
                d2spef(j,i,k) = zro
                d2spepf(j,i,k) = zro
                d2spesp(j,i,k) = zro
             endif
          end do
       end do
       do j = 1,14
          spe(j,k) = zro
          do i = 1,10
             d2spe(j,i,k) = zro
          end do
       end do
    end do

    do i = 1,4
       pmcoa(i) = zro
       dpmcoa(i) = zro
       ycoal(i) = zro
       dycoal(i) = zro
       emeanco(i) = zro
       temp = output%energyBinSubStep(i)
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          print *, 'divide by zero error in initial.f90 line 311'
       end if
       nti(i) = nint((output%energyBins(i)%upperBound - &
          & output%energyBins(i)%lowerBound)/temp)
    end do
    do j = 1,10
       if (output%angleBins(j)%lowerBound < zro) then
          ntet = j - 1
          go to 10
       endif
    end do
    ntet = 10

10  do k = 1,20
       do j = 1,9
          ang(j,k) = zro
       end do
       do j = 1,14
          an(j,k) = zro
          if (j <= 6) then
             anf(j,k) = zro
             angsp(j,k) = zro
             angpf(j,k) = zro
          endif
          if (j <= 4) anco(j,k) = zro
       end do
    end do

    do k = 1,218
       chan(k) = zro
    end do

    do k = 1,151
       do j = 1,351
          dadz(j,k) = zro
          do l=1,4
             dadz4(l,j,k)=zro
          enddo
       end do
    end do

    do i=1,5
       do j=1,5
          do k=1,250
             rdis(i,j,k) = zro
          enddo
       enddo
    enddo
    dex = 10.d0
    dpm = 10.d0
    if (proj%kinEnergy > 0.7d0) dpm = 20.d0
    if (proj%kinEnergy > 1.5d0) dpm = 40.d0
    if (proj%kinEnergy > 2.5d0) dpm = 50.d0
    if (proj%kinEnergy > 1.1d0) dex = 20.d0
    do i=1,7
       do j=1,185
          opan(i,j) = zro
       enddo
    enddo
    dth12 = one

    do  i=1,6
       do  j=1,155
          disnm(i,j) = zro
       enddo
    enddo

    do i = 1,14
       pmulc(i) = zro
       dpmulc(i) = zro
       yc(i) = zro
       dyc(i) = zro
       emeanc(i) = zro
    end do

    do i = 1,6
       pmulf(i) = zro
       dpmulf(i) = zro
       yf(i) = zro
       dyf(i) = zro
       emeanf(i) = zro
       pmulpf(i) = zro
       dpmulpf(i) = zro
       ypf(i) = zro
       dypf(i) = zro
       emeanpf(i) = zro
       pmulsp(i) = zro
       dpmulsp(i) = zro
       ysp(i) = zro
       dysp(i) = zro
       emeansp(i) = zro
    end do

    tkmax = proj%kinEnergy * thousand + 30.d0
    if (proj%particleFlag == pionProjFlag) tkmax = tkmax + 110.d0
    do j = 1,199
       j1 = j - nti(1)
       j2 = j1 - nti(2)
       j3 = j2 - nti(3)
       if (j1 <= 0) then
          te1(j) = output%energyBins(1)%lowerBound + &
               & output%energyBinSubStep(1)*dble(j-1)
          te2(j) = te1(j) + output%energyBinSubStep(1)
       elseif (j2 <= 0) then
          te1(j) = output%energyBins(2)%lowerBound + &
               & output%energyBinSubStep(2)*dble(j1-1)
          te2(j) = te1(j) + output%energyBinSubStep(2)
       elseif (j3 <= 0) then
          te1(j) = output%energyBins(3)%lowerBound + &
               & output%energyBinSubStep(3)*dble(j2-1)
          te2(j) = te1(j) + output%energyBinSubStep(3)
       else
          te1(j) = output%energyBins(4)%lowerBound + &
               & output%energyBinSubStep(4)*dble(j3-1)
          te2(j) = te1(j) + output%energyBinSubStep(4)
       endif
       dtt(j) = te2(j) - te1(j)

!  KKG 10/13/04; AJS 02/15/05
       log1 = te1(j) >= tkmax
       temp = 0.0_real64
       if (proj%particleFlag == bremsProjFlag) temp = proj%brems%tMax()
       log2 = te1(j) >= (temp * thousand + 30.d0)
       if ((proj%particleFlag /= bremsProjFlag .and. log1) .or. &
            & (proj%particleFlag == bremsProjFlag .and. log2)) then
          ntt = min(j, 198)
          nt2 = ntt + 1
          nt3 = ntt + 2
          go to 20
       endif
    end do
    ntt = 198
    nt2 = 199
    nt3 = 200

20  output%fisonly = .false.
    if (output%multiplicities == 2) output%fisonly = .true.
    if (output%energySpectra == 2) output%fisonly = .true.
    if (output%doubleDiffSpectra == 2) output%fisonly = .true.
    if (output%angularSpectra == 2) output%fisonly = .true.
!    iz0 = 0
    istp = 0
    neq = 0
    npreq = 0
    nres = 0
    nret = 0
!    azro = zro
    bftot = zro
    bfsq = zro
    bfmin = 100.d0
    bfmax = zro
!    elzro = zro
    estart = zro
    estarsq = zro
    estrmn = thousand
    estrmx = zro
    eestt = zro
    eestsq = zro
    eestrmn = thousand
    eestrmx = zro
    epstt = zro
    epstsq = zro
    epstrmn = thousand
    epstrmx = zro
    erkt = zro
    erksq = zro
    erkmn = thousand
    erkmx = zro
    atot = zro
    atsq = zro
    ammax = zro
    ammin = 400.d0
    aemin = 400.d0
    aemax = zro
    aeqtot = zro
    aeqsq = zro
    apmin = 400.d0
    apmax = zro
    apqtot = zro
    apqsq = zro
    armin = 400.d0
    armax = zro
    artot = zro
    arsq = zro
    artmin = 400.d0
    artmax = zro
    arttot = zro
    artsq = zro
    elsq = zro
    eltot = zro
    elmmax = zro
    elmmin = 100.d0
    eletot = zro
    elesq = zro
    elemin = 100.d0
    elemax = zro
    elptot = zro
    elpsq = zro
    elpmin = 100.d0
    elpmax = zro
    elrtot = zro
    elrsq = zro
    elrmin = 100.d0
    elrmax = zro
    elrttot = zro
    elrtsq = zro
    elrtmin = 100.d0
    elrtmax = zro
    totke = zro
    ztot = zro
    ztsq = zro
    zmmax = zro
    zmmin = 100.d0
    zeqtot = zro
    zeqsq = zro
    zemin = 100.d0
    zemax = zro
    zpqtot = zro
    zpqsq = zro
    zpmin = 100.d0
    zpmax = zro
    zrtot = zro
    zrsq = zro
    zrmin = 100.d0
    zrmax = zro
    zrttot = zro
    zrtsq = zro
    zrtmin = 100.d0
    zrtmax = zro

    !$OMP end critical

    return
  end subroutine initial
