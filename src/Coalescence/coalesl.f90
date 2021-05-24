
  subroutine  coalesl (coalObj, results, maxA2, maxA3, maxA4, maxA6, maxA7)

! ======================================================================
!
!     Written by K. Gudima; uses the coalescence model
!     to "create" high-energy d, t, He-3, and He-4 from 
!     emitted cascade nucleons;
!     Initial number of cascade particle kstart-1 is 
!     recalculated.
!
!     Version for CEM2K2  20.03.2002 (rijmin=2.2)
!    "Last" change: 14-AUG-2003 by NVMokhov
!     Modified by A. J. Sierk, LANL T-16, October, 2003.
!     Modified to change p0 values for 300 < E < 1000 MeV
!     S. G. Mashnik, July, 2005.
!     Error corrected, AJS November, 2005.
!     Edited by AJS, LANL T-2, December, 2011.
!     Edited LMK, XCP-3, July 2013 (included error protection).
!     Modified by LMK, 2015, expanded coalescence to Be7.
!     Modified by CMJ, XCP-3, July 2018 (Coalescence class creation).
!     Modified by CMJ, XCP-3, Jan. 2019 (utilize less memory, ensure arrays don't overflow)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use coalescenceParams, only : one, two, four, avg_mass

    implicit none
    class(Coalescence),       intent(inout) :: coalObj
    type(coalescenceResults), intent(inout) :: results
    integer(int32),           intent(in   ) :: maxA2
    integer(int32),           intent(in   ) :: maxA3
    integer(int32),           intent(in   ) :: maxA4
    integer(int32),           intent(in   ) :: maxA6
    integer(int32),           intent(in   ) :: maxA7

    integer(int32) :: i, iq, j, &
         & i1, i2, i3, i4, i5, i6, i7, ia, ib, ib1, id, it, itt, &
         & numDeut, numDeut1, numTrit, numAlpha, numHe6, numLi6, numA7
    real(real64)   :: &
         & cf1, cf2, cf3, cf4, cf5, cf6, cf7, cfj, &
         & ct1, ct2, ct3, ct4, ct5, ct6, ct7, ctj, &
         &  e1,  e2,  e3,  e4,  e5,  e6,  e7,  ej, es, &
         & sf1, sf2, sf3, sf4, sf5, sf6, sf7, sfj, &
         & st1, st2, st3, st4, st5, st6, st7, stj, &
         & t1s, t2s, t3s, t4s, t5s, t6s, t7s, tjs, &
         & temp, wi, wj
    integer(int32), dimension(7) :: ind=0_int32
    real(real64),   dimension(3) :: &
         & p1l=zro, p1s=zro, &
         & p2l=zro, p2s=zro, &
         & p3l=zro, p3s=zro, &
         & p4l=zro, p4s=zro, &
         & p5l=zro, p5s=zro, &
         & p6l=zro, p6s=zro, &
         & p7l=zro, p7s=zro, &
         & pjl=zro, pjs=zro, &
         &   v=zro

! ======================================================================

    ! Variables for calculation of coalescence particles
    real(real64),   parameter   :: he6Reduction = four       ! Reduce He6 production by a factor of 4 (abnormally high rates otherwise)
    real(real64),   parameter   :: p0Factor = two * avg_mass ! Factor for nucleon comparisons
! rmin = 1/p0 = 1/5.06/(0.090GeV/c):
!    real(real64), parameter :: rijmin= (1 / 5.06) / defaultCRDeut(1)


    ! Coalescence of [D] particles
    integer(int32), dimension(maxA2) :: id1, id2
    ! Coalescence of [T/He3] particles
    integer(int32), dimension(maxA3) :: it1, it2, it3
    ! Coalescence of [He4] particles
    integer(int32), dimension(maxA4) :: ia1, ia2, ia3, ia4
    ! [ ----- EXPANDED COALESCENCE ----- ]
    ! Coalescence of [He6] particles
    integer(int32), dimension(maxA6) :: ihe61, ihe62, ihe63, ihe64, &
         & ihe65, ihe66
    ! Coalescence of [Li6] particles 
    integer(int32), dimension(maxA6) :: ili61, ili62, ili63, ili64, &
         & ili65, ili66
    ! Coalescence of [A=7] particles
    integer(int32), dimension(maxA7) :: i71, i72, i73, i74, i75, &
         & i76, i77
! ======================================================================

! ----------------
! Variable setup
! ---------------
    i = 1       ! Index of first nucleon to consider
    numDeut = 0      ! Number of formed deuterons
! LMK, 12/2014
    numA7 = 0      ! Number of formed A=7 particles
    numLi6 = 0    ! Number of formed Li6 fragments
    numHe6 = 0    ! Number of formed He6 fragments


    ! ---------------------------------------
    ! Start of main calculation:
    ! ---------------------------------------
10  continue
    wi = one

! (LMK) Only coalesce nucleons (filters all other INC fragments out)
    if (results%partBnk(i)%numBaryons  /=  1 .or. &
         & results%partBnk(i)%coalesceNum /= 0 .or. &
         & results%partBnk(i)%strangeness /= 0 .or. wi < one) then
       if (i > (results%numParticles-1)) go to 30
       i = i + 1
       go to 10
    endif
    j = i + 1

! Find allowed interacting nucleon
20  if (j > results%numParticles) then
       if (i > (results%numParticles-1)) go to 30
       i = i + 1
       go to 10
    endif
    wj = one
! (LMK) Only coalesce nucleons
    if (results%partBnk(j)%numBaryons /= 1 .or. results%partBnk(j)%coalesceNum /= 0 &
         & .or. results%partBnk(j)%strangeness /= 0 .or. wj < one) then
       j = j + 1
       go to 20
    endif

! Obtain nucleon information
    i1 = i
    e1 = results%partBnk(i1)%kinEnergy + avg_mass
    ej = results%partBnk(j)%kinEnergy + avg_mass
    pjl(1) = results%partBnk(j)%linearMomX
    pjl(2) = results%partBnk(j)%linearMomY
    pjl(3) = results%partBnk(j)%linearMomZ
    es = ej + e1
    p1l(1) = results%partBnk(i1)%linearMomX
    p1l(2) = results%partBnk(i1)%linearMomY
    p1l(3) = results%partBnk(i1)%linearMomZ
    temp = es
    if (temp < div0Limit .and. temp > -div0Limit) then
       temp = div0Limit
       write(coalObj%io%message, 2000) "249-251"
       call coalObj%io%print(4, 3, coalObj%io%message)
    end if
    v(1) = -(p1l(1) + pjl(1))/temp
    v(2) = -(p1l(2) + pjl(2))/temp
    v(3) = -(p1l(3) + pjl(3))/temp
    call coalObj%kinema (pjl, v, pjs, ctj, stj, cfj, sfj, tjs, avg_mass)
!    call rijm (i, j, rmin) !! Coalescence will ignore separation
!    if (rmin > rijmin .or. & !! in r-space; only consider p-space!
!         & sqrt(tjs*(tjs + p0Factor)) > coalObj%data%coalesRadiiDeut()) then !! SGM 07/07/05
!       j = j + 1
    temp = tjs*(tjs + p0Factor)
    if (temp < 0.0d0) then
       temp = sqrRootCorrection
       write(coalObj%io%message, 2100) "262"
       call coalObj%io%print(4, 3, coalObj%io%message)
    end if

    ! Check if partner's momentum is small enough
    if (sqrt(temp) > coalObj%data%coalesRadiiDeut()) then
       ! Participant won't coalesce, find new partner
       j = j + 1
       go to 20
    endif

    results%numCoalesced = results%numCoalesced + 2 ! 2 nucleons coalesed into a 'd' (n+p)
    ! if space remains in the internal array, add a deuteron (coalesced)
    if (numDeut < maxA2 ) then
! (LMK) Coalesce a deuteron
       numDeut = numDeut + 1
       id1(numDeut) = i
       id2(numDeut) = j
       results%partBnk(i)%coalesceNum = numDeut
       results%partBnk(j)%coalesceNum = numDeut
       if (i <= results%numParticles-1) then
          ! Find another set of nucleons to coalesce
          i = i + 1
          go to 10
       endif
    else
       ! Arrays exceeded!
       write(coalObj%io%message, 1500) "deuteron"
       call coalObj%io%print(1, 3, coalObj%io%message)
    endif


! Check for coalesced alpha particles (from the coalesced deuterons, d+d)
30  numAlpha = 0
    if (numDeut >= 2) then
       numDeut1 = numDeut - 1
       do id = 1,numDeut1
          ! Obtain index of deuteron's nucleons (not yet coalesced)
          i1 = id1(id)
          i2 = id2(id)

          if (i1 == 0) go to 50 ! No more deuterons to coalesce

          ! Look for another deuteron check if can coalesce with the previously found one
          ib1 = id + 1
          do ib = ib1,numDeut
             ! Obtain deuteron index
             i3 = id1(ib)
             j = id2(ib)

             if (i3 == 0) go to 40 ! No more deuterons to coalesce with

             ! Obtain the fragment charge, check if charge is that of helium
             iq = results%partBnk(i1)%charge + results%partBnk(i2)%charge + &
                  & results%partBnk(i3)%charge + results%partBnk(j)%charge
             if (iq /= 2) go to 40
             ! Obtain energy and momentum of each nucleons
             e1 = results%partBnk(i1)%kinEnergy + avg_mass
             e2 = results%partBnk(i2)%kinEnergy + avg_mass
             e3 = results%partBnk(i3)%kinEnergy + avg_mass
             ej = results%partBnk(j)%kinEnergy + avg_mass
             p1l(1) = results%partBnk(i1)%linearMomX
             p1l(2) = results%partBnk(i1)%linearMomY
             p1l(3) = results%partBnk(i1)%linearMomZ
             p2l(1) = results%partBnk(i2)%linearMomX
             p2l(2) = results%partBnk(i2)%linearMomY
             p2l(3) = results%partBnk(i2)%linearMomZ
             p3l(1) = results%partBnk(i3)%linearMomX
             p3l(2) = results%partBnk(i3)%linearMomY
             p3l(3) = results%partBnk(i3)%linearMomZ
             pjl(1) = results%partBnk(j)%linearMomX
             pjl(2) = results%partBnk(j)%linearMomY
             pjl(3) = results%partBnk(j)%linearMomZ
             es = ej + e1 + e2 + e3 ! Obtain total energy
             temp = es
             if (temp < div0Limit .and. temp > -div0Limit) then
                temp = div0Limit
                write(coalObj%io%message, 2000) "317-319"
                call coalObj%io%print(4, 3, coalObj%io%message)
             end if
             v(1) = -(p1l(1) + p2l(1) + p3l(1) + pjl(1))/temp
             v(2) = -(p1l(2) + p2l(2) + p3l(2) + pjl(2))/temp
             v(3) = -(p1l(3) + p2l(3) + p3l(3) + pjl(3))/temp
             call coalObj%kinema (p1l, v, p1s, ct1, st1, cf1, sf1, t1s, avg_mass)
             call coalObj%kinema (p2l, v, p2s, ct2, st2, cf2, sf2, t2s, avg_mass)
             call coalObj%kinema (p3l, v, p3s, ct3, st3, cf3, sf3, t3s, avg_mass)
             call coalObj%kinema (pjl, v, pjs, ctj, stj, cfj, sfj, tjs, avg_mass)

! Check if nucleons can coalesce
! (LMK, Altered comparison to eliminate LHS sqrt)
             if (t1s*(t1s + p0Factor) > coalObj%data%coalesRadiiAlphaSqrd()) go to 40
             if (t2s*(t2s + p0Factor) > coalObj%data%coalesRadiiAlphaSqrd()) go to 40
             if (t3s*(t3s + p0Factor) > coalObj%data%coalesRadiiAlphaSqrd()) go to 40
             if (tjs*(tjs + p0Factor) > coalObj%data%coalesRadiiAlphaSqrd()) go to 40

             if (numAlpha >= maxA4 ) then
                ! Arrays exceeded!
                write(coalObj%io%message, 1500) "alpha"
                call coalObj%io%print(1, 3, coalObj%io%message)
             else 
! (LMK) Coalesce an alpha particle
                numAlpha = numAlpha + 1
                ia1(numAlpha) = i1
                ia2(numAlpha) = i2
                ia3(numAlpha) = i3
                ia4(numAlpha) = j
                ! Set flag for coalesced nucleons (for alpha)
                results%partBnk(i1)%coalesceNum = numAlpha
                results%partBnk(i2)%coalesceNum = numAlpha
                results%partBnk(i3)%coalesceNum = numAlpha
                results%partBnk(j)%coalesceNum = numAlpha
                ! Remove deuteron flag
                id1(id) = 0
                id2(id) = 0
!  Misprint found by AJS 11/18/05
!                id1(id) = 0
!                id2(id) = 0
                id1(ib) = 0
                id2(ib) = 0
                go to 50 ! Deuterons have coalesced; look for a new pair
             endif
40           continue
          end do
50        continue
       end do
    endif

    ! Coalesce tritons/He-3 (d+n or d+p)
    numTrit = 0
    if (numDeut == 0) go to 150
    do id = 1,numDeut
       i1 = id1(id)
       i2 = id2(id)
       if (i1 /= 0) then
          do j = 1,results%numParticles
             wj = one

             if (results%partBnk(j)%numBaryons /= 1     .or. &
                  & results%partBnk(j)%coalesceNum /= 0 .or. &
                  & results%partBnk(j)%strangeness /= 0 .or. wj < one) go to 70

             iq = results%partBnk(i1)%charge + results%partBnk(i2)%charge + &
                  results%partBnk(j)%charge

             ! Get rid of 3n or 3p coalesced groups, LMK
             if (iq == 0 .or. iq == 3) go to 70

             e1 = results%partBnk(i1)%kinEnergy + avg_mass  
             e2 = results%partBnk(i2)%kinEnergy + avg_mass  
             ej = results%partBnk(j)%kinEnergy + avg_mass  
             p1l(1) = results%partBnk(i1)%linearMomX
             p1l(2) = results%partBnk(i1)%linearMomY
             p1l(3) = results%partBnk(i1)%linearMomZ
             p2l(1) = results%partBnk(i2)%linearMomX
             p2l(2) = results%partBnk(i2)%linearMomY
             p2l(3) = results%partBnk(i2)%linearMomZ
             pjl(1) = results%partBnk(j)%linearMomX
             pjl(2) = results%partBnk(j)%linearMomY
             pjl(3) = results%partBnk(j)%linearMomZ
             es = ej + e1 + e2
             temp = es
             if (temp < div0Limit .and. temp > -div0Limit) then
                temp = div0Limit
                write(coalObj%io%message, 2000) "389-391"
                call coalObj%io%print(4, 3, coalObj%io%message)
             end if
             v(1) = -(p1l(1) + p2l(1) + pjl(1))/temp
             v(2) = -(p1l(2) + p2l(2) + pjl(2))/temp
             v(3) = -(p1l(3) + p2l(3) + pjl(3))/temp
             call coalObj%kinema (p1l, v, p1s, ct1, st1, cf1, sf1, t1s, avg_mass)
             call coalObj%kinema (p2l, v, p2s, ct2, st2, cf2, sf2, t2s, avg_mass)
             call coalObj%kinema (pjl, v, pjs, ctj, stj, cfj, sfj, tjs, avg_mass) 
! LMK, Altered comparison to eliminate LHS sqrt
             if (t1s*(t1s + p0Factor) > coalObj%data%coalesRadiiTritSqrd()) go to 70
             if (t2s*(t2s + p0Factor) > coalObj%data%coalesRadiiTritSqrd()) go to 70
             if (tjs*(tjs + p0Factor) > coalObj%data%coalesRadiiTritSqrd()) go to 70

             ! One nucleon coalesced w/ a d (d+n or d+p)
             results%numCoalesced = results%numCoalesced + 1

             if (numTrit >= maxA3 ) then
                write(coalObj%io%message, 1500) "triton or helion"
                call coalObj%io%print(1, 3, coalObj%io%message)
             else
! (LMK) Coalesce a triton or 3He
                numTrit = numTrit + 1
                it1(numTrit) = i1
                it2(numTrit) = i2
                it3(numTrit) = j
                results%partBnk(i1)%coalesceNum = numTrit
                results%partBnk(i2)%coalesceNum = numTrit
                results%partBnk(j)%coalesceNum = numTrit
                id1(id) = 0
                id2(id) = 0
                go to 80
             endif
70           continue
          end do
       endif
80     continue
    end do

    ! Coalesce more He4 (t+p or He3+n)
    if (numTrit > 0) then
       do it = 1,numTrit
          do j = 1,results%numParticles
             wj = one
             if (results%partBnk(j)%numBaryons /= 1     .or. &
                  & results%partBnk(j)%coalesceNum /= 0 .or. &
                  & results%partBnk(j)%strangeness /= 0 .or. wj < one) go to 90
             i1 = it1(it)
             i2 = it2(it)
             i3 = it3(it)
             iq = results%partBnk(i1)%charge + results%partBnk(i2)%charge + & 
                  & results%partBnk(i3)%charge + results%partBnk(j)%charge
             if (iq /= 2) go to 90
! Make sure combination is t+p, or 3He+n, to coalesce 4He, LMK  ^^^^^
             e1 = results%partBnk(i1)%kinEnergy + avg_mass  
             e2 = results%partBnk(i2)%kinEnergy + avg_mass  
             e3 = results%partBnk(i3)%kinEnergy + avg_mass  
             ej = results%partBnk(j)%kinEnergy + avg_mass  
             p1l(1) = results%partBnk(i1)%linearMomX
             p1l(2) = results%partBnk(i1)%linearMomY
             p1l(3) = results%partBnk(i1)%linearMomZ
             p2l(1) = results%partBnk(i2)%linearMomX
             p2l(2) = results%partBnk(i2)%linearMomY
             p2l(3) = results%partBnk(i2)%linearMomZ
             p3l(1) = results%partBnk(i3)%linearMomX
             p3l(2) = results%partBnk(i3)%linearMomY
             p3l(3) = results%partBnk(i3)%linearMomZ
             pjl(1) = results%partBnk(j)%linearMomX
             pjl(2) = results%partBnk(j)%linearMomY
             pjl(3) = results%partBnk(j)%linearMomZ
             es = ej + e1 + e2 + e3
             temp = es
             if (temp < div0Limit .and. temp > -div0Limit) then
                temp = div0Limit
                write(coalObj%io%message, 2000) "453-455"
                call coalObj%io%print(4, 3, coalObj%io%message)
             end if
             v(1) = -(p1l(1) + p2l(1) + p3l(1) + pjl(1))/temp
             v(2) = -(p1l(2) + p2l(2) + p3l(2) + pjl(2))/temp
             v(3) = -(p1l(3) + p2l(3) + p3l(3) + pjl(3))/temp
             call coalObj%kinema (p1l, v, p1s, ct1, st1, cf1, sf1, t1s, avg_mass)
             call coalObj%kinema (p2l, v, p2s, ct2, st2, cf2, sf2, t2s, avg_mass)
             call coalObj%kinema (p3l, v, p3s, ct3, st3, cf3, sf3, t3s, avg_mass)
             call coalObj%kinema (pjl, v, pjs, ctj, stj, cfj, sfj, tjs, avg_mass)
! LMK, Altered comparison to eliminate LHS sqrt
             if (t1s*(t1s + p0Factor) > coalObj%data%coalesRadiiAlphaSqrd()) go to 90
             if (t2s*(t2s + p0Factor) > coalObj%data%coalesRadiiAlphaSqrd()) go to 90
             if (t3s*(t3s + p0Factor) > coalObj%data%coalesRadiiAlphaSqrd()) go to 90
             if (tjs*(tjs + p0Factor) > coalObj%data%coalesRadiiAlphaSqrd()) go to 90

             ! 1 Nucleon coalesced (t+p or He3+n)
             results%numCoalesced = results%numCoalesced + 1

             if (numAlpha >= maxA4 ) then
                write(coalObj%io%message, 1500) "alpha"
                call coalObj%io%print(1, 3, coalObj%io%message)
             else
! Coalesce 4He from t+p or 3He+n, LMK
                numAlpha = numAlpha + 1
                ia1(numAlpha) = i1
                ia2(numAlpha) = i2
                ia3(numAlpha) = i3
                ia4(numAlpha) = j
                results%partBnk(i1)%coalesceNum = numAlpha
                results%partBnk(i2)%coalesceNum = numAlpha
                results%partBnk(i3)%coalesceNum = numAlpha
                results%partBnk(j)%coalesceNum = numAlpha
                it1(it) = 0
                it2(it) = 0
                it3(it) = 0
                go to 100
             endif
90           continue
          end do
100       continue
       end do
    endif

! (LMK, 12/2014) Add coalescence of 6He, 6,7Li, 7Be
! [He4 + t/He3, He4 + d/2n]
    if ( coalObj%options%expandedCoalescence > 0 ) then
       do ia = 1, numAlpha
          if (ia4(ia) < 0.5) cycle   ! Check for 4He already coalesced... 
          i1 = ia1(ia)
          i2 = ia2(ia)
          i3 = ia3(ia)
          i4 = ia4(ia)
          e1 = results%partBnk(i1)%kinEnergy + avg_mass
          e2 = results%partBnk(i2)%kinEnergy + avg_mass
          e3 = results%partBnk(i3)%kinEnergy + avg_mass
          e4 = results%partBnk(i4)%kinEnergy + avg_mass       
          p1l(1) = results%partBnk(i1)%linearMomX
          p1l(2) = results%partBnk(i1)%linearMomY
          p1l(3) = results%partBnk(i1)%linearMomZ
          p2l(1) = results%partBnk(i2)%linearMomX
          p2l(2) = results%partBnk(i2)%linearMomY
          p2l(3) = results%partBnk(i2)%linearMomZ
          p3l(1) = results%partBnk(i3)%linearMomX
          p3l(2) = results%partBnk(i3)%linearMomY
          p3l(3) = results%partBnk(i3)%linearMomZ
          p4l(1) = results%partBnk(i4)%linearMomX
          p4l(2) = results%partBnk(i4)%linearMomY
          p4l(3) = results%partBnk(i4)%linearMomZ
          do it = 1, numTrit ! Coalesce 7Be, 7Li, a+t
             if (it3(it) < 0.5) cycle    ! Check for t/3He already coalesced... 
             i5 = it1(it)
             i6 = it2(it)
             i7 = it3(it)
             e5 = results%partBnk(i5)%kinEnergy + avg_mass
             e6 = results%partBnk(i6)%kinEnergy + avg_mass
             e7 = results%partBnk(i7)%kinEnergy + avg_mass
             p5l(1) = results%partBnk(i5)%linearMomX
             p5l(2) = results%partBnk(i5)%linearMomY
             p5l(3) = results%partBnk(i5)%linearMomZ
             p6l(1) = results%partBnk(i6)%linearMomX
             p6l(2) = results%partBnk(i6)%linearMomY
             p6l(3) = results%partBnk(i6)%linearMomZ
             p7l(1) = results%partBnk(i7)%linearMomX
             p7l(2) = results%partBnk(i7)%linearMomY
             p7l(3) = results%partBnk(i7)%linearMomZ
             es = e1 + e2 + e3 + e4 + e5 + e6 + e7
             if (es > 0) then
                v(1) = -(p1l(1) + p2l(1) + p3l(1) + p4l(1) + p5l(1) + p6l(1) &
                     & + p7l(1))/es
                v(2) = -(p1l(2) + p2l(2) + p3l(2) + p4l(2) + p5l(2) + p6l(2) &
                     + p7l(2))/es
                v(3) = -(p1l(3) + p2l(3) + p3l(3) + p4l(3) + p5l(3) + p6l(3) &
                     + p7l(3))/es
             else
                v(:) = 0.0d0
                write(coalObj%io%message, 2000) "538"
                call coalObj%io%print(4, 3, coalObj%io%message)
             end if
             call coalObj%kinema (p1l, v, p1s, ct1, st1, cf1, sf1, t1s, avg_mass)
             call coalObj%kinema (p2l, v, p2s, ct2, st2, cf2, sf2, t2s, avg_mass)
             call coalObj%kinema (p3l, v, p3s, ct3, st3, cf3, sf3, t3s, avg_mass)
             call coalObj%kinema (p4l, v, p4s, ct4, st4, cf4, sf4, t4s, avg_mass)
             call coalObj%kinema (p5l, v, p5s, ct5, st5, cf5, sf5, t5s, avg_mass)
             call coalObj%kinema (p6l, v, p6s, ct6, st6, cf6, sf6, t6s, avg_mass)
             call coalObj%kinema (p7l, v, p7s, ct7, st7, cf7, sf7, t7s, avg_mass)
             if (t1s*(t1s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t2s*(t2s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t3s*(t3s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t4s*(t4s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t5s*(t5s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t6s*(t6s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t7s*(t7s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 

! Coalesced 7Be or 7Li
             numA7 = numA7 + 1
             if (numA7 <= maxA7 ) then
                i71(numA7) = i1
                i72(numA7) = i2
                i73(numA7) = i3
                i74(numA7) = i4
                i75(numA7) = i5
                i76(numA7) = i6
                i77(numA7) = i7
                results%partBnk(i1)%coalesceNum = numA7
                results%partBnk(i2)%coalesceNum = numA7
                results%partBnk(i3)%coalesceNum = numA7
                results%partBnk(i4)%coalesceNum = numA7
                results%partBnk(i5)%coalesceNum = numA7
                results%partBnk(i6)%coalesceNum = numA7
                results%partBnk(i7)%coalesceNum = numA7
                ia1(ia) = 0
                ia2(ia) = 0
                ia3(ia) = 0
                ia4(ia) = 0
                it1(it) = 0
                it2(it) = 0
                it3(it) = 0
             else
                write(coalObj%io%message, 1500) "Li7 or Be7"
                call coalObj%io%print(1, 3, coalObj%io%message)
             end if
             exit 
          end do ! numTrit loop
          if (ia4(ia) < 0.5) cycle      ! If alpha was coalesced, skip to next alpha

! Coalesce 6Li or 6He, a+d
          do id = 1, numDeut
             if (id2(id) < 0.5) cycle    ! Check for d already coalesced... 
             i5 = id1(id)
             i6 = id2(id)
             ! Check for coalesced 6Be and cycle
             if ( (results%partBnk(i1)%charge + results%partBnk(i2)%charge + &
                  & results%partBnk(i3)%charge + results%partBnk(i4)%charge + &
                  & results%partBnk(i5)%charge + results%partBnk(i6)%charge) > 3.5) &
                  & cycle
             e5 = results%partBnk(i5)%kinEnergy + avg_mass
             e6 = results%partBnk(i6)%kinEnergy + avg_mass
             p5l(1) = results%partBnk(i5)%linearMomX
             p5l(2) = results%partBnk(i5)%linearMomY
             p5l(3) = results%partBnk(i5)%linearMomZ
             p6l(1) = results%partBnk(i6)%linearMomX
             p6l(2) = results%partBnk(i6)%linearMomY
             p6l(3) = results%partBnk(i6)%linearMomZ
             es = e1 + e2 + e3 + e4 + e5 + e6
             if (es > 0) then
                v(1) = -(p1l(1) + p2l(1) + p3l(1) + p4l(1) + p5l(1) + p6l(1))/es
                v(2) = -(p1l(2) + p2l(2) + p3l(2) + p4l(2) + p5l(2) + p6l(2))/es
                v(3) = -(p1l(3) + p2l(3) + p3l(3) + p4l(3) + p5l(3) + p6l(3))/es
             else
                v(:) = 0.0d0
                write(coalObj%io%message, 2000) "605"
                call coalObj%io%print(4, 3, coalObj%io%message)
             end if
             call coalObj%kinema (p1l, v, p1s, ct1, st1, cf1, sf1, t1s, avg_mass)
             call coalObj%kinema (p2l, v, p2s, ct2, st2, cf2, sf2, t2s, avg_mass)
             call coalObj%kinema (p3l, v, p3s, ct3, st3, cf3, sf3, t3s, avg_mass)
             call coalObj%kinema (p4l, v, p4s, ct4, st4, cf4, sf4, t4s, avg_mass)
             call coalObj%kinema (p5l, v, p5s, ct5, st5, cf5, sf5, t5s, avg_mass)
             call coalObj%kinema (p6l, v, p6s, ct6, st6, cf6, sf6, t6s, avg_mass)
             if (t1s*(t1s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t2s*(t2s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t3s*(t3s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t4s*(t4s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t5s*(t5s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t6s*(t6s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
! Coalesced 6Li or 6He
             if ( (results%partBnk(i1)%charge + results%partBnk(i2)%charge + &
                  & results%partBnk(i3)%charge + results%partBnk(i4)%charge + &
                  & results%partBnk(i5)%charge + results%partBnk(i6)%charge) < 2.5) then
                ! Reduce amount of 6He, take only 1/he6Reduction [originally =4]
                if ( mod(dble(i2), he6Reduction) < 2.5) cycle

             end if
             numLi6 = numLi6 + 1
             if (numLi6 <= maxA6 ) then
                ili61(numLi6) = i1
                ili62(numLi6) = i2
                ili63(numLi6) = i3
                ili64(numLi6) = i4
                ili65(numLi6) = i5
                ili66(numLi6) = i6
                results%partBnk(i1)%coalesceNum = numLi6
                results%partBnk(i2)%coalesceNum = numLi6
                results%partBnk(i3)%coalesceNum = numLi6
                results%partBnk(i4)%coalesceNum = numLi6
                results%partBnk(i5)%coalesceNum = numLi6
                results%partBnk(i6)%coalesceNum = numLi6
                ia1(ia) = 0
                ia2(ia) = 0
                ia3(ia) = 0
                ia4(ia) = 0
                id1(id) = 0
                id2(id) = 0
             else
                write(coalObj%io%message, 1500) "Li6"
                call coalObj%io%print(1, 3, coalObj%io%message)
             end if
             exit 
          end do ! numDeut loop
       end do ! numAlpha loop
       do it = 1, numTrit-1 ! Coalesce 6He or 6Li, t+t. t+He3
          if (it3(it) < 0.5) cycle      ! Check for t/3He already coalesced... 
          i1 = it1(it)
          i2 = it2(it)
          i3 = it3(it)
          e1 = results%partBnk(i1)%kinEnergy + avg_mass
          e2 = results%partBnk(i2)%kinEnergy + avg_mass
          e3 = results%partBnk(i3)%kinEnergy + avg_mass
          p1l(1) = results%partBnk(i1)%linearMomX
          p1l(2) = results%partBnk(i1)%linearMomY
          p1l(3) = results%partBnk(i1)%linearMomZ
          p2l(1) = results%partBnk(i2)%linearMomX
          p2l(2) = results%partBnk(i2)%linearMomY
          p2l(3) = results%partBnk(i2)%linearMomZ
          p3l(1) = results%partBnk(i3)%linearMomX
          p3l(2) = results%partBnk(i3)%linearMomY
          p3l(3) = results%partBnk(i3)%linearMomZ
          do itt = it+1, numTrit     
             if (it3(itt) < 0.5) cycle   ! Check for t/3He already coalesced... 
             i4 = it1(itt)
             i5 = it2(itt)
             i6 = it3(itt) 
             ! Check for coalesced 6Be and cycle
            if ( (results%partBnk(i1)%charge + results%partBnk(i2)%charge + &
                  & results%partBnk(i3)%charge + results%partBnk(i4)%charge + &
                  & results%partBnk(i5)%charge + results%partBnk(i6)%charge) &
                  & > 3.5) cycle
             e4 = results%partBnk(i4)%kinEnergy + avg_mass
             e5 = results%partBnk(i5)%kinEnergy + avg_mass
             e6 = results%partBnk(i6)%kinEnergy + avg_mass
             p4l(1) = results%partBnk(i4)%linearMomX
             p4l(2) = results%partBnk(i4)%linearMomY
             p4l(3) = results%partBnk(i4)%linearMomZ
             p5l(1) = results%partBnk(i5)%linearMomX
             p5l(2) = results%partBnk(i5)%linearMomY
             p5l(3) = results%partBnk(i5)%linearMomZ
             p6l(1) = results%partBnk(i6)%linearMomX
             p6l(2) = results%partBnk(i6)%linearMomY
             p6l(3) = results%partBnk(i6)%linearMomZ
             es = e1 + e2 + e3 + e4 + e5 + e6
             if (es > 0) then
                v(1) = -(p1l(1) + p2l(1) + p3l(1) + p4l(1) + p5l(1) + p6l(1))/es
                v(2) = -(p1l(2) + p2l(2) + p3l(2) + p4l(2) + p5l(2) + p6l(2))/es
                v(3) = -(p1l(3) + p2l(3) + p3l(3) + p4l(3) + p5l(3) + p6l(3))/es
             else
                v(:) = 0.0d0
                write(coalObj%io%message, 2000) "693"
                call coalObj%io%print(4, 3, coalObj%io%message)
             end if
             call coalObj%kinema (p1l, v, p1s, ct1, st1, cf1, sf1, t1s, avg_mass)
             call coalObj%kinema (p2l, v, p2s, ct2, st2, cf2, sf2, t2s, avg_mass)
             call coalObj%kinema (p3l, v, p3s, ct3, st3, cf3, sf3, t3s, avg_mass)
             call coalObj%kinema (p4l, v, p4s, ct4, st4, cf4, sf4, t4s, avg_mass)
             call coalObj%kinema (p5l, v, p5s, ct5, st5, cf5, sf5, t5s, avg_mass)
             call coalObj%kinema (p6l, v, p6s, ct6, st6, cf6, sf6, t6s, avg_mass)
             if (t1s*(t1s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t2s*(t2s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t3s*(t3s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t4s*(t4s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t5s*(t5s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
             if (t6s*(t6s + p0Factor) > coalObj%data%coalesRadiiLFragSqrd()) cycle 
! Coalesced 6He or 6Li
             if ( (results%partBnk(i1)%charge + results%partBnk(i2)%charge + &
                  & results%partBnk(i3)%charge + results%partBnk(i4)%charge + &
                  & results%partBnk(i5)%charge + results%partBnk(i6)%charge) < 2.5) &
                  & then
                ! Reduce amount of 6He, take only 1/he6Reduction
                if (mod(dble(i2), he6Reduction) < 2.5) cycle
             end if
             numHe6 = numHe6 + 1
             if (numHe6 <= maxA6 ) then
                ihe61(numHe6) = i1
                ihe62(numHe6) = i2
                ihe63(numHe6) = i3
                ihe64(numHe6) = i4
                ihe65(numHe6) = i5
                ihe66(numHe6) = i6
                results%partBnk(i1)%coalesceNum = numHe6
                results%partBnk(i2)%coalesceNum = numHe6
                results%partBnk(i3)%coalesceNum = numHe6
                results%partBnk(i4)%coalesceNum = numHe6
                results%partBnk(i5)%coalesceNum = numHe6
                results%partBnk(i6)%coalesceNum = numHe6
                it1(it) = 0
                it2(it) = 0
                it3(it) = 0
                it1(itt) = 0
                it2(itt) = 0
                it3(itt) = 0
             else
                write(coalObj%io%message, 1500) "He6"
                call coalObj%io%print(1, 3, coalObj%io%message)
             end if
             exit 
          end do ! numTrit loop
       end do ! numTrit loop

       do i = 1, numA7
          i1 = i71(i)
          i2 = i72(i)
          i3 = i73(i)
          i4 = i74(i)
          i5 = i75(i)
          i6 = i76(i)
          i7 = i77(i)
          ind(1) = i1
          ind(2) = i2
          ind(3) = i3
          ind(4) = i4
          ind(5) = i5
          ind(6) = i6
          ind(7) = i7
          call coalObj%codir( results, ind, 7)
          results%partBnk(i1)%charge = results%partBnk(i1)%charge + &
               & results%partBnk(i2)%charge + results%partBnk(i3)%charge + &
               & results%partBnk(i4)%charge + results%partBnk(i5)%charge + &
               results%partBnk(i6)%charge + results%partBnk(i7)%charge
          results%partBnk(i1)%coalesceNum = 0
          results%partBnk(i1)%numBaryons = 7
       end do
       do i = 1, numLi6
          i1 = ili61(i)
          i2 = ili62(i)
          i3 = ili63(i)
          i4 = ili64(i)
          i5 = ili65(i)
          i6 = ili66(i)
          ind(1) = i1
          ind(2) = i2
          ind(3) = i3
          ind(4) = i4
          ind(5) = i5
          ind(6) = i6
          call coalObj%codir( results, ind, 6)
          results%partBnk(i1)%charge = results%partBnk(i1)%charge + &
               & results%partBnk(i2)%charge + results%partBnk(i3)%charge + &
               & results%partBnk(i4)%charge + results%partBnk(i5)%charge + &
               & results%partBnk(i6)%charge
          results%partBnk(i1)%coalesceNum = 0
          results%partBnk(i1)%numBaryons = 6
       end do
       do i = 1, numHe6
          i1 = ihe61(i)
          i2 = ihe62(i)
          i3 = ihe63(i)
          i4 = ihe64(i)
          i5 = ihe65(i)
          i6 = ihe66(i)
          ind(1) = i1
          ind(2) = i2
          ind(3) = i3
          ind(4) = i4
          ind(5) = i5
          ind(6) = i6
          call coalObj%codir( results, ind, 6)
          results%partBnk(i1)%charge = results%partBnk(i1)%charge + &
               & results%partBnk(i2)%charge + results%partBnk(i3)%charge + &
               & results%partBnk(i4)%charge + results%partBnk(i5)%charge + &
               & results%partBnk(i6)%charge
          results%partBnk(i1)%coalesceNum = 0
          results%partBnk(i1)%numBaryons = 6
       end do
    end if
! end expanded coalescence



! Coalesce remaining flagged d, t/He3, and He4
    if (numAlpha > 0) then
       ! Coalesce all alpha fragments
       do ia = 1,numAlpha
          if (ia1(ia) > 0) then
             i1 = ia1(ia)
             i2 = ia2(ia)
             i3 = ia3(ia)
             i4 = ia4(ia)
             ind(1) = i1
             ind(2) = i2
             ind(3) = i3
             ind(4) = i4
             call coalObj%codir( results, ind, 4)
             results%partBnk(i1)%charge = 2 
             results%partBnk(i1)%coalesceNum = 0
             results%partBnk(i1)%numBaryons = 4
          end if
       end do
    endif
    if (numTrit /= 0) then
! --------------------------------------
       ! Coalesce all tritium fragments
       do it = 1,numTrit
          if (it1(it) > 0) then
             i1 = it1(it)
             i2 = it2(it)
             i3 = it3(it)
             ind(1) = i1
             ind(2) = i2
             ind(3) = i3
             call coalObj%codir( results, ind, 3)
             results%partBnk(i1)%charge = results%partBnk(i1)%charge + &
                  & results%partBnk(i2)%charge + results%partBnk(i3)%charge
             results%partBnk(i1)%coalesceNum = 0
             results%partBnk(i1)%numBaryons = 3
          endif
       end do
    endif
! --------------------------------------
    do id = 1,numDeut
       ! Coalesce all deuterium fragments
       if (id1(id) > 0) then
          i1 = id1(id)
          i2 = id2(id)
          if ((results%partBnk(i1)%charge + results%partBnk(i2)%charge) == 1) then
             ind(1) = i1
             ind(2) = i2
             call coalObj%codir( results, ind, 2)
             results%partBnk(i1)%charge = 1
             results%partBnk(i1)%coalesceNum = 0
             results%partBnk(i1)%numBaryons = 2
          else
             results%partBnk(i1)%coalesceNum = 0
             results%partBnk(i2)%coalesceNum = 0
          endif
       endif
    end do


! Done coalescing particles; shift reminaing nucleons down
    i = 1
130 if (results%partBnk(i)%coalesceNum == 0) then
       ! Particle did not coalesce; look at next particle
       if (i >= results%numParticles) go to 140 ! Reached end of particle list, stop checking
       i = i + 1
       go to 130
    endif
    ! Found a particle that was coalesced earlier, shift top particles down
    results%partBnk(i)%linearMomX   = results%partBnk(results%numParticles)%linearMomX
    results%partBnk(i)%linearMomY   = results%partBnk(results%numParticles)%linearMomY
    results%partBnk(i)%linearMomZ   = results%partBnk(results%numParticles)%linearMomZ
    results%partBnk(i)%coalesceFlag = results%partBnk(results%numParticles)%coalesceFlag
    results%partBnk(i)%kinEnergy    = results%partBnk(results%numParticles)%kinEnergy
    results%partBnk(i)%restMass     = results%partBnk(results%numParticles)%restMass
    results%partBnk(i)%charge       = results%partBnk(results%numParticles)%charge
    results%partBnk(i)%coalesceNum  = results%partBnk(results%numParticles)%coalesceNum
    results%partBnk(i)%strangeness  = results%partBnk(results%numParticles)%strangeness
    results%partBnk(i)%numBaryons   = results%partBnk(results%numParticles)%numBaryons
    results%numParticles = results%numParticles - 1
    if (i <= results%numParticles) go to 130
140 continue
150 continue


    ! All particles coalesced;
    ! 'Reset' portion of array that was emptied from coalesced particles
    do i = (results%numParticles + 1), results%numParticles
       results%partBnk(i)%linearMomX   = zro
       results%partBnk(i)%linearMomY   = zro
       results%partBnk(i)%linearMomZ   = zro
       results%partBnk(i)%coalesceFlag = zro
       results%partBnk(i)%kinEnergy    = zro
       results%partBnk(i)%restMass     = zro
       results%partBnk(i)%charge       = zro
       results%partBnk(i)%coalesceNum  = zro
       results%partBnk(i)%strangeness  = zro
       results%partBnk(i)%numBaryons   = zro
    end do


    ! Set number of "new" fragments (those created from coalescence)
    results%numFormedFragments = results%numCoalesced - &
         & ( results%numParticles - results%numParticles )


    return

! ======================================================================
1500 format("More ", A, " particle(s) were created due to ", &
          & "coalescence than can be stored in memory. Coalescence of ", &
          & "these particles will no longer be considered.")
2000 format("Divide by zero error prevented in ", &
          & "'coalesl.f90', line(s) ", A)
2100 format("Square root error prevented in ", &
          & "'coalesl.f90', line(s) ", A)
! ======================================================================

  end subroutine coalesl
