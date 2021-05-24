
  subroutine fisgem (evapObj, compound, fisinh, results)

! ======================================================================
!
!****************This routine was originally in LAHET     **************
!  FIS
!      Pick post-fission parameters such as mass, charge, kinetic energy
!      and excitation energy
!      Replaced inefficient GAUSSN routine with GAUSS2, AJS 10/14/03.
!
! <variables>
!     a    :   mass of nucleus before fission                     (IN)
!     z    :   charge  of nucleus before fission                  (IN)
!     u    :   excitation energy of nucleus before fission        (IN)
!  fisinh  : true=fission occur, false=no fission                 (IN)
!  zfis    :  charge of the fission fragment                     (OUT)
!  afis    :  mass of the fission fragment                       (OUT)
!  ufis    :  excitation energy of the fission fragment          (OUT)
!   er     :  recoil energy of the fission fragment              (OUT)
!  betf    :  recoil direction of the fission fragment           (OUT)
!
!   "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (Evap. class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, two, amu, twpi, ato3rd

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    type(evapCompound), intent(in   ) :: compound
    logical,            intent(inout) :: fisinh
    type(evaporationResults), intent(inout) :: results

    integer(int32) :: i, iaf, iafis1, iafis2, in, isym, izf, &
         & j, ja, jz, na, nck, nck2
    real(real64)   :: a1, a1twthr, a2twthr, amc1, amc2, amcf, &
         & amdiff, amean, arg, be0, be1, be2, bf, ef, gam, pcm, pcms, &
         & pcmx, pcmy, pcmz, ph, pres, proba, redm, rm0, rm1, rm2, &
         & sigmak, sigmas, sigz, spcm1i, spcm2i, svi, temp, temp2, &
         & th, totke, totkm, upr, v1sq, v1x, v1y, v1z, v2sq, v2x, &
         & v2y, v2z, vres, x, xx, z1, z2

    ! Fission fragment information
    real(real64), dimension(2)   :: zfis = zro   != Z of fission fragment
    real(real64), dimension(2)   :: afis = zro   != A of fission fragment
    real(real64), dimension(2)   :: ufis = zro   != (U or Ek) of fission fragment
    real(real64), dimension(2)   :: er   = zro   != E* of fission fragment
    real(real64), dimension(2,3) :: betf = zro   != (px/y/z or angle info) of fission fragment

! ======================================================================

!   How many times to try to find fission fragments
    integer(int32), parameter :: getFissFragmentTries = 20_int32

!   Following data for asymmetric vs symmetric choice:
    real(real64), parameter, dimension(4) :: evodba = &
         & [18.8d0, 18.1d0, 18.1d0, 18.5d0]

!   Asym gaussian width and high-mass mean:
    real(real64), parameter :: sigmaa = 6.5d0
    real(real64), parameter :: aamean = 140.d0

!   Mass excess for neutron, diff. of p & n masses and n mass:
    real(real64), parameter :: afact     = 8.071323d0  ! Neutron mass excess
    real(real64), parameter :: zfact     = 0.782354d0  ! Difference in n and p masses
    real(real64), parameter :: enmass    = 939.56563d0 ! neutron mass

! ======================================================================

    jz = nint(compound%numProtons)
    ja = nint(compound%numBaryons)
    if (.not.fisinh) then

!***********************************************************************
!
!  Logic error FISGEM called with fisinh unset:
!
!***********************************************************************

       write (evapObj%io%message, 60) ja, jz, compound%kinEnergy, compound%recEnergy
       call evapObj%io%print(1, 1, evapObj%io%message)
       return   ! No fission fragments have been created yet
    endif

!***********************************************************************
!
!   Pick the masses:
!
!***********************************************************************

    nck2 = 0
10  continue
    nck2 = nck2 + 1
    if (nck2 >  getFissFragmentTries) then

!***********************************************************************
!
!     Fission failure (invalid fission fragments found 'getFissFragmentTries' times):
!
!***********************************************************************

       write (evapObj%io%message, 50) ja, jz, compound%kinEnergy, compound%recEnergy, &
            & zfis(1), afis(1), zfis(2), afis(2)
       call evapObj%io%print(1, 1, evapObj%io%message)
       fisinh = .false.

       go to 100 ! Stores fission fragment information to the class and returns

    endif
    isym = 0
    temp = compound%numBaryons
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
        write(evapObj%io%message,1000) "112"
        call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    temp = compound%numProtons * compound%numProtons /temp
    if (temp <= 35.d0) then
       if (jz <= 88) go to 20

!***********************************************************************
!
!   High-Z fission mass distribution. Competition for sym vs. asym
!      is a simple symmetric to asymmetric data fit:
!
!***********************************************************************

    elseif (compound%kinEnergy <= 62.d0) then
       arg = -0.36d0*compound%kinEnergy
       arg = 4.87d3*exp(arg)
       temp = one + arg
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) "131"
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       proba = arg/(temp)
       if (evapObj%rang() <= proba) then

!***********************************************************************
!
!    Asymmetric fission:
!
!***********************************************************************

          a1 = evapObj%gauss2 (aamean, sigmaa)
          go to 30
       endif
    endif

!***********************************************************************
!
!     Find asymmetric barrier for width computation;
!     asymmetric barrier from Seaborg and Vandenbosch
!     Phys. Rev. 88, 507 (1952) Phys. Rev. 110, 507 (1958)
!     Find if ee, eo, oe, or oo nucleus:
!     na=1 odd-odd, 2 even-odd, 3 odd-even, 4 even-even
!
!***********************************************************************

    in = ja - jz
    na = 1
    if (jz == 2*(jz/2)) na = na + 1
    if (in == 2*(in/2)) na = na + 2
    temp = compound%numBaryons
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) "164"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    temp = compound%numProtons * compound%numProtons /temp
    ef = evodba(na) - 0.36d0*temp
20  continue
!  Furihata ... May 4, 1999:
    if (evapObj%options%fissParameter.ne.defaultFissParameter) then
       upr = compound%kinEnergy - ef
       upr = min (upr, 1.d+02)
       sigmas = 0.425d0*(upr*(one - 0.005d0*upr) + 9.35d0)
    else
       temp = compound%numBaryons
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) "178"
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       xx = compound%numProtons * compound%numProtons /temp
       bf = evapObj%efms (compound%numProtons, compound%numBaryons)
       if (bf < zro) bf = ef
       upr = min(compound%kinEnergy - bf, 450.0d0)
       sigmas = 0.122d0*xx**2 - 7.77d0*xx + 134.0d0  + 3.32d-2*upr
    endif

!***********************************************************************
!
!     Sigmas is symmetric fission mass width; Taken from systematics of
!     Neuzil & Fairhall, Phys. Rev. 129, 2705 (1963)
!
!      Low-Z fission is always symmetric
!
!      High-Z fission is sometimes sym. Loop back here from 1 loop if
!      symmetric fission predicted.
!
!***********************************************************************
!  Furihata ... May 4, 1999:

    amean = 0.5d0*compound%numBaryons
    a1 = evapObj%gauss2 (amean, sigmas)
    isym = 1
30  continue

!***********************************************************************
!
!   Loop for asymmetrin fission returns here:
!
!***********************************************************************

    afis(1) = anint(a1)

!***********************************************************************
!
!   Check for low final A:
!
!***********************************************************************
!    afis(1) = max(5.d0, afis(1))
!   Change to prevent fission fragments with A < 13: (AJS 8/09/08)

    afis(1) = max(13.d0, afis(1))
    afis(2) = compound%numBaryons - afis(1)
    if (afis(2) < 13.d0) then
       afis(1) = compound%numBaryons - 13.d0
       afis(2) = 13.d0
    endif
!***********************************************************************
!
!    Pick the charge:
!
!***********************************************************************

    iafis1 = nint(afis(1))
    iafis2 = nint(afis(2))
    a1twthr = ato3rd(iafis1)**2
    a2twthr = ato3rd(iafis2)**2
    temp = 131.d0 + a1twthr
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) "240"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    z1 = 65.5d0*afis(1)/(temp)
    temp = 131.d0 + a2twthr
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) "246"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    z2 = 65.5d0*afis(2)/(temp)
    z1 = z1 + 0.5d0*(compound%numProtons - z1 - z2)

!***********************************************************************
!
!   We use constant charge density with a 2 unit gaussian smearing:
!   (Changed to 0.5 by Furihata:)
!
!***********************************************************************
!   sigma = 0.75;   z1 = gaussn (z1, dph)
!   sigma = 0.75;   z1 = evapObj%gauss2 (z1, dph)     (AJS)
!                                by Furihata 18/DEC/1997
!***********************************************************************

    if (evapObj%options%fissParameter.ne.defaultFissParameter) then
       sigz = two
    else
       sigz = 0.75d0
    endif
    z1 = evapObj%gauss2 (z1, sigz)
    zfis(1) = anint(z1)
    zfis(2) = compound%numProtons - zfis(1)

!***********************************************************************
!
!   check for reasonable Z, A combinations.
!
!***********************************************************************

    if (zfis(1) >= afis(1)) go to 10
    if (zfis(2) >= afis(2)) go to 10
    if (zfis(1) < one) go to 10
    if (zfis(2) < one) go to 10

!***********************************************************************
!
!      Compute binding energy and actual masses of fragments
!
!***********************************************************************

    be0 = afact*compound%numBaryons  - zfact*compound%numProtons - evapObj%energy (jz, ja)
    rm0 = enmass*compound%numBaryons - zfact*compound%numProtons - be0
    iaf = nint(afis(1))
    izf = nint(zfis(1))
    be1 = afact*afis(1) - zfact*zfis(1) - evapObj%energy (izf, iaf)
    rm1 = enmass*afis(1) - zfact*zfis(1) - be1
    iaf = nint(afis(2))
    izf = nint(zfis(2))
    be2 = afact*afis(2) - zfact*zfis(2) - evapObj%energy (izf, iaf)
    rm2 = enmass*afis(2) - zfact*zfis(2) - be2

!***********************************************************************
!
!   Pick recoil kinetic energy; use the systematics of 
!   Unik et al. Proc. 3rd IAEA Symp. on Phys. & Chem. of Fission,
!   Rochester, 1973 Vol. II, pg. 19.
!
!***********************************************************************

    temp = ato3rd(ja)
    if (evapObj%options%fissParameter.ne.defaultFissParameter) then
       totkm = 0.13323d0*compound%numProtons*compound%numProtons/temp - 11.4d0
    else
!   4 May 1999 by Furihata
       x = compound%numProtons*compound%numProtons/temp
       if (x <= 900.d0) then
          totkm = 0.131d0*x
       else if (x <= 1800.d0) then
          totkm = 0.104d0*x + 24.3d0
       else
          write(evapObj%io%message,2000)
          call evapObj%io%print(1, 2, evapObj%io%message)
       endif
    endif

!***********************************************************************
!
!   Use a width of 15% value at half height (RAL):
!
!***********************************************************************

    if (evapObj%options%fissParameter.ne.defaultFissParameter) then
       sigmak = 0.084d0*totkm
!   4 May 1999 by Furihata
    else
       if (x < 1000.d0) then
          sigmak = 86.5d0
       else if (x < 1800.d0) then
          sigmak = 5.70d-4*(x - 1000.d0)**2 + 86.5d0
       else
          sigmak = 0
          write(evapObj%io%message,2100)
          call evapObj%io%print(1, 2, evapObj%io%message)
       endif
       temp = sigmak
       if (temp < 0.0d0) then
          temp = 0.01d0
          write(evapObj%io%message,1100) "346"
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       sigmak = sqrt(temp)
    endif

!***********************************************************************
!
!   Check that event is energetically possible:
!   (Cut off gaussian KE distribution at maximum allowable energy.)
!
!***********************************************************************

    temp2 = compound%kinEnergy + be1 + be2 - be0
    nck = 0
40  continue
    totke = evapObj%gauss2 (totkm, sigmak)
    if (nck > 10) go to 10
    nck = nck + 1
    if (totke > temp2) go to 40
    
!***********************************************************************
!
!   Pick excitation from equidistribution of original plus 
!   energy balance.
!
!***********************************************************************

    temp = compound%numBaryons
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) "376-378"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    temp = (temp2 - totke)/temp
    ufis(1) = afis(1)*temp
    ufis(2) = afis(2)*temp
!   Find total masses, including excitation energies, at evap time:
    amcf = rm0 + compound%kinEnergy
    amc1 = rm1 + ufis(1)
    amc2 = rm2 + ufis(2)
    amdiff = amcf - amc1 - amc2
    if (amdiff < zro) go to 40

!***********************************************************************
!
!     amdiff = ekin should be satisfied
!
!***********************************************************************

!   Velocity of pre-fission nucleus : vres*bet0:
!    pres = sqrt(erec*2 + two*a*amu*erec) error found by KKG 10 Feb 09.
    temp = compound%recEnergy**2 + two*compound%numBaryons*amu*compound%recEnergy
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(evapObj%io%message,1100) "399"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    pres = sqrt(temp)
    temp = compound%recEnergy + compound%numBaryons*amu
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) "405"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    vres = pres/(temp)
    
!   Momentum of CM system:
    temp = compound%numBaryons
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) "413"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    redm = amu*(afis(1)*afis(2))/temp
    pcms = totke**2 + two*redm*totke
    temp = pcms
    if (temp < 0.0d0) then
       temp = 0.01d0
       write(evapObj%io%message,1100) "420"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    pcm = sqrt(temp)
    th  = acos(two*evapObj%rang() - one)
    ph  = twpi*evapObj%rang()
    pcmx = pcm*sin(th)*cos(ph)
    pcmy = pcm*sin(th)*sin(ph)
    pcmz = pcm*cos(th)
    temp = pcms + (amu*afis(1))**2
    if (temp < div0Lim) then
       temp = 0.01d0
       write(evapObj%io%message,1200) "431"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    spcm1i = one/sqrt(temp)
    temp = pcms + (amu*afis(2))**2
    if (temp < div0Lim) then
       temp = 0.01d0 
       write(evapObj%io%message,1200) "437"
       call evapObj%io%print(4, 3, evapObj%io%message)
   end if
    spcm2i = one/sqrt(temp)

!   Velocity of fission fragment 1 in CM system:
    v1x = spcm1i*pcmx
    v1y = spcm1i*pcmy
    v1z = spcm1i*pcmz
!   Velocity of fission fragment 2 in CM system:
    v2x = -spcm2i*pcmx
    v2y = -spcm2i*pcmy
    v2z = -spcm2i*pcmz

!   Boost to Lab system:
    v1x = v1x + vres*compound%linearMomFrac(1)
    v1y = v1y + vres*compound%linearMomFrac(2)
    v1z = v1z + vres*compound%linearMomFrac(3)
    v2x = v2x + vres*compound%linearMomFrac(1)
    v2y = v2y + vres*compound%linearMomFrac(2)
    v2z = v2z + vres*compound%linearMomFrac(3)

!   Kinetic Energy and velocity of fission fragment 1 in Lab system:
    v1sq = v1x**2 + v1y**2 + v1z**2
    temp = one - v1sq
    if (temp < div0Lim) then
       temp = 0.01d0
       write(evapObj%io%message,1200) "463"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    gam = one/sqrt(temp)
    er(1) = afis(1)*amu*(gam - one)
    temp = v1sq
    if (temp < div0Lim) then
       temp = 0.01d0
       write(evapObj%io%message,1200) "470"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    svi = one/sqrt(temp)
    betf(1,1) = svi*v1x
    betf(1,2) = svi*v1y
    betf(1,3) = svi*v1z

!   Kinetic Energy and velocity of fission fragment 2 in Lab system:
    v2sq = v2x**2 + v2y**2 + v2z**2
    temp = one - v2sq
    if (temp < div0Lim) then
       temp = 0.01d0
       write(evapObj%io%message,1200) "482"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    gam = one/sqrt(temp)
    er(2) = afis(2)*amu*(gam - one)
    temp = v2sq
    if (temp < div0Lim) then
       temp = 0.01d0
       write(evapObj%io%message,1200) "489"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    svi = one/sqrt(temp)
    betf(2,1) = svi*v2x
    betf(2,2) = svi*v2y
    betf(2,3) = svi*v2z


100 continue

    ! Get fission fragments (fission fragments successfully created)
    do i = 1, 2
       results%numFissionFragments = results%numFissionFragments + 1
       results%fissionBnk(results%numFissionFragments)%numBaryons  = afis(i)
       results%fissionBnk(results%numFissionFragments)%numProtons  = zfis(i)
       results%fissionBnk(results%numFissionFragments)%kinEnergy   = ufis(i)
       results%fissionBnk(results%numFissionFragments)%rotEnergy   = er(i)
       do j = 1, 3
          results%fissionBnk(results%numFissionFragments)%linearMomFrac(j) = betf(i,j)
       end do
       ! Flag that fission fragment was created
       results%fissionBnk(results%numFissionFragments)%physicsFlag = evaporationFlag
    end do

    return

! ======================================================================
50  format ('---> Fission failed: ja = ',i5,'  jz = ',i5/'-     u = ', &
         e10.3,'      erec = ',e10.3/'     zfis1 = ',e10.3,'     afis1 = ', &
         e10.3/'     zfis2 = ',e10.3, '     afis2 = ',e10.3)
60  format (//'  Logic error in fission. FISGEM called with fisinh ', &
         'flag unset.',2i10,2f10.5)
! ----------------------------------------------------------------------
1000 format("Divide by zero error prevented in 'fisgem.f90', line ", A)
1100 format("Square root error prevented in 'fisgem.f90', line ", A)
1200 format("Divide by zero/square root error prevented ", &
          & "in 'fisgem.f90', line ", A)
2000 format("Error regarding the total kinetic energy in subroutine 'fisgem'.")
2100 format("(Z^2)/(A^(1/3)) too large in subroutine ", &
          & "'fisgem'.")
! ======================================================================
  end subroutine fisgem
