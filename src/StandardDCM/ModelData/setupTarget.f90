
  function setupTarget (dataObj) result(errorFlag)

! ======================================================================
!
!     Neutron density, proton density, nucleon density .
!     Same functional form, different normalizations.
!     Determines zone boundaries in nucleus.
!     Finds average: Fermi momentum, Coulomb potential, nucleon density
!     in each zone of nucleus.
!
!   Rsm(ia) is the value of the radius at which the Woods-Saxon
!   nucleon density reaches the value of a(ia) * (the central value).
!   Hin is the integral from rsm(ia-1) (or 0) to rsm(ia) of
!   the density * r**2.
!   Sumhin is the integral of r**2/(1 + exp((r-rn)/bn) from 0 to
!   rsm(n)  [1/(4*pi) * (unnormalized) total mass inside of rsm(n)].
!   Rhon0 is the density normalization factor so all nucleons are inside
!   rsm(n).
!   By an integration by parts, the Coulomb potential at r is
!   proportional to:  1/r *Int[0,r]{r'**2 * rhop(r')}dr' +
!   Int[r,rsm(n)]{r' * rhop(r')}dr' .
!   Fis is proportional to this quantity * r**2.
!   Rhon is the average nucleon density between rsm(ia-1) and rsm(ia);
!   af is the average Coulomb potential (for a unit charge) in zone ia.
!   Tfp and tfn are the average Fermi momenta for protons and neutrons
!   in the nuclear zone ia.
!   Rbig is the normalized [to rsm(n)] radius of the outer boundary of
!   each zone.
!
! NOTE: Modified labeling according to:
!    (in options data type)
!       n    => dataObj%options%numZones
!       r0n  => dataObj%options%r0
!       bn   => dataObj%options%expDenom
!       rmax => dataObj%options%maxRad
!       a    => dataObj%options%aveDen
!    (in target data type)
!       rsm  => zoneBoundR
!       rbig => zoneBRFrac
!       rhop => protonDensity
!       rhon => neutronDensity
!       af   => coulombPote
!       tfp  => protFermiMom
!       tfn  => neutFermiMom
!
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk   LANL  T-2,  February-March, 1996.
!   "Last" change: 14-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMDataParams, only: zro, one, thr, five, twpi, fourpi, &
         & esq, ato3rd, twthrd

    implicit none
    class(StandardDCMData), intent(inout) :: dataObj
    integer(int32)                        :: errorFlag

    integer(int32) :: ia, k, outerBound1
    real(real64)   :: endn, ennuc, fac, fi0, pden, rhon0, rn, &
         & sumhin, temp, tempR
    real(real64), dimension(maxZonesAllowed) :: hin=zro, fi=zro

! ======================================================================

    procedure(IntegralInterface), pointer :: sfisPtr => sfis
    procedure(IntegralInterface), pointer :: ss2Ptr => ss2

! ======================================================================

    errorFlag = properlyConstructed

    ! Obtain parameters for the density equation
    ia = nint(dataObj%target%numBaryons)
    rn = dataObj%options%r0*ato3rd(ia)

    ! Obtain fractions of each nucleon type in the nucleus
    ennuc = dataObj%target%numBaryons - dataObj%target%numProtons   ! Num neutrons
    pden = dataObj%target%numProtons/dataObj%target%numBaryons      ! Proton  fraction in nucleus
    endn = ennuc/dataObj%target%numBaryons      ! Neutron fraction in nucleus


    ! Obtain radial boundaries of each zone [fm]
    ! NOTE: Protection for IEEE errors is provided for in the constructor.
    do ia = 1, dataObj%options%numZones
       dataObj%target%zoneBoundR(ia) = rn + dataObj%options%expDenom * &
            & log(one/dataObj%options%aveDen(ia) - one)
    end do
    ! Last zone bound extends to the max allowed
    dataObj%target%zoneBoundR(dataObj%options%numZones+1) = dataObj%options%maxRad


    ! Obtain integral values regarding radial bounds for neutron density (hin) and fermi energy (fi)
    do ia = 1,dataObj%options%numZones
       if (ia == 1) then
          call sfint (sfisPtr, dataObj%target%zoneBoundR(ia), zro, &
               & dataObj%options%expDenom, rn, fi(ia), dataObj)
          call sfint (ss2Ptr, dataObj%target%zoneBoundR(ia),  zro, &
               & dataObj%options%expDenom, rn, hin(ia), dataObj)
       else
          call sfint (sfisPtr, dataObj%target%zoneBoundR(ia), dataObj%target%zoneBoundR(ia-1), &
               & dataObj%options%expDenom, rn, fi(ia), dataObj)
          call sfint (ss2Ptr, dataObj%target%zoneBoundR(ia), dataObj%target%zoneBoundR(ia-1), &
               & dataObj%options%expDenom, rn, hin(ia), dataObj)
       endif
    end do

    ! Sum all portions of integral
    sumhin = zro
    do ia = 1,dataObj%options%numZones
       sumhin = sumhin + hin(ia)
    end do

    temp = fourpi*sumhin   ! Assumed isotropic (hence *4pi)
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(dataObj%io%message,1000) "108"
       call dataObj%io%print(4, 3, dataObj%io%message)
    end if
    ! Obtain base line neutron density
    rhon0 = thr*dataObj%target%numBaryons/(temp)
    temp = sumhin
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(dataObj%io%message,1000) "114"
       call dataObj%io%print(4, 3, dataObj%io%message)
    end if
    ! Obtain base line neutron fermi energy
    fi0 = thr*esq*dataObj%target%numProtons/temp

    ! For all zones, establish neutron densities and fermi energies
    tempR = dataObj%target%zoneBoundR(dataObj%options%numZones)
    if ( tempR < div0Lim .and. temp > -div0Lim ) then
       write(dataObj%io%message,1000) "137"
       call dataObj%io%print(4, 3, dataObj%io%message)
       tempR = div0Lim
    end if
    do ia = 1,dataObj%options%numZones

       ! Obtain a volume component [1 / (deltaR**3)]
       if (ia == 1) then
          temp = dataObj%target%zoneBoundR(ia)
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(dataObj%io%message,1000) "122"
             call dataObj%io%print(4, 3, dataObj%io%message)
          end if
          fac = one/temp**3
       else
          temp = dataObj%target%zoneBoundR(ia)**3 - dataObj%target%zoneBoundR(ia-1)**3
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(dataObj%io%message,1000) "129"
             call dataObj%io%print(4, 3, dataObj%io%message)
          end if
          fac = one/(temp)
       endif

       ! Obtain neutron density
       dataObj%target%neutronDensity(ia) = hin(ia)*rhon0*fac
       dataObj%target%coulombPote(ia) = fi0*fi(ia)*fac

       ! Obtain fermi energies based on the respective nucleon fraction
       dataObj%target%protFermiMom(ia) = 0.198595d0 * &
            & (dataObj%target%neutronDensity(ia)*pden)**twthrd
       dataObj%target%neutFermiMom(ia) = 0.198322d0 * &
            & (dataObj%target%neutronDensity(ia)*endn)**twthrd

       ! Obtain fraction of zone bound to that of the last zone bound (excluding max)
       dataObj%target%zoneBRFrac(ia) = dataObj%target%zoneBoundR(ia)/tempR

    end do
    ! Obtain fraction of maximum zone bound to that of the last zone bound (excluding max)
    dataObj%target%zoneBRFrac(dataObj%options%numZones+1) = dataObj%options%maxRad/tempR


    ! Set outermost boundary and remaining bins for nucleon densities and Fermi energies
    outerBound1 = dataObj%options%numZones + 1
    dataObj%target%coulombPote(outerBound1)    = &
         & dataObj%target%numProtons*esq/dataObj%options%maxRad
    dataObj%target%protFermiMom(outerBound1)  = zro
    dataObj%target%neutFermiMom(outerBound1)  = zro
    dataObj%target%neutronDensity(outerBound1) = zro
    dataObj%target%protonDensity(outerBound1)  = zro
    k = dataObj%options%numZones + 2
    do ia = k, maxZonesAllowed
       dataObj%target%zoneBoundR(ia)     = zro
       dataObj%target%zoneBRFrac(ia)     = zro
       dataObj%target%neutronDensity(ia) = zro
       dataObj%target%protonDensity(ia)  = zro
       dataObj%target%coulombPote(ia)    = zro
       dataObj%target%protFermiMom(ia)  = zro
       dataObj%target%neutFermiMom(ia)  = zro
    end do

    ! Scale neutron/proton zone densities according to the neutron/proton fraction
    do ia = 1, dataObj%options%numZones
       dataObj%target%protonDensity(ia) = dataObj%target%neutronDensity(ia)*pden
       dataObj%target%neutronDensity(ia) = dataObj%target%neutronDensity(ia)*endn
    end do

    ! Obtain geometric cross section:
    dataObj%target%geomCrossSection = twpi * five * &
         & ( dataObj%target%zoneBoundR( dataObj%options%numZones )**2 )

    return
! ======================================================================
1000 format("Divide by zero error prevented in 'setuptarget.f90', ", &
          & "line(s) ", A)
! ======================================================================
  end function setupTarget
