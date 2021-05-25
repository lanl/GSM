
  subroutine prinp (gsmObj, proj, targ, output, outData, icase)

! ======================================================================
!
!     Generates printout of a description of the current calculation
!     onto unit 14 (cem95.out) and unit 31 (cem95.res).
!     For an excitation function, only prints the description once,
!     but prints the parameters for each energy.
!
!    CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Modified by AJS, November, 1997.
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by K. K. Gudima, June, 2004.
!    Edited by A. J. Sierk,  January, 2005.
!    Modified by A. J. Sierk, December, 2005, following suggestion from
!      R. E. Prael.
!    Edited by AJS, LANL T-2, December, 2011 - January, 2012.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64, int64
    use molnixClass, only: molnixOptions
    use fissionBarrierClass, only: fissionBarrierOptions

    implicit none
    class(GSM),           intent(in   ) :: gsmObj
    class(GSMProjectile), intent(in   ) :: proj
    class(GSMTarget),     intent(in   ) :: targ
    class(GSMOutput),     intent(in   ) :: output
    class(OutputData),    intent(in   ) :: outData
    integer(int32),       intent(in   ) :: icase

    integer(int32) ::  icase0, ifirst, ih, j

! ======================================================================

    save ifirst, icase0

    data ifirst, icase0 /1, 0/

    type(molnixOptions) :: molOpt
    type(fissionBarrierOptions) :: fbOpt

! ======================================================================

    if (ifirst == 1) then
       write (31, 600)
       write (31, 650)
       write (31, 700)
       write (31, 800) trim(output%date)
    endif
    if (icase.ne.icase0) then
       icase0 = icase
       if (output%numComments > 0) then
          do ih = 1, output%numComments
             write (31, 900) trim(adjustl(output%comments(ih)))
             write (16, 900) trim(adjustl(output%comments(ih)))
          end do
       endif
       write (31, 2000) gsmObj%options%evapOpts%numEvapType
       write (31, 2050) gsmObj%options%preeqOpts%numPreeqType
       write (31, 1000) proj%restMass, proj%kinEnergy, &
            & targ%numBaryons, targ%numProtons, &
            & proj%numProtons, proj%numBaryons, &
            & output%numInelasticEvents, output%accountFission
       molOpt = gsmObj%genData%molnixEnergies%queryOptions()
       fbOpt  = gsmObj%genData%fissBarrier%queryOptions()
       write (31, 1200) proj%dKinEnergy, &
            & 1000.0_real64 * proj%kinEnergyMax, &
            & output%deltaTheta, &
            & output%energySpectra, &
            & output%multiplicities, &
            & output%channelCrossSection, &
            & output%nuclideCrossSection, &
            & output%doubleDiffSpectra,   &
            & output%angularSpectra, &
            & output%minEjectileRange, output%maxEjectileRange, &
            & fbOpt%r0m, molOpt%cevap
       if (output%doubleDiffSpectra >= 1) then
          write (31, 1300) (output%angleBins(j)%lowerBound, &
               & output%angleBins(j)%upperBound, j=1,10)
       endif
       if (output%energySpectra >= 1 .or. output%doubleDiffSpectra >= 1) then
          write (31, 1400) (output%energyBins(j)%lowerBound, output%energyBins(j)%upperBound, &
               & output%energyBinSubStep(j), j=1,4)
       endif
       write (31, 1500) output%maxEventAttempts
    endif
    write (31, 1100) outData%sigom
    if (output%accountFission == 0 .and. ifirst == 1) then
       write (31, 4300)
    elseif (output%accountFission > 0 .and. ifirst == 1) then
       write (31, 1600)
       write (31, 2100)
    endif

    if (ifirst == 1) then
       ifirst = 0
       write (31, 4400)
       write (31, 5800)
       write (31, 6000)
    endif
    return

! ======================================================================
600 format (/30x,'NOTICE', //3x,'This software and ancillary ', &
         &       'information (herein called "software") named', /3x,'GSM-1', &
         &       '.01 is made available under the terms described here. ', &
         &       'The software has', /3x,'been approved for release with ', &
         &       'associated LA-CC number LA-CC-04-085.', //8x, &
         &       'Copyright (2012). Los Alamos National Security, LLC.', &
         &       /3x,'This material was produced under U.S. Government ', &
         &       'contract', /3x,'DE-AC52-06NA25396 for Los Alamos ', &
         &       'National Laboratory, which is operated', /3x, &
         &       'by Los Alamos National Security, LLC, ', &
         &       'for the U.S. Department of Energy.', /3x, 'The Government ', &
         &       'is granted for itself and others acting on its behalf a ', &
         &      /3x,'paid-up, nonexclusive, irrevocable worldwide license ', &
         &       'in this material to', /3x,'reproduce, prepare derivative ', &
         &       'works, and perform publicly and display publicly.', //3x, &
         &       'NEITHER THE UNITED STATES NOR THE UNITED STATES ', &
         &       'DEPARTMENT OF ENERGY,', /3x,'NOR LOS ALAMOS NATIONAL ', &
         &       'SECURITY LLC, NOR ANY OF THEIR EMPLOYEES,', /3x,'MAKES ANY', &
         &       ' WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL ', &
         &       'LIABILITY')
650 format (3x,'OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, ', &
         &       'OR USEFULNESS OF ANY', /3x,'INFORMATION, APPARATUS, ', &
         &       'PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS', /3x, &
         &       'THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.')
700 format (/8x,'Additionally, the program is free software; you can', &
         &        ' redistribute', /3x,'it and/or modify it under the terms', &
         &        ' of the GNU General Public', /3x,'License as published by', &
         &        ' the Free Software Foundation; either', /3x,'version 2 of', &
         &        ' the License, or (at your option) any later version.', /3x, &
         &        'Accordingly, this program is distributed in the hope ', &
         &        'that it will', /3x,'be useful, but WITHOUT ANY WARRANTY; ', &
         &        'without even the implied', /3x,'warranty of MERCHANTABI', &
         &        'LITY or FITNESS FOR A PARTICULAR PURPOSE.', /3x,'See the ', &
         &        'GNU General Public License for more details.', //11x,'The', &
         &        ' primary authors of CEM03.03 are: S. G. Mashnik (LANL),', &
         &        /6x,'K. K. Gudima (IAP), and A. J. Sierk (LANL); with ', &
         &        'important contributions', /6x,'from R. E. Prael (LANL), ', &
         &        'M. I. Baznat (IAP), ', &
         &        'and N. V. Mokhov (Fermilab).', /6x,'(IAP = Institute of ', &
         &        'Applied Physics, Academy of Science of Moldova.)', //2x, &
         &        '-------------------------------------------------------', &
         &        '-----------------------', /)
800 format (1x,a24/)
900 format (a)
1000 format (/1x,'   M      T0     A    Z   Q  B   limc   idel', &
          &    /es9.3, ", ", f7.4, ", ", &
          &       f5.0, ", ", f5.0, ", ", &
          &       i3,   ", ", i3,   ", ", i8, 2x, i2/)
1100 format (/1x,'Geometrical cross section =',f8.2,' mb.')
1200 format (1x,'dt0 = ',f6.1,', t0max = ',f6.1,', dTheta = ', f4.1, &
          &       //'  mspec mpyld  mchy  misy mdubl ', &
          &       'mang  ipar1 ipar2 ' &
          &       /1x,8(2x,i2,2x)//1x, &
          &       'r0m =',f5.1,', & cevap =',f5.1,'.')
1300 format (/1x,'   Theta1      Theta2      Theta3      Theta4     ', &
          &       ' Theta5      Theta6'/1x,f5.1,11f6.1//1x, &
          &       '   Theta7      Theta8      Theta9     Theta10'/ &
          &         1x,f5.1,7f6.1/)
1400 format (/1x,'Tmin, Tmax, dT{1};Tmin, Tmax, dT{2};Tmin, Tmax, ', &
          &       'dT{3};Tmin, Tmax, dT{4}.'/2f6.1,f6.2,2f6.1,f6.2, &
          &        2f6.0,f6.1,3f6.0)
1500 format (/1x,'lim = ',i10,' .')
1600 format (1x,'Calculation takes into account fission and ', &
          &           'evaporation processes')
! 1800 format (/1x,'Preequilibrium emission is included; sigpre = ',f5.2)
! 1900 format (1x,'Preequilibrium emission is excluded.')
2000 format (1x,'Number of types of evaporated particles = ',i2)
2050 format (' Number of types of preequilibrium particles = ', i2)    !LMK 07/2012
2100 format (1x,'using Furihata-s GEM2 code.')
4300 format (1x,'Calculated without taking into account fission ', &
          &       'processes.'/)
4400 format(1x,'The following level density parameters were used ', &
          &      'in'/1x,'the preequilibrium part of this calculation:')
5800 format (1x,'a(Z,N,E) was calculated with Moller, Nix, Myers ', &
          &       '& Swiatecki microscopic'/1x,'corrections; ', &
          &       '[Atomic Data Nucl. Data Tables, 59, 185 (1995)];')
6000 format (1x,'Level density is from a shifted Fermi-gas formula, ', &
          &       'with the shift given by'/1x,' 0, delta-p, delta-n, or ', &
          &       'delta-p + delta-n, for odd-odd; odd-n, even-p;'/1x, &
          &       'odd-p, even-n; and even-even nuclei, respectively for ', &
          &       'the compound'/1x,'nucleus, and similarly using ', &
          &       'deltaM-p and deltaM-n for the saddle point.'/1x, &
          &       'delta-n and delta-p are tabulated by Moller, Nix & Kratz', &
          &       ' and'/1x,'deltaM-n and deltaM-p are 4.80 MeV * Bs * ', &
          &       '{1/N**(1/3) or 1/Z**(1/3)}.'/1x,'Bs is the surface area ', &
          &       'of the saddle-point shape with respect to a sphere.')

! ======================================================================
  end subroutine prinp
