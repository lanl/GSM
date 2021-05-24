
  subroutine pointe (sDCM, clientTarg, partin, ipatin, partne, ipatne, nabs, nout, &
       & v, u, tin1, sigp, sign, sigabs, t0gev)

! ======================================================================
!
!     Determination of interaction point inside nucleus.
!
!     Called by: CASCAD
!
!     Calls: ABSORP GEOM8 PARTN REFRAC SIGMAT8 SLQEK TINVU
!
!    CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk,  LANL  T-2  February-March, 1996.
!    Edited by AJS, July-August, 1997
!    Modified by AJS, March, 1999.
!   "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================
!
!  nout = 1 when cascade particle is exiting the nucleus.
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
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, thr
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    class(StandardDCMData), intent(inout) :: clientTarg
    real(real64),           intent(inout) :: partin(9)
    integer(int32),         intent(inout) :: ipatin(5)
    real(real64),           intent(inout) :: partne(9)
    integer(int32),         intent(inout) :: ipatne(5)
    integer(int32),         intent(  out) :: nabs
    integer(int32),         intent(  out) :: nout
    real(real64),           intent(  out) :: v(3)
    real(real64),           intent(  out) :: u
    real(real64),           intent(  out) :: tin1
    real(real64),           intent(  out) :: sigp
    real(real64),           intent(  out) :: sign
    real(real64),           intent(  out) :: sigabs
    real(real64),           intent(in   ) :: t0gev

    integer(int32) :: ig, itemp, ksin, ksip, l, mb, me, ms
    real(real64)   :: deltsi, delx, piks, plambi, s, sizen, sk, temp, &
         & temp1, temp2, temp3, temp5

! ======================================================================

    sk = thr
    temp1 = log(sDCM%rang())
    ig = 0
!   s is the distance along the current trajectory to the first
!   intersection of the trajectory with a boundary of the current
!   nuclear zone, ipatin(5).
10  s = sDCM%geom8 (clientTarg, partin, ipatin)
!   Rerun this cascade if very rare error occurs in GEOM8.
    if (ipatin(1) == 2) return
20  ig = ig + 1
    if (ig > 50 .and. ipatin(4) > 0) then
!   Absorb nucleon if it has been here for too long.
       nabs = 1
       return
    endif
!  Randomly select a partner for interaction by isospin and Fermi motion
    call sDCM%partn (clientTarg, partin, ipatin, partne, ipatne)
!  Find tin1, the kinetic energy of the cascade particle in the frame
!  in which the partner is at rest. v is the velocity of the center of
!  mass of the 2 interacting particles with respect to the nuclear
!  rest (lab) frame. U is the total energy in the 2-particle C-M frame.
    call sDCM%tinvu (partin(8), partin(9), partne(8), partne(9), partin(4), &
         & partin(5), partin(6), partin(7), partne(4), partne(5), &
         & partne(6), partne(7), v, u, tin1)
!  Determine the index ksin for cascade particle interacting with a
!  neutron.
    call slqek (l, ms, mb, ksin, me, ipatin(2), ipatin(3), ipatin(4), &
         & ipatin(1), ipatne(2), ipatne(3), ipatne(4), 0)
!  Find the total cross section of the cascade particle (pi+, pi0, pi-,
!  p, n, or gamma) with a neutron in the nucleus.
    sign = sDCM%sigmat8 (l, ms, mb, ksin, 0, tin1)
!  Determine the index ksip for cascade particle interacting with a
!  proton.
    call slqek (l, ms, mb, ksip, me, ipatin(2), ipatin(3), ipatin(4), &
         & ipatin(1), ipatne(2), ipatne(3), ipatne(4), 1)

!  Find the total cross section of the cascade particle (pi+, pi0, pi-,
!  p, n, or gamma) with a proton in the nucleus.
    sigp = sDCM%sigmat8 (l, ms, mb, ksip, 0, tin1)

!  Find absorption cross section of cascade particle on two nucleons.
!  Only nonzero for pions and gammas!
    sigabs = sDCM%sigmat8 (l, ms, mb, ksip, 3, partin(8))
    itemp = ipatin(5)
!  Particle mean free path (in units of clientTarg%zoneBoundR(n); sigma in microbarns;
!  rho in fm^-3).
    temp = (clientTarg%neutronDensity(itemp) * sign + &
         & clientTarg%protonDensity(itemp)*sigp + &
         & clientTarg%protonDensity(itemp)*sigabs) * &
         & clientTarg%zoneBoundR( clientTarg%numZones() )
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "114"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    plambi = 10.d0/(temp)
!  -plambi*temp1 is the random distance [in units of clientTarg%zoneBoundR(n)] that
!  the particle will traverse before its next interaction.
    delx = -plambi*temp1
!  FCG; photonuclear extension:
    sizen = clientTarg%neutronDensity(1) / &
         & clientTarg%neutronDensity( clientTarg%numZones() )
    if (ipatin(2).ne.0) sizen = 1.0d10
    if (delx > 6.d0*sizen) then
!  The interaction distance is much greater than the nuclear size:
       partin(8) = t0gev
       nout = 1
    else
       temp2 = sk
       if (temp2 < div0Lim .and. temp2 > -div0Lim) then
          temp2 = div0Lim
          write(sDCM%io%message,1000) "135"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       deltsi = plambi/temp2
       temp5 = min(s, deltsi)
       if (plambi < div0Lim .and. plambi > -div0Lim) then
          plambi = div0Lim
          write(sDCM%io%message,1000) "141"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       piks = temp1 + temp5/plambi
       if (piks > zro) then
!   Interaction point is before intersection with zone boundary or
!   lambda/3.  Advance particle along trajectory; interaction will
!   occur at this point; still in the same zone. Then return to CASCAD.
          partin(1) = partin(1) + delx*partin(4)*partin(7)
          partin(2) = partin(2) + delx*partin(4)*partin(6)
          partin(3) = partin(3) + delx*partin(5)
          partne(1) = partin(1)
          partne(2) = partin(2)
          partne(3) = partin(3)
       else
!   Projected interaction point is either beyond the boundary of
!   the nuclear zone, or more than lambda/3.
!   Magnitude of temp1 is reduced; as is s.
          temp1 = piks
          s = s - temp5
!   Advance particle along trajectory either to intersection with
!   zone boundary, or a distance of lambda/3.
          partin(1) = partin(1) + temp5*partin(4)*partin(7)
          partin(2) = partin(2) + temp5*partin(4)*partin(6)
          partin(3) = partin(3) + temp5*partin(5)
!   If boundary not reached, go back to 20 with new value of temp1;
!   distance to the boundary  i  is reduced by the amount moved.
          if (s.ne.zro) go to 20
!   Particle moves to boundary of zone:
          call sDCM%refrac (clientTarg, partin, ipatin, nabs, nout)
          if (partin(8) <= 1.d-5) then
!   Absorb very low energy particle:
!   FCG photonuclear extension:
             if (ipatin(4) == 1 .or. ipatin(2).ne.0) then
                nabs = 1
                return
             endif
          endif
          if (ipatin(5) < clientTarg%numZones()+1) then
!   Particle is not outside nucleus;
!   Continue trajectory in new zone unless particle is absorbed:
             if (nabs <= 0) go to 10
          else
!   Particle is in zone n+1; will exit nucleus at clientTarg%zoneBoundR(n+1) = rmax
             temp3 = sDCM%geom8 (clientTarg, partin, ipatin)
!   Flag for rare problem in GEOM8.  AJS 3/99
             if (ipatin(1) == 2) return
             partin(1) = partin(1) + temp3*partin(4)*partin(7)
             partin(2) = partin(2) + temp3*partin(4)*partin(6)
             partin(3) = partin(3) + temp3*partin(5)
             call sDCM%refrac (clientTarg, partin, ipatin, nabs, nout)
             nout = 1
          endif
       endif
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'pointe.f90' line(s) ", A)
! ======================================================================
  end subroutine pointe


  subroutine pointe1 ( sDCM, clientTarg, partin, ipatin, partne, &
       & ipatne, nabs, nout, irefl, v, u, tin1, sigp, sign, sigabs, t0gev)

! ======================================================================
!
!     Determination of interaction point inside nucleus.
!     This version is used when refraction is enables.
!
!     Called by: CASCAD
!
!     Calls: ABSORP GEOM8 PARTN REFRAC1 SIGMAT8 SLQEK TINVU
!
!    CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk,  LANL  T-2  February-March, 1996.
!    Edited by AJS, July-August, 1997
!    Modified by AJS, March, 1999.
!   "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by K. K. Gudima,  Feb., 2004. (refrac==>refrac1)
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================
!
!  nout = 1 when cascade particle is exiting the nucleus.
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
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, thr
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    class(StandardDCMData), intent(inout) :: clientTarg
    real(real64),           intent(inout) :: partin(9)
    integer(int32),         intent(inout) :: ipatin(5)
    real(real64),           intent(inout) :: partne(9)
    integer(int32),         intent(inout) :: ipatne(5)
    integer(int32),         intent(  out) :: nabs
    integer(int32),         intent(  out) :: nout
    integer(int32),         intent(  out) :: irefl
    real(real64),           intent(  out) :: v(3)
    real(real64),           intent(  out) :: u
    real(real64),           intent(  out) :: tin1
    real(real64),           intent(  out) :: sigp
    real(real64),           intent(  out) :: sign
    real(real64),           intent(  out) :: sigabs
    real(real64),           intent(in   ) :: t0gev

    integer(int32) :: ig, itemp, ksin, ksip, l, mb, me, ms
    real(real64)   :: deltsi, delx, piks, plambi, s, sizen, sk, temp, &
         & temp1, temp2, temp3, temp5

! ======================================================================

    sk = thr
    temp1 = log(sDCM%rang())
    ig = 0

!   s is the distance along the current trajectory to the first
!   intersection of the trajectory with a boundary of the current
!   nuclear zone, ipatin(5).

10  s = sDCM%geom8 (clientTarg, partin, ipatin)
!   Rerun this cascade if very rare error occurs in GEOM8.
    if (ipatin(1) == 2) return
20  continue
    ig = ig + 1
    if (ig > 50 .and. ipatin(4) > 0) then
!   Absorb nucleon if it has been here for too long.
       nabs = 1
       return
    endif

!  Randomly select a partner for interaction by isospin and Fermi motion
    call sDCM%partn (clientTarg, partin, ipatin, partne, ipatne)

!  Find tin1, the kinetic energy of the cascade particle in the frame
!  in which the partner is at rest. v is the velocity of the center of
!  mass of the 2 interacting particles with respect to the nuclear
!  rest (lab) frame. U is the total energy in the 2-particle C-M frame.

    call sDCM%tinvu (partin(8), partin(9), partne(8), partne(9), partin(4), &
         & partin(5), partin(6), partin(7), partne(4), partne(5), &
         & partne(6), partne(7), v, u, tin1)

!  Determine the index ksin for cascade particle interacting with a
!  neutron.
    call slqek (l, ms, mb, ksin, me, ipatin(2), ipatin(3), ipatin(4), &
         & ipatin(1), ipatne(2), ipatne(3), ipatne(4), 0)

!  Find the total cross section of the cascade particle (pi+, pi0, pi-,
!  p, or n) with a neutron in the nucleus.
    sign = sDCM%sigmat8 (l, ms, mb, ksin, 0, tin1)

!  Determine the index ksip for cascade particle interacting with a
!  proton.
    call slqek (l, ms, mb, ksip, me, ipatin(2), ipatin(3), ipatin(4), &
         & ipatin(1), ipatne(2), ipatne(3), ipatne(4), 1)

!  Find the total cross section of the cascade particle (pi+, pi0, pi-,
!  p, or n) with a proton in the nucleus.
    sigp = sDCM%sigmat8 (l, ms, mb, ksip, 0, tin1)

!  Find absorption cross section of cascade particle on two nucleons.
!  Only nonzero for pions!
    sigabs = sDCM%sigmat8 (l, ms, mb, ksip, 3, partin(8))

    itemp = ipatin(5)
!  Particle mean free path (in units of clientTarg%zoneBoundR(n); sigma in microbarns;
!  rho in fm^-3).
    temp = (clientTarg%neutronDensity(itemp) * sign + &
         & clientTarg%protonDensity(itemp)*sigp + &
         & clientTarg%protonDensity(itemp)*sigabs) * &
         & clientTarg%zoneBoundR( clientTarg%numZones() )
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "320"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    plambi = 10.d0/(temp)

!  -plambi*temp1 is the random distance [in units of clientTarg%zoneBoundR(n)] that
!  the particle will traverse before its next interaction.
    delx = -plambi*temp1

!  FCG; photonuclear extension:
    sizen = clientTarg%neutronDensity(1) / &
         & clientTarg%neutronDensity( clientTarg%numZones() )
    if (ipatin(2).ne.0) sizen = 1.0d10
    if (delx > 6.d0*sizen) then
!  The interaction distance is much greater than the nuclear size:
       partin(8) = t0gev
       nout = 1
    else
       temp2 = sk
       if (temp2 < div0Lim .and. temp2 > -div0Lim) then
          temp2 = div0Lim
          write(sDCM%io%message,1000) "343"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       deltsi = plambi/temp2
       temp5 = min(s, deltsi)
       if (plambi < div0Lim .and. plambi > -div0Lim) then
          plambi = div0Lim
          write(sDCM%io%message,1000) "349"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       piks = temp1 + temp5/plambi
       if (piks > zro) then

!   Interaction point is before intersection with zone boundary or
!   lambda/3.  Advance particle along trajectory; interaction will
!   occur at this point; still in the same zone. Then return to CASCAD.
          partin(1) = partin(1) + delx*partin(4)*partin(7)
          partin(2) = partin(2) + delx*partin(4)*partin(6)
          partin(3) = partin(3) + delx*partin(5)
          partne(1) = partin(1)
          partne(2) = partin(2)
          partne(3) = partin(3)
       else
!   Projected interaction point is either beyond the boundary of
!   the nuclear zone, or more than lambda/3.
!   Magnitude of temp1 is reduced; as is s.
          temp1 = piks
          s = s - temp5

!   Advance particle along trajectory either to intersection with
!   zone boundary, or a distance of lambda/3.
          partin(1) = partin(1) + temp5*partin(4)*partin(7)
          partin(2) = partin(2) + temp5*partin(4)*partin(6)
          partin(3) = partin(3) + temp5*partin(5)

!   If boundary not reached, go back to 20 with new value of temp1;
!   distance to the boundary  i  is reduced by the amount moved.
          if (s.ne.zro) go to 20
!   Particle moves to boundary of zone:
          call sDCM%refrac1 (clientTarg, partin, ipatin, nabs, nout, irefl)
          if (nabs == 1) return

          if (partin(8) <= 1.d-5) then
!  FCG photonuclear extension:
             if (ipatin(4) == 1 .or. ipatin(2).ne.0) then
                nabs = 1
                return
             endif
          endif
          if (ipatin(5) <  clientTarg%numZones() + 1) then
!   Particle is not outside nucleus;
!   Continue trajectory in new zone unless particle is absorbed:
             if (nabs <= 0) go to 10
          else
!   Particle is in zone n+1; will exit nucleus at clientTarg%zoneBoundR(n+1) = rmax
             temp3 = sDCM%geom8 (clientTarg, partin, ipatin)
!   Flag for rare problem in GEOM8.  AJS 3/99
             if (ipatin(1) == 2) return
             partin(1) = partin(1) + temp3*partin(4)*partin(7)
             partin(2) = partin(2) + temp3*partin(4)*partin(6)
             partin(3) = partin(3) + temp3*partin(5)
             call sDCM%refrac1 (clientTarg, partin, ipatin, nabs, nout, irefl)
             if (nabs == 1) return
             nout = 1
          endif
       endif
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'pointe.f90' line(s) ", A)
! ======================================================================
  end subroutine pointe1
