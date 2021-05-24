
  subroutine elex (sDCM, v, u, tin1, partin, ipatin, partne, ipatne, &
       & mv, np, l, mb, ksi, me, sigex, sigelex, photoData, results)

! ======================================================================
!
!     Calculation of particle characteristics in elastic and
!     single charge exchange scattering.
!     Corrected for photon interactions by AJS, 10/21/03.
!       "elastic" means pion produced with no charge change of nucleon;
!       "exchange" means pion produced with charge change of nucleon.
!     Modified version to use Duarte's angular distribution approx.,
!     03/2004. This is for N + N at energies less than or equal to 
!     2 GeV.
!
!   Called by: TYPINT
!
!   Calls: ABEL CDUARTE COSEL COSEX
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk  LANL  T-2  February, 1996.
!    Edited by A. J. Sierk  LANL  T-2  December, 1997.
!    Modified by A. J. Sierk  LANL  T-2  February, 1999.
!    Modified by SGM on 12/05/01 on KKG's suggestion
!   "Last" change: 12-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by K. K. Gudima, Feb., 2004.
!    Edited by A. J. Sierk  LANL  T-16  March, 2004.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================
!
!  np is the number of particles in the final state (always = 2, unless
!  storage block /memory/ is overfilled).
!  mv is the index of the last storage cell in /memory/ filled prior to
!  the call to ELEX.
!
!  Definition of partin (partne similar for partner particle):
!                       partin(1); x coordinate of particle
!                       partin(2); y coordinate of particle
!                       partin(3); z coordinate of particle
!                       partin(4); sin(theta), direction of momentum
!                       partin(5); cos(theta), direction of momentum
!                       partin(6); sin(phi), direction of momentum
!                       partin(7); cos(phi), direction of momentum
!                       partin(8); kinetic energy of particle
!                       partin(9); rest mass of particle
!
!  Definition of ipatin (ipatne similar for partner particle):
!                       ipatin(1); charge of particle
!                       ipatin(2); non-zero for photon interactions
!                       ipatin(3); strangeness of particle
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two, twpi, massPiPM, massPi0, &
         & emneut, emprot

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    real(real64),   intent(in   ) :: v(3)
    real(real64),   intent(in   ) :: u
    real(real64),   intent(in   ) :: tin1
    real(real64),   intent(in   ) :: partin(9)
    integer(int32), intent(in   ) :: ipatin(5)
    real(real64),   intent(in   ) :: partne(9)
    integer(int32), intent(in   ) :: ipatne(5)
    integer(int32), intent(in   ) :: mv
    integer(int32), intent(  out) :: np
    integer(int32), intent(in   ) :: l
!    integer(int32), intent(in   ) :: ms
    integer(int32), intent(in   ) :: mb
    integer(int32), intent(in   ) :: ksi
    integer(int32), intent(in   ) :: me
    real(real64),   intent(in   ) :: sigex
    real(real64),   intent(in   ) :: sigelex
    type(sDCMPhotonCrossSections), intent(in   ) :: photoData
    class(StandardDCMResults), intent(inout) :: results

    integer(int32) :: ie, ne
    real(real64)   :: abc, b1, b2, betaex, cmi, cmn, ctsti, fisti, r1, temp
    real(real64), dimension(3) :: pist=zro, pnst=zro

! ======================================================================

    if (ipatin(2).ne.0) then
!  "elastic" for photons means produce pion with nucleon unchanged:
       cmi = massPi0
    else
       cmi = partin(9)
    endif
    cmn = partne(9)
!   Ratio of charge exchange to (elastic + SCX)
!   (Will be zero for N-N channel; non-zero for pi+ + p; pi- + n
!    and gamma + N):
    temp = sigelex
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "95"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    betaex = sigex/temp
    b1 = sDCM%rang()
    if (b1 >= betaex) then
!  Elastic scattering:
       ie = ipatin(1)
       ne = ipatne(1)
!  KKG 02/03/04:
!       ctsti = cosel (l, mb, ksi, tin1, partin(9))
       if (mb == 2 .and. tin1 <= two)  then
          ctsti = sDCM%cduarte (tin1, ipatin(1), ipatne(1), results)
       else
          ctsti = sDCM%cosel (l, mb, ksi, tin1, partin(9), photoData)
       endif
    else
!  Charge exchange scattering:
       if (ipatin(2).ne.0) then
!  Photon:
          if (ipatne(1).ne.0) then
!  gamma + p --> n + pi+
             ie =  1
             ne =  0
          else
!  gamma + n --> p + pi-
             ie = -1
             ne =  1
          endif
       else
          if (ipatin(1).ne.0) then
!   Charged pion --> neutral pion:
             ie = 0
!   Charge on outgoing nucleon = total charge in system:
             ne = me
          else
!   Neutral pion --> charged pion:
!   pi0 p --> pi+ n; or pi0 n --> pi- p:
             ne = 1 - ipatne(1)
             ie = me - ne
          endif
       endif
       abc = abs(dble(ie))
       cmi = massPi0*(one - abc) + massPiPM*abc
       cmn = emneut*(one - dble(ne)) + emprot*dble(ne)
       r1 = sDCM%rang()
       ctsti = sDCM%cosex (l, tin1, partin(9), r1, photoData)
    endif
!   Both options continue from here:
    b2 = sDCM%rang()
    fisti = twpi*b2
    call sDCM%abel (partin, v, u, pist, pnst, ctsti, fisti, cmi, cmn)
    if ( mv > results%maxProgenyM3 ) then
       np = 0
!       write (16, 2000)
!       write ( *, 2000)
    else
       results%pmemo(1,mv+3) = pist(1)
       results%pmemo(2,mv+3) = pist(2)
       results%pmemo(3,mv+3) = pist(3)
       results%pmemo(9,mv+3) = cmi
       results%imemo(1,mv+3) = ie
       results%imemo(2,mv+3) = 0
       results%imemo(3,mv+3) = ipatin(3)
       results%imemo(4,mv+3) = ipatin(4)
       results%pmemo(1,mv+1) = pnst(1)
       results%pmemo(2,mv+1) = pnst(2)
       results%pmemo(3,mv+1) = pnst(3)
       results%pmemo(9,mv+1) = cmn
       results%imemo(1,mv+1) = ne
       results%imemo(2,mv+1) = 0
       results%imemo(3,mv+1) = 0
       results%imemo(4,mv+1) = 1
       np = 2
    endif

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'elex.f90' line(s) ", A)
! 2000 format("Interacting bank filled in 'elex.f90'. Progeny will not be stored.")
! ======================================================================
  end subroutine elex
