
  subroutine pinpn (sDCM, clientProj, clientTarg, partin, ipatin, &
       & nout, am0, p0, t3)

! ======================================================================
!
!     Calculation of entry point of particle into nucleus, initial
!     angular momentum of nucleus.
!
!     Called by: CASCAD
!
!     Calls: POTEN
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Modified by A. J. Sierk,  LANL  T-2  March, 1996.
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!
! ======================================================================
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
!                       ipatin(1); charge of particle
!                       ipatin(2); non-zero for photon reactions
!                       ipatin(3); strangeness of particle
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two, twpi
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    type(sDCMProjectile),   intent(inout) :: clientProj
    class(StandardDCMData), intent(inout) :: clientTarg
    real(real64),           intent(  out) :: partin(9)
    integer(int32),         intent(  out) :: ipatin(5)
    integer(int32),         intent(  out) :: nout
    real(real64),           intent(  out) :: am0(3)
    real(real64),           intent(  out) :: p0
    real(real64),           intent(  out) :: t3

    real(real64) :: temp1, temp2, temp3

! ======================================================================

!   sin(theta) direction of momentum
    partin(4) = zro
!   cos(theta) direction of momentum
    partin(5) = one
!   cos(phi) direction of momentum
    partin(6) = zro
!   sin(phi) direction of momentum
    partin(7) = one
    partin(9) = clientProj%restMass
    ipatin(1) = clientProj%numProtons
!  clientProj%gammaFlag is nonzero for photons:
    ipatin(2) = clientProj%gammaFlag
    ipatin(3) = clientProj%decayNumber
    ipatin(4) = clientProj%numBaryons
    ipatin(5) = clientTarg%numZones()
    nout = 0
    temp1 = sDCM%rang()
    temp2 = twpi*sDCM%rang()
    temp3 = clientTarg%zoneBRFrac( clientTarg%numZones() )*sqrt(abs(temp1))

!   x, y, and z coordinates (in units of rsm(n) of point where particle
!   enters nucleus, measured with respect to the center of the nucleus
!   with the z-axis being parallel to the beam:

    partin(1) = temp3*cos(temp2)
    partin(2) = temp3*sin(temp2)
    partin(3) = -clientTarg%zoneBRFrac( clientTarg%numZones() ) &
         & * sqrt(abs(one - temp1))
    p0 = sqrt(  abs( clientProj%kinEnergy * &
         & (clientProj%kinEnergy + two*clientProj%restMass) )  )
!  Cartesian components of angular momentum of incoming particle:
    am0(1) =  p0*partin(2)
    am0(2) = -p0*partin(1)
    am0(3) =  zro
!  Incident kinetic energy:
    t3 = clientProj%kinEnergy
!  Total kinetic energy of particle in outside zone of the potential
!   well (tin + pot):
    partin(8) = clientProj%kinEnergy + &
         & poten (clientTarg, clienttarg%numZones(), ipatin)
    return

! ======================================================================
  end subroutine pinpn


  subroutine pinpn1 (sDCM, clientProj, clientTarg, partin, ipatin, &
       & nout, am0, p0, t3)

! ======================================================================
!
!     Calculation of entry point of particle into nucleus, initial
!     angular momentum of nucleus.
!     This version is used with refraction turned on.
!
!     Called by: CASCAD
!
!     Calls: GEOM8 POTEN REFRAC1
!
!  Original PINPN:
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Modified by A. J. Sierk,  LANL  T-2  March, 1996.
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!  New PINPN1:
!    Modified by K. K. Gudima,  Feb., 2004. (refra! ==> refrac1)
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================
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
!                       ipatin(1); charge of particle
!                       ipatin(2); non-zero for photon reactions
!                       ipatin(3); strangeness of particle
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two, twpi
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    type(sDCMProjectile),   intent(inout) :: clientProj
    class(StandardDCMData), intent(inout) :: clientTarg
    real(real64),           intent(  out) :: partin(9)
    integer(int32),         intent(  out) :: ipatin(5)
    integer(int32),         intent(  out) :: nout
    real(real64),           intent(  out) :: am0(3)
    real(real64),           intent(  out) :: p0
    real(real64),           intent(  out) :: t3

    integer(int32) :: irefl, nabs
    real(real64)   :: costin, gkapa, rin, sintin, temp1, temp2, temp3

! ======================================================================

    partin(9) = clientProj%restMass
    ipatin(1) = clientProj%numProtons
!  clientProj%gammaFlag is nonzero for photons:
    ipatin(2) = clientProj%gammaFlag
    ipatin(3) = clientProj%decayNumber
    ipatin(4) = clientProj%numBaryons
10  continue
!   sin(theta) direction of momentum
    partin(4) = zro
!   cos(theta) direction of momentum
    partin(5) = one
!   cos(phi) direction of momentum
    partin(6) = zro
!   sin(phi) direction of momentum
    partin(7) = one
    nout = 0
    ipatin(5) = clientTarg%numZones() + 2
    temp1 = sDCM%rang()
    temp2 = twpi*sDCM%rang()
    temp3 = clientTarg%zoneBRFrac( clientTarg%numZones()+1 )*sqrt(abs(temp1))

!   x, y, and z coordinates (in units of rsm(n) of point where particle
!   enters nucleus, measured with respect to the center of the nucleus
!   with the z-axis being parallel to the beam:

    partin(1) = temp3*cos(temp2)
    partin(2) = temp3*sin(temp2)
    partin(3) = -clientTarg%zoneBRFrac( clientTarg%numZones()+1 ) * &
         & sqrt(abs(one - temp1))
    p0 = sqrt(  abs( clientProj%kinEnergy * &
         & (clientProj%kinEnergy + two*clientProj%restMass) )  )
!  Cartesian components of angular momentum of incoming particle:
    am0(1) =  p0*partin(2)
    am0(2) = -p0*partin(1)
    am0(3) =  zro
!  Incident kinetic energy:
    t3 = clientProj%kinEnergy
!  Total kinetic energy of particle in outside zone of the potential
!   well (tin + pot):
    partin(8) = clientProj%kinEnergy
    irefl = 0
    call sDCM%refrac1 (clientTarg, partin, ipatin, nabs, nout, irefl)
    rin = sqrt(partin(1)**2 + partin(2)**2 + partin(3)**2)
    if (rin < div0Lim .and. rin > -div0Lim) then
       rin = div0Lim
       write(sDCM%io%message,1000) "214"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    costin = (partin(1)*partin(4)*partin(7) + &
         & partin(2)*partin(4)*partin(6) + partin(3)*partin(5))/rin
    sintin = sqrt(abs(one - costin**2))
    if ( (clientTarg%zoneBRFrac( clientTarg%numZones()+1 ) * sintin) &
         &  >=  clientTarg%zoneBRFrac( clientTarg%numZones() ) )  go to 10
    gkapa = sDCM%geom8(clientTarg, partin, ipatin )
    partin(1) = partin(1) + gkapa*partin(4)*partin(7)
    partin(2) = partin(2) + gkapa*partin(4)*partin(6)
    partin(3) = partin(3) + gkapa*partin(5)
    irefl = 0
    call sDCM%refrac1 (clientTarg, partin, ipatin, nabs, nout, irefl)
    if ( ipatin(5) > clientTarg%numZones() ) go to 10

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'pinpn.f90' line(s) ", A)
! ======================================================================
  end subroutine pinpn1
