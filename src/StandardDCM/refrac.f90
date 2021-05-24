
  subroutine refrac (sDCM, clientTarg, partin, ipatin, nabs, nout)

! ======================================================================
!
!     Calculation of energy and direction change of a particle
!     when entering another zone of the nucleus.
!
!    NOTE: NO DIRECTION CHANGE CALCULATED!!!
!
!    Called by: POINTE
!
!    Calls: POTEN
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February-March, 1996.
!   "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================
!
!    nout = 1 if particle is in a zone outside the nucleus.
!           [ipatin(5) > n+1]
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
    use standardDCMParams, only: zro
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    class(StandardDCMData), intent(inout) :: clientTarg
    real(real64),           intent(inout) :: partin(9)
    integer(int32),         intent(inout) :: ipatin(5)
    integer(int32),         intent(  out) :: nabs
    integer(int32),         intent(  out) :: nout

    integer(int32) ::  ip5, jr1
    real(real64)   :: ctinjb, rin, tinj1

! ======================================================================

    nabs = 0
    nout = 0
    ip5 = ipatin(5)
    rin = sqrt(partin(1)**2 + partin(2)**2 + partin(3)**2)
    if (rin < div0Lim .and. rin > -div0Lim) then
       rin = div0Lim
       write(sDCM%io%message,1000) "68"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    ctinjb = (partin(4)*(partin(1)*partin(7) + &
         & partin(2)*partin(6)) + &
         & partin(3)*partin(5))/rin
!   ctinjb < 0 means inward-going particle
    if (ctinjb < zro) then
       jr1 = ip5 - 1
    else
       jr1 = ip5 + 1
    endif
    tinj1 = partin(8) + poten (clientTarg, jr1, ipatin) &
         & - poten (clientTarg, ip5, ipatin)
    if (tinj1 <= zro) nabs = 1
    partin(8) = tinj1
    ipatin(5) = jr1
    if (jr1 - clientTarg%numZones()  >=  2) nout = 1

    return
! ======================================================================
1000 format("Divide by zero error prevented in 'refrac.f90', line(s) ", A)
! ======================================================================
  end subroutine refrac


  subroutine refrac1 ( sDCM, clientTarg, partin, ipatin, nabs, nout, irefl)

! ======================================================================
!
!     Calculation of energy and direction change of a particle
!     when entering another zone of the nucleus.
!
!    NOTE: DIRECTION IS CHANGED !!!
!
!    Called by: PINPN1 POINTE1
!
!    Calls: POTEN ROTOR
!
!   Original no refraction version, REFRAC:
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February-March, 1996.
!   "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!   Modified version with refraction, REFRAC1:
!    Modified by K. K. Gudima,  Feb., 2004.
!    Edited by A. J. Sierk, LANL T-16, April, 2004.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================
!
!    nout = 1 if particle is in a zone outside the nucleus.
!              [ipatin(5) > n+1]
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
    use standardDCMParams, only: zro, one, two
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    class(StandardDCMData), intent(inout) :: clientTarg
    real(real64),           intent(inout) :: partin(9)
    integer(int32),         intent(inout) :: ipatin(5)
    integer(int32),         intent(  out) :: nabs
    integer(int32),         intent(  out) :: nout
    integer(int32),         intent(inout) :: irefl

    integer(int32) :: ip5, jr1
    real(real64)   :: ctibj1, ctinjb, pinj, pinj1, rin, stibj1, &
         & stinjb, temp, tinj1
    real(real64), dimension(3) :: a = zro, b = zro, cstar = zro, c = 0

! ======================================================================

    nabs = 0
    nout = 0
    ip5 = ipatin(5)
    rin = sqrt(partin(1)**2 + partin(2)**2 + partin(3)**2)
    if (rin < div0Lim .and. rin > -div0Lim) then
       rin = div0Lim
       write(sDCM%io%message,1000) "157"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    ctinjb = (partin(4)*(partin(1)*partin(7) + &
         & partin(2)*partin(6)) + &
         & partin(3)*partin(5))/rin
!   ctinjb < 0 means inward-going particle
    if (ctinjb < zro) then
       jr1 = ip5 - 1
    else
       jr1 = ip5 + 1
    endif
    tinj1 = partin(8) + &
         & poten (clientTarg, jr1, ipatin) - &
         & poten (clientTarg, ip5, ipatin)
    if (partin(8).ne.tinj1) then
       if (tinj1 <= zro) then
          nabs = 1
          go to 10
       endif
       stinjb = sqrt(one - ctinjb**2)
       pinj  = sqrt(abs(partin(8)*(partin(8) + two*partin(9))))
       pinj1 = sqrt(abs(tinj1*(tinj1 + two*partin(9))))
       if (pinj1 < div0Lim .and. pinj1 > -div0Lim) then
          pinj1 = div0Lim
          write(sDCM%io%message,1000) "178" 
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       stibj1 = (pinj/pinj1)*stinjb
       a(1) = partin(1)
       a(2) = partin(2)
       a(3) = partin(3)
       b(1) = pinj*partin(4)*partin(7)
       b(2) = pinj*partin(4)*partin(6)
       b(3) = pinj*partin(5)
       if (stibj1 <= one) then
          ctibj1 = sign(one,ctinjb)*sqrt(one - stibj1**2)
          cstar(1) = pinj1*stibj1
          cstar(2) = zro
          cstar(3) = pinj1*ctibj1
       else
          tinj1 = partin(8)
          jr1 = ip5
          pinj1 = pinj
          cstar(1) = pinj1*stinjb
          cstar(2) = zro
          cstar(3) =-pinj1*ctinjb
          irefl = irefl + 1
          if (irefl > 3)  then
             irefl = 0
             nabs = 1
             return
          endif
       endif
       call sDCM%rotor (a, b, cstar, c)
       if (pinj1 < div0Lim .and. pinj1 > -div0Lim) then
          pinj1 = div0Lim
          write(sDCM%io%message,1000) "209"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       partin(5) = c(3)/pinj1
       if (abs(partin(5)) >= one) then
          partin(4) = zro
       else
          partin(4) = sqrt(abs(one - partin(5)**2))
       endif
       if (partin(4).ne.zro) then
          temp = pinj1*partin(4)
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "223, 224"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          partin(7) = c(1)/(temp)
          partin(6) = c(2)/(temp)
       else
          partin(7) = one
          partin(6) = zro
       endif
10     partin(8) = max(zro, tinj1)
    endif
    ipatin(5) = jr1
    if (jr1 - clientTarg%numZones()  >=  2) nout = 1
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'refrac.f90', line(s) ", A)
! ======================================================================
  end subroutine refrac1
