
  function geom8 (sDCM, clientTarg, partin, ipatin)

! ======================================================================
!
!     Calculation of geometrical path of particle.
!     Geom is distance from current point along trajectory to the first
!     intersection of the trajectory with the edge of the current
!     nuclear zone.
!
!    Called by: POINTE REFRAC1
!
!    CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk   LANL  T-2,  February, 1996.
!    Edited by AJS, December, 1997.
!    Modified by AJS to remove a very rare problem for very small
!    ctinjb, March, 1999.
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
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
!
!  Definition of ipatin:
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    class(StandardDCMData), intent(inout) :: clientTarg
    real(real64),           intent(in   ) :: partin(9)
    integer(int32),         intent(inout) :: ipatin(5)
    real(real64)                          :: geom8

    integer(int32) :: iq, kg
    real(real64)   :: ctinjb, deltag, g, rin, rl, stinjb

! ======================================================================

    rin = sqrt(partin(1)**2 + partin(2)**2 + partin(3)**2)
    if (abs(rin) < div0Lim) then
       write (sDCM%io%message, 2000) rin
       call sDCM%io%print(3, 3, sDCM%io%message)
       rin = div0Lim
    end if
    if (rin.ne.zro) then
       ctinjb = (partin(1)*partin(4)*partin(7) + partin(2)*partin(4)* &
            & partin(6) + partin(3)*partin(5)) / rin
       if (abs(ctinjb) > one) ctinjb = sign (one, ctinjb)
    else
!  Particle at origin, necessarily moving outward.
       ctinjb = one
    endif
    stinjb = sqrt(abs(one - ctinjb**2))

    if (ctinjb < zro) then
!  Particle moving toward center of nucleus (r decreasing).
       g = -one
       kg = ipatin(5) - 1
       if (kg <= 0) then
          rl = zro
       else
          rl = clientTarg%zoneBRFrac(kg)
       endif
    else
!  Particle moving away from center of nucleus (r increasing).
       g = one
       kg = ipatin(5)
       rl = clientTarg%zoneBRFrac(kg)
    endif
    iq = 0
10  deltag = rl**2 - (rin*stinjb)**2
    iq = iq + 1
    if (iq > 2 .and. deltag < zro) then
       if (abs(ctinjb) <= div0Lim) then
          if (abs(rl - rin) <= div0Lim) then
!   Particle is at outside boundary of zone, moving nearly tangentially.
             ipatin(5) = ipatin(5) + 1
             kg = ipatin(5)
             rl = clientTarg%zoneBRFrac(kg)
             go to 10
          elseif (abs(clientTarg%zoneBRFrac(kg+1) - rin) <= div0Lim) then
!   Particle is at outside boundary of next outer zone, moving nearly
!   tangentially.
             ipatin(5) = ipatin(5) + 2
             kg = ipatin(5)
             rl = clientTarg%zoneBRFrac(kg)
             go to 10
          endif
       else
          ipatin(1) = 2
          geom8 = zro
          return
       endif
    endif
    if (deltag < zro) then
!   Does not intersect inner edge of zone, look for intersection with
!   outer edge.
       g = one
       kg = ipatin(5)
       rl = clientTarg%zoneBRFrac(kg)
       go to 10
    else
!  For g = -1, distance along trajectory to intersection with inner
!  edge of zone; for g = +1, distance to intersection with outer edge.
       geom8 = g*sqrt(abs(deltag)) - rin*ctinjb
    endif
    return

! ======================================================================
2000 format ("Particle origin is r=", 1pe12.4, " [fm]. Preventing ", &
          & "divide by zero error in 'geom8.f90'.")
! ======================================================================
  end function geom8
