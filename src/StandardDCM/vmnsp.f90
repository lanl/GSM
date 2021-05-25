
  subroutine vmnsp (sDCM, partin, ipatin, u, mv, np, ith, mb, tin1, &
       & lp, results)

! ======================================================================
!
!     Calculate secondary particle number and determine absolute
!     values of momenta in inelastic (pion production) interaction.
!     For gamma reactions, involves production of two pions.
!
!   Called by: BINEL
!
!   Calls: JTYPB PMOM
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Edited by AJS, July-August, 1997.
!   Modified by AJS, September, 1998.
!   Edited by AJS, December, 1998.
!   Modified by AJS, March, 1999.
!   Edited by AJS, LANL T-2, February, 2009.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================
!
!  np is the number of particles in the final state; set to 0 if the
!     block "memory" is overfilled.
!  lp is set to 2 if the subroutine iterates 100 times without
!     finding a satisfactory solution.
!
!  Definition of partin (pmemo similar):
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
!   pmemo(4,ltemp) is temporarily the total momentum of the particle.
!
!  Definition of ipatin (imemo similar):
!                       ipatin(1); charge of particle
!                       ipatin(2); non-zero for photon interactions
!                       ipatin(3); strangeness of particle
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, two, massPiPM, emnucg

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    real(real64),   intent(in   ) :: partin(9)
    integer(int32), intent(in   ) :: ipatin(5)
    real(real64),   intent(in   ) :: u
    integer(int32), intent(in   ) :: mv
    integer(int32), intent(  out) :: np
    integer(int32), intent(in   ) :: ith
    integer(int32), intent(in   ) :: mb
    real(real64),   intent(in   ) :: tin1
    integer(int32), intent(  out) :: lp
    class(StandardDCMResults), intent(inout) :: results

    integer(int32) :: i, jb, kh, lambda, ltemp
    real(real64)   :: deltu, el, pmax, r1, sigma, temp, u1

! ======================================================================

    kh = 0
    lp = 0
10  continue
    u1 = u
    lambda = 1
20  continue
    ltemp = mv + lambda
    if (ltemp > results%maxProgeny) then
       np = 0
       return
    endif
    if (lambda == 1) then
!  The first particle is ALWAYS a nucleon (from the Fermi sea)
!  (Possibly charge exchanged?) Charge and exact mass aren't determined
!  until AFTER calling CHINEL!
       results%pmemo(9,ltemp) = emnucg
       results%imemo(2,ltemp) = 0
       results%imemo(3,ltemp) = 0
       results%imemo(4,ltemp) = 1
    elseif (lambda == 3) then
!  "Projectile" (possibly charge exchanged(?); either pion or nucleon)
       results%pmemo(9,ltemp) = partin(9)
       results%imemo(2,ltemp) = 0
       results%imemo(3,ltemp) = ipatin(3)
       results%imemo(4,ltemp) = ipatin(4)
    else
!  Produced pion:
       results%pmemo(9,ltemp) = massPiPM
       results%imemo(2,ltemp) = 0
       results%imemo(3,ltemp) = 0
       results%imemo(4,ltemp) = 0
    endif
    jb = jtypb (ith, mb, lambda)
!  results%pmemo(4,...) is temporarily filled with the 3-momentum
!  of the particle.
    r1 = sDCM%rang()
    results%pmemo(4,ltemp) = pmom (jb, tin1, r1)
    el = sqrt(results%pmemo(4,ltemp)**2 + results%pmemo(9,ltemp)**2)
!  el is the total energy of the particle; deltu is remaining
!  total energy, after the particle is removed.
    deltu = u1 - el
    if (lambda == 2) then
       if (deltu <= partin(9)) then
!  Not enough energy left for "projectile"; rerun the subroutine.
          kh = kh + 1
          if (kh < 100) go to 10
!  After 100 iterations, rerun TYPINT
          lp = 2
          np = 1
          return
       endif
       if (ith.ne.0) then
!   Single pion produced:
          results%pmemo(4,mv+3) = sqrt(abs(deltu**2 - partin(9)**2))
          results%pmemo(9,mv+3) = partin(9)
          results%imemo(2,mv+3) = 0
          results%imemo(3,mv+3) = ipatin(3)
          results%imemo(4,mv+3) = ipatin(4)
          temp = results%pmemo(4,mv+1)-results%pmemo(4,mv+2)-results%pmemo(4,mv+3)
          if ( temp  <=  zro) then
             temp = abs( results%pmemo(4,mv+2)-results%pmemo(4,mv+3) )
             if (results%pmemo(4,mv+1)  >  temp ) then
                np = 3
                return
             endif
          endif
          kh = kh + 1
          if (kh < 100) go to 10
!  After 100 iterations, rerun TYPINT
          lp = 2
          np = 1
          return
       endif
       u1 = deltu
       lambda = lambda + 1
       go to 20
    endif
    if (deltu > massPiPM) then
       u1 = deltu
       lambda = lambda + 1
       go to 20
    endif
    if (lambda <= 1 .or. lambda == 3) then
       kh = kh + 1
       if (kh < 100) go to 10
!  After 100 iterations, rerun TYPINT
       lp = 2
       return
    endif
    el = deltu + el
!  Reset momentum of last pion emitted to use up the residual
!  of total energy.
    results%pmemo(4,ltemp) = sqrt(abs(el**2 - results%pmemo(9,ltemp)**2))
    np = lambda
!     i = 1
!  Find maximum individual particle 3-momentum (magnitude).

    pmax = results%pmemo(4,mv+1)
    sigma = results%pmemo(4,mv+1)
    do i = 2,np
       pmax = max(pmax, results%pmemo(4,mv+i))
       sigma = sigma + results%pmemo(4,mv+i)
    end do

    if (two*pmax >= sigma) then
       kh = kh + 1
       if (kh < 100) go to 10
!  After 100 iterations, rerun TYPINT
       lp = 2
       np = 1
    endif
    return

! ======================================================================
  end subroutine vmnsp
