
  subroutine isValidParticle(a, z, e_k, cm0, valid)

! ====================================================================
!
! Created to ensure that only valid residual particles are created
! during the IntraNuclear Cascade of LAQGSM (time dependent DCM)
!
! Created by CMJ, XCP-3, Oct. 2017 for LAQGSM Expansion
! Modifed by CMJ, XCP-3, OCt. 2017 for bad residual allowance
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMParams, only: zro, one

    implicit none
    integer(int32), intent(in   ) :: a       ! atomic mass of fragment
    integer(int32), intent(in   ) :: z       ! mass number of fragment (charge)
    real(real64),   intent(in   ) :: e_k     ! kinetic energy of fragment (gev)
    real(real64),   intent(in   ) :: cm0     ! rest mass of fragment (gev)
    logical,        intent(  out) :: valid   ! =true if fragment is physical (valid)

! ====================================================================

    valid = .false.

    if ( a >= zro .and. z >= -one ) then
       if ( a > z .or. (z <= one .and. a <= one) ) then
          if ( e_k > zro .and. cm0 >= zro ) then
             ! Physical residual produced, tally it
             valid = .true.
          else
             ! Negative kinetic energy or rest mass, unphysical
             write(*,1100) e_k, cm0
          endif
       else
          ! Z > A, or Z too large
          write(*,1200) a, z
       endif
    else
       ! Negative values (electron or unphysical residual)
       write(*,1200) a, z
    endif

    return
! ====================================================================
1100 format (3x, 'warning: a fragment was produced in the dcm', &
          & ' with a negative kinetic energy (', es11.9, ') or rest mass(', &
          & es11.9, ').  discarding fragment.')
1200 format (3x, 'warning: an unphysical fragment (a = ', i3, &
          & ', z = ', i3, ') was produced in the dcm.  discarding fragment.')
1500 format (5x, 'comment: allowing residual per user options.')
! ====================================================================
  end subroutine isValidParticle


  subroutine tallyMDCMProgeny(ip, p)

! ====================================================================
!
! Created to simlify how the parz, spt, and sptr arrays are assigned
! in the LAQGSM cascade.  Desire to do so in one command ( call getLA...)
!
! Created by CMJ, XCP-3, July 2017 for LAQGSM Expansion
! CMJ, XCP-3, 11/2017 (Removed passing of parz/spt/sptr)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMParams, only: zro, hlf, one, two, twpi
    use modifiedDCMClass, only: results

    implicit none
    integer(int32), intent(in   ) :: ip(5)
    real(real64),   intent(inout) :: p(9)

    integer(int32) :: indx, numProg
    real(real64)   :: parCharge, fi1, stheta, ctheta, sphi, cphi, totLinearMom, &
         & xyLinearMom, theta

! ====================================================================

    ! uses 'rm2' ONLY!
    real(real64) ::  anucl1, anucl2, znucl1, znucl2, t0, eps1, eps2, vpi, &
         & a1, a2, c1, c2, d1, d2, r0n1, r0n2, tf01, tf02, rm1, rm2
    common /hcasc/ anucl1,anucl2,znucl1,znucl2,t0,eps1,eps2,vpi,a1,a2, &
         & c1,c2,d1,d2,r0n1,r0n2,tf01,tf02,rm1,rm2

      ! For anti-lab correction
    integer(int32) :: kobr
    real(real64)   :: blab, glab
    common/bglab/ blab, glab, kobr

      ! For X/Y/Z normalization
    real(real64)   :: tint
    common/actim/ tint
    real(real64)   :: tprod
    common/tprod/ tprod(5999)

! ====================================================================

    ! Tallying residual
    ! Ensure k in limits (checked from previous routines)
    if ( results%numProgeny > results%maxProgenyM1 ) return
    results%numProgeny = results%numProgeny + 1

    ! Normalize X/Y/Z
    do indx = 1, 3
       !     x = x - px*dT / E_tot, same for y/z (change x's to y, z)
       p(indx) = p(indx) - p(indx+3) * (tint - tprod(results%numProgeny)) / &
            & ( p(8) + p(9) )
    end do

    ! Correct for anti-lab system of INC simulation
    if ( kobr == 1 ) then
       p(6) = - glab * ( p(6) - blab*( p(8) + p(9) )  )
    end if

    ! Obtain angles from momenta in P() array, and remove momenta data w/ angular information
    totLinearMom = sqrt( p(4)**2 + p(5)**2 + p(6)**2 ) ! total momentum
    ! Obtain theta information
    ctheta = p(6) / totLinearMom ! cos(theta) = pz/ptot
    stheta = sqrt( 1 - ctheta**2 ) ! sin(theta) = sqrt[ 1 - cos(theta)^2 ]
    ! Obtain phi information
    xyLinearMom = totLinearMom * stheta ! ptot * sin(theta)
    cphi = p(4) / ( xyLinearMom ) ! cos(phi) = px / [ptot * sin(tehta)]
    sphi = p(5) / ( xyLinearMom ) ! sin(phi) = py / [ptot * sin(tehta)]

    theta = atan2( stheta, ctheta )
    fi1   = atan2( sphi,   cphi   )
    if (fi1 < zro) fi1 = twpi + fi1


    ! Add progeny to the results object:
    numProg = results%numProgeny
    ! Particle ID:
    results%progenyBnk(numProg)%numBaryons    = ip(4)
    results%progenyBnk(numProg)%numProtons    = ip(1)
    results%progenyBnk(numProg)%strangeness   = ip(3)
    results%progenyBnk(numProg)%photonFlag    = ip(2)
    ! Energy:
    results%progenyBnk(numProg)%kinEnergy     = p(8)
    results%progenyBnk(numProg)%restMass      = p(9)
    ! Position:
    results%progenyBnk(numProg)%position(1)   = p(1)
    results%progenyBnk(numProg)%position(2)   = p(2)
    results%progenyBnk(numProg)%position(3)   = p(3)
    ! Linear Momentum:
    results%progenyBnk(numProg)%linearMom(1)  = p(4)
    results%progenyBnk(numProg)%linearMom(2)  = p(5)
    results%progenyBnk(numProg)%linearMom(3)  = p(6)
    ! Emission Angles:
    results%progenyBnk(numProg)%theta         = theta
    results%progenyBnk(numProg)%phi           = fi1
    ! Time values:
    results%progenyBnk(numProg)%formationTime = p(7)
    results%progenyBnk(numProg)%timeParameter = ip(5)

    return
! ====================================================================
1100 format (/3x,'warning: number of particles (', i3, &
          & ') almost exceeded in cascaw, trial number ',i7,'.  ', i4, &
          & ' created.')
1200 format (/5x,'error: number of particles (', i3, ') exceeded ', &
          & 'in cascaw, trial number ',i7,'.  ', i4, ' created.')
1400 format (5x,'fragment ', i4, ' of ', i6, ': a = ', i3, ', z = ', &
          & i3, ', e_k = ', es9.3, ' mev, and e_0 (rest mass) = ', es9.3, &
          & ' gev.')
1500 format (3x, 'comment: ', a)
! ====================================================================
  end subroutine tallyMDCMProgeny


  subroutine resetmainvars()

! ====================================================================
!
! Resets parz and spt arrays that are passed in to routine to 0.
!
!
! Made by CMJ, XCP-3, 2016-2017, LAQGSM Expansion
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: real64
    implicit none

! ====================================================================

    real(real64) ::  parz, spt
    common /blok77/ spt(5,150)
    common /zapp/   parz(6,150)

! ====================================================================

    parz(:,:) = 0.0_real64
    spt(:,:)  = 0.0_real64

    return
! ====================================================================
  end subroutine resetmainvars
