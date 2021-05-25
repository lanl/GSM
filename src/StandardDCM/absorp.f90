
  subroutine absorp (sDCM, clientTarg, partin, ipatin, partne, &
       & ipatne, mv, np, v, results)

! ======================================================================
!
!   Calculates outgoing nucleon characteristics in pion or gamma
!   absorption.
!
!   Called by: TYPINT
!
!   Calls: ABEL CHABS PARTN TINVU
!
!   CEM95 written by S. G. Mashnik
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Modified by AJS, July-August, 1997.
!   Modified by AJS, March, 1999.
!   Edited by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by AJS, LANL T-2, February, 2009.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, LANL XCP-3, July 2013 (included error protection).
!
! ======================================================================
!
!  np is the number of particles in the final state (= 2, unless
!  dimensions of storage array are exceeded).
!
!  Definition of partin:
!                       partin(1); Normalized x coordinate of pion
!                       partin(2); Normalized y coordinate of pion
!                       partin(3); Normalized z coordinate of pion
!                       partin(4); sin(theta), direction of momentum
!                       partin(5); cos(theta), direction of momentum
!                       partin(6); sin(phi), direction of momentum
!                       partin(7); cos(phi), direction of momentum
!                       partin(8); kineti! energy of pion
!                       partin(9); rest mass of pion
!
!  partne is similar to partin for first nucleon of the pair which
!  absorbs the pion; ipatne is similar to ipatin.
!
!  Definition of ipatin:
!                       ipatin(1); charge of pion
!                       ipatin(2); non-zero for photon interactions
!                       ipatin(3); strangeness of pion (= 0!)
!                       ipatin(4); particle baryon number (= 0!)
!                       ipatin(5); zone number of nucleus where
!                                  interaction occurs.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two, twpi, emneut, emprot
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    class(StandardDCMData), intent(inout) :: clientTarg
    real(real64),           intent(in   ) :: partin(9)
    integer(int32),         intent(in   ) :: ipatin(5)
    real(real64),           intent(in   ) :: partne(9)
    integer(int32),         intent(inout) :: ipatne(5)
    integer(int32),         intent(in   ) :: mv
    integer(int32),         intent(  out) :: np
    real(real64),           intent(  out) :: v(3)
    class(StandardDCMResults), intent(inout) :: results

    integer(int32) :: ie
    real(real64)   :: b1, b2, cfn1, ctn1, ctst, dot, e1, e2, em1, em2, em2n, &
         & fist, pafm, pn, pn1, r1, sfn1, stn1, taf, temp, temp1, temp2, tin1, &
         & twom1, twom2, u
    real(real64),   dimension(9) :: par1 = zro
    integer(int32), dimension(5) :: ipa1 = 0_int32
    real(real64),   dimension(3) :: paf = zro, pist = zro, pnst = zro

! ======================================================================

    em1 = emneut*(one - ipatne(1)) + emprot*ipatne(1)
    twom1 = two*em1
    temp = partne(8)*(partne(8) + twom1)
    if (temp < 0.0d0 ) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "77"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    pn = sqrt(temp)
!   Pick 2nd nucleon partner for absorption.
    call sDCM%partn (clientTarg, partin, ipatin, par1, ipa1)
    em2 = emneut*(one - ipa1(1)) + emprot*ipa1(1)
    twom2 = two*em2
    temp = par1(8)*(par1(8) + twom2)
    if (temp < 0.0d0 ) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "87"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    pn1 = sqrt(temp)
!   Momentum components of 2-nucleon system:
    paf(1) = pn*partne(4)*partne(7) + pn1*par1(4)*par1(7)
    paf(2) = pn*partne(4)*partne(6) + pn1*par1(4)*par1(6)
    paf(3) = pn*partne(5) + pn1*par1(5)
    temp = paf(1)**2 + paf(2)**2 + paf(3)**2
    if (temp < 0.0d0 ) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "97"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    pafm = sqrt(temp)
!   Cosine of angle between velocity of 2-nucleon CM system and beam
!   direction:
    temp = pafm
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "105"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    ctn1 = paf(3)/temp
    temp1 = one - ctn1**2
    if (temp1 <= zro) then
       ctn1 = one
       stn1 = zro
       sfn1 = zro
       cfn1 = one
    else
       stn1 = sqrt(temp1)
       temp2 = pafm*stn1
       if (temp2 < div0Lim .and. temp2 > -div0Lim) then
          temp2 = div0Lim
          write(sDCM%io%message,1000) "119, 120"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       sfn1 = paf(2)/temp2
       cfn1 = paf(1)/temp2
    endif
!   Find invariant mass of 2-nucleon system:
    temp = pn1**2 + em1**2
    if (temp < 0.0d0 ) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "128"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    e1 = sqrt(temp)
    temp = pn**2 + em2**2
    if (temp < 0.0d0 ) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "134"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    e2 = sqrt(temp)
    dot = pn*pn1*(partne(4)*par1(4)*(partne(7)*par1(7) + &
         & partne(6)*par1(6)) + partne(5)*par1(5))
    temp = em1**2 + em2**2 + two*(e1*e2 - dot)
    if (temp < 0.0d0 ) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "142"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    em2n = sqrt(temp)
!   Kinetic energy of 2N system w. r. t. lab frame.
    temp = pafm**2 + em2n**2
    if (temp < 0.0d0 ) then
       temp = 0.01d0
       write(sDCM%io%message,1100) "149"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    taf = sqrt(temp) - em2n
!   Find kinetic energy tin1 of pion in 2N CM frame, and velocity (v)
!   of 3-particle CM frame w. r. t. lab.
    call sDCM%tinvu (partin(8), partin(9), taf, em2n, partin(4), partin(5), &
         & partin(6), partin(7), stn1, ctn1, sfn1, cfn1, v, u, tin1)
!   Find isospin of outgoing channel (p-p, n-p, or n-n)
!   New method; keep track of original partner's isospin:
    r1 = sDCM%rang()
    call sDCM%chabs (ipatin(2), ipatin(1), ipatne(1), ie, &
         & clientTarg%numBaryons(), clientTarg%numProtons(), r1)
    b1 = sDCM%rang()
    ctst = one - two*b1
!   Low-energy gamma case:
    if (ipatin(2).ne.0 .and. tin1 <= 0.455d0) ctst = sDCM%costa (19, tin1, b1)
    em1 = emneut*(one - ipatne(1)) + ipatne(1)*emprot
    em2 = emneut*(one - ie) + ie*emprot
!   Define a random direction in space. (Direction of 2-particle
!   final state momenta in their rest frame, CM')
    b2 = sDCM%rang()
    fist = twpi*b2
    call sDCM%abel (partin, v, u, pist, pnst, ctst, fist, em1, em2)
!   In ABEL, we find pist and pnst; the oppositely directed momenta
!   of the 2 nucleons in the final state, in their center-of-mass frame.
    if (mv <=  results%maxProgenyM3 ) then
       results%pmemo(1, mv+3) = pist(1)
       results%pmemo(2, mv+3) = pist(2)
       results%pmemo(3, mv+3) = pist(3)
       results%pmemo(9, mv+3) = em2
       results%imemo(1, mv+3) = ie
       results%imemo(2, mv+3) = 0
       results%imemo(3, mv+3) = 0
       results%imemo(4, mv+3) = 1
       results%pmemo(1, mv+1) = pnst(1)
       results%pmemo(2, mv+1) = pnst(2)
       results%pmemo(3, mv+1) = pnst(3)
       results%pmemo(9, mv+1) = em1
       results%imemo(1, mv+1) = ipatne(1)
       results%imemo(2, mv+1) = 0
       results%imemo(3, mv+1) = 0
       results%imemo(4, mv+1) = 1
       np = 2
    else
       np = 0
    endif
    return

! ======================================================================
1000 format("Divide by zero error preveted in 'absrob.f90' line(s) ", A)
1100 format("Square root error preveted in 'absrob.f90' line(s) ", A)
! ======================================================================
  end subroutine absorp
