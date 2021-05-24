
  subroutine partn (sDCM, clientTarg, partin, ipatin, partne, ipatne)

! ======================================================================
!
!     Partner selection.
!     partin and ipatin refer to the cascade particle;
!     partne and ipatne refer to the partner nuclear particle which is
!        potentially interacting with the cascade particle.
!
!    Called by: ABSORP POinTE
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Edited by AJS, August, 1997.
!    "Last" change: 13-AUG-2003 by NVMokhov
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
!
!  Definition of ipatin:
!                       ipatin(1); charge of particle
!                       ipatin(2); non-zero for photons
!                       ipatin(3); strangeness of particle; always 0!
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
!  partne, ipatne are corresponding characteristics of the partner.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: one, two, twthrd, emneut, emprot, twpi
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    class(StandardDCMData), intent(inout) :: clientTarg
    real(real64),           intent(in   ) :: partin(9)
    integer(int32),         intent(in   ) :: ipatin(5)
    real(real64),           intent(  out) :: partne(9)
    integer(int32),         intent(  out) :: ipatne(5)

    integer(int32) :: j1
    real(real64)   :: b1, b2, beta1, beta2, bpartn, phin, tnjr

! ======================================================================

    j1 = ipatin(5)
    bpartn = (clientTarg%numBaryons() - clientTarg%numProtons() ) / &
         & clientTarg%numBaryons()
    beta1 = sDCM%rang()
    beta2 = sDCM%rang()**twthrd
!   Choose whether interacting particle is a neutron or a proton:
    if (beta1 >= bpartn) then
!   Proton
       ipatne(1) = 1
       partne(9) = emprot
!   Kinetic energy randomly chosen from Fermi sea.
       tnjr = clientTarg%protFermiMom(j1)*beta2
    else
!  Neutron
       ipatne(1) = 0
       partne(9) = emneut
!   Kinetic energy randomly chosen from Fermi sea.
       tnjr = clientTarg%neutFermiMom(j1)*beta2
    endif
!   The partner has a randomly selected direction of motion
    b1 = sDCM%rang()
    partne(5) = one - two*b1
    b2 = sDCM%rang()
    phin = twpi*b2
    partne(4) = sqrt(abs(one - partne(5)**2))
    partne(7) = cos(phin)
    partne(6) = sin(phin)
    partne(1) = partin(1)
    partne(2) = partin(2)
    partne(3) = partin(3)
    partne(8) = tnjr
    ipatne(5) = ipatin(5)
    ipatne(2) = 0
    ipatne(3) = 0
    ipatne(4) = 1

    return
! ======================================================================
  end subroutine partn
