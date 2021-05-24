
  subroutine paulip (sDCM, clientTarg, partin, ipatin, v, mv, np, irefl, &
       & ipa, pauliData, results)

! ======================================================================
!
!     Pauli principle subroutine. Converts momenta from interaction
!     CM frame to lab frame, and determines if any nucleon final
!     state is Pauli blocked.
!
!    Called by: CASCAD
!
!    Calls: KINEMA
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    More editing, AJS  July, 1997.
!    Modified by AJS, August, 2002.
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================
!
!   The C. M. momentum components are initially stored in pmemo(1:3).
!   They are transformed into the lab frame in KINEMA; then tested to
!   see if the lab kinetic energies are greater than the Fermi kinetic
!   energy.
!
!  Definition of partin:
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
!  Definition of ipatin:
!                       ipatin(1); charge of particle
!                       ipatin(2); non-zero for photon
!                       ipatin(3); strangeness of particle
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: one, massPiPM
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCM),       intent(inout) :: sDCM
    class(StandardDCMData),   intent(inout) :: clientTarg
    real(real64),             intent(inout) :: partin(9)
    integer(int32),           intent(inout) :: ipatin(5)
    real(real64),             intent(in   ) :: v(3)
    integer(int32),           intent(inout) :: mv
    integer(int32),           intent(in   ) :: np
    integer(int32),           intent(inout) :: irefl
    integer(int32),           intent(inout) :: ipa
    type(sDCMPauliInfo),      intent(inout) :: pauliData
    type(StandardDCMResults), intent(inout) :: results

    integer(int32) :: j1, l, mtemp, ntemp
    real(real64)   :: cfi, ct, sfi, st, t, temp, tf
    real(real64)   :: p(3), pstar(3)

! ======================================================================

    j1 = ipatin(5)
    l = 1
10  if (np == 2 .and. l == 2) l = l + 1
    mtemp = mv + l
    pstar(1) = results%pmemo(1,mtemp)
    pstar(2) = results%pmemo(2,mtemp)
    pstar(3) = results%pmemo(3,mtemp)
    call sDCM%kinema (pstar, v, p, ct, st, cfi, sfi, t, results%pmemo(9,mtemp))
    results%pmemo(1,mtemp) = partin(1)
    results%pmemo(2,mtemp) = partin(2)
    results%pmemo(3,mtemp) = partin(3)
    results%pmemo(4,mtemp) = st
    results%pmemo(5,mtemp) = ct
    results%pmemo(6,mtemp) = sfi
    results%pmemo(7,mtemp) = cfi
    results%pmemo(8,mtemp) = t
    results%imemo(5,mtemp) = ipatin(5)
    pauliData%ngen(mtemp) = pauliData%ing + 1
    if (results%imemo(4,mtemp).ne.1 .or. results%imemo(3,mtemp).ne.0) then
!  Not a baryon:
       if (l < np) then
          l = l + 1
          go to 10
       elseif (l >= np .and. ipa == 0) then
          ipa = 2
          return
       else
          go to 20
       endif
    else
!  A baryon:
       temp = results%imemo(1,mtemp)
       tf = clientTarg%protFermiMom(j1) * temp + &
            & clientTarg%neutFermiMom(j1) * (one - temp)
       if (results%pmemo(8,mtemp) > tf) then
          if (l < np) then
!   If all particles haven't been checked, go back.
             l = l + 1
             go to 10
          elseif (l >= np .and. ipa == 0) then
!   Energy was less than cutof3 in CASCAD, the 2 imaginary
!   potentials were more than 30% different, and both nucleons
!   exceeded the Fermi energy:
             ipa = 2
             return
          endif
       else
!   Kinetic energy <= Fermi energy; return to CASCAD with ipa unchanged.
          return
       endif
    endif
20  if (l >= np .and. ipa.ne.0) then
       if ((abs(partin(9)-massPiPM) < 0.01d0 .and. np == 2) .and. &
            & (results%pmemo(9,mv+1) > 0.9d0 .and. &
            & results%pmemo(9,mv+3) > 0.9d0)) pauliData%indi = 3
       partin(4) = results%pmemo(4,mv+3)
       partin(5) = results%pmemo(5,mv+3)
       partin(6) = results%pmemo(6,mv+3)
       partin(7) = results%pmemo(7,mv+3)
       partin(8) = results%pmemo(8,mv+3)
       partin(9) = results%pmemo(9,mv+3)
       ipatin(1) = results%imemo(1,mv+3)
       ipatin(2) = results%imemo(2,mv+3)
       ipatin(3) = results%imemo(3,mv+3)
       ipatin(4) = results%imemo(4,mv+3)
       pauliData%ing = pauliData%ngen(mv+3)
       results%excitons%numExcHoles = results%excitons%numExcHoles + 1
       if (np > 2) then
          ntemp = mv + np
          results%pmemo(4,mv+3) = results%pmemo(4,ntemp)
          results%pmemo(5,mv+3) = results%pmemo(5,ntemp)
          results%pmemo(6,mv+3) = results%pmemo(6,ntemp)
          results%pmemo(7,mv+3) = results%pmemo(7,ntemp)
          results%pmemo(8,mv+3) = results%pmemo(8,ntemp)
          results%pmemo(9,mv+3) = results%pmemo(9,ntemp)
          results%imemo(1,mv+3) = results%imemo(1,ntemp)
          results%imemo(2,mv+3) = results%imemo(2,ntemp)
          results%imemo(3,mv+3) = results%imemo(3,ntemp)
          results%imemo(4,mv+3) = results%imemo(4,ntemp)
          pauliData%meso = 1
       endif
       mv = mv + np - 1
       irefl = 0
    endif
    return

! ======================================================================
  end subroutine paulip
