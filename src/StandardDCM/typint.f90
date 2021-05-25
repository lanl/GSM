
  subroutine typint (sDCM, clientTarg, partin, ipatin, partne, &
       & ipatne, v, u, tin1, sigp, sign, sigabs, mv, np, results)

! ======================================================================
!
!     Determine interaction type and calculate
!     secondary particles' characteristics.
!
!   Called by: CASCAD
!
!   Calls: ABSORP BINEL ELEX SIGMAT8 SLQEK
!
!   CEM95 written by S. G. Mashnik
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Edited by AJS, July-August, 1997.
!   Modified by AJS, December, 1998.
!   Modified by AJS, February, 1999.
!   Modified by AJS, April, 2000.
!   "Last" change: 13-AUG-2003 by NVM
!   Modified by A. J. Sierk, October, 2003.
!   Modified for bremsstrahlung gammas, AJS, February, 2005.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================
!
!  Definition of ipatin (ipatne similar for 2nd particle):
!                      ipatin(1); charge of particle
!                      ipatin(2); non-zero for gammas
!                      ipatin(3); strangeness of particle
!                      ipatin(4); particle baryon number
!                      ipatin(5); nuclear zone number
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams,    only: zro, one
    use standardDCMData,      only: pionProdThresh
    use standardDCMDataClass, only: sDCMData => StandardDCMData

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    class(sDCMData),        intent(inout) :: clientTarg
    real(real64),           intent(in   ) :: partin(9)
    integer(int32),         intent(in   ) :: ipatin(5)
    real(real64),           intent(in   ) :: partne(9)
    integer(int32),         intent(inout) :: ipatne(5)
    real(real64),           intent(  out) :: v(3)
    real(real64),           intent(in   ) :: u
    real(real64),           intent(in   ) :: tin1
    real(real64),           intent(in   ) :: sigp
    real(real64),           intent(in   ) :: sign
    real(real64),           intent(in   ) :: sigabs
    integer(int32),         intent(in   ) :: mv
    integer(int32),         intent(  out) :: np
    class(StandardDCMResults), intent(inout) :: results

    integer(int32) :: id, ksi, l, mtemp, mb, me, ms, nin
    real(real64)   :: b1, b2, betabs, betael, r1, sig1pi, sigdif, sigel, &
         & sigelex, sigex, sigtot, temp

! ======================================================================

    ! Create photon cross section data for the event, if needed
    !    (replace the /isecgpn/ block)
    type(sDCMPhotonCrossSections) :: photoData

!  Thresholds for pi0 induced charge exchange reactions (1st 2);
!  "Effective" thresholds for pion- and nucleon-induced 2 pion
!  production (3 & 4). AJS, 4/13/00, 12/22/98.
! NOTE: 'thresh' moved to sDCMParams module
    real(real64), parameter, dimension(2) :: &
         & cxthr =  [ 0.003783_real64, 0.006755_real64 ]


!  "Effective" threshold for 1 or 2 pion production by gammas:
    real(real64), parameter, dimension(2) :: &
         & gamthr = [ 0.152_real64, 0.321_real64 ]

! ======================================================================

    mtemp = ipatin(5)
!   Betabs will be zero unless cascade particle is a pion or gamma.
    temp = clientTarg%protonDensity(mtemp) * sigp + &
         & clientTarg%neutronDensity(mtemp) * sign + &
         & clientTarg%protonDensity(mtemp) * sigabs
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "72"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    betabs = clientTarg%protonDensity(mtemp)*sigabs/(temp)
10  b1 = sDCM%rang()
    if (b1 <= betabs) then
       call sDCM%absorp (clientTarg, partin, ipatin, partne, ipatne, &
            & mv, np, v, results)
    else
       call slqek (l, ms, mb, ksi, me, ipatin(2), ipatin(3), ipatin(4), &
            & ipatin(1), ipatne(2), ipatne(3), ipatne(4), ipatne(1))
       sigtot = sDCM%sigmat8 (l, ms, mb, ksi, 0, tin1)
       sigex  = sDCM%sigmat8 (l, ms, mb, ksi, 2, tin1)
       sigel  = sDCM%sigmat8 (l, ms, mb, ksi, 1, tin1)
!  For photons, sigex and sigel are single pion production with and
!  without charge change of the nucleon.
       if (ipatin(1) == 0 .and. ipatin(4) == 0 .and. ipatin(2) == 0) then
!  pi-0 induced charge exchange thresholds:
          if (tin1 <= cxthr(ipatne(1)+1)) then
             sigex = zro
             sigel = sigtot
          endif
       endif
       sigel = max(sigel, zro)
       sigelex = sigex + sigel
       sig1pi = max (zro, sDCM%sigmat8 (l, ms, mb, ksi, 7, tin1))


       !  Apply photon thresholds:
       if (ipatin(2).ne.0) then
          if (tin1 < gamthr(1)) then
             ! Below threshold, no reaction (no cross section)
             sigtot  = zro
             sigex   = zro
             sigel   = zro
             sigelex = zro
             sig1pi  = zro
             np = 0
             return
          endif
          if (tin1 < gamthr(2)) then
             sig1pi  = zro
             sigelex = sigtot
          endif
       endif

       ! Obtain cross section difference
       sigdif = sigtot - sigelex


!  Betael = ratio of sum of (elastic cross section + SCX) to total:
!  For photons, ratio of single pion production to total:
       id = ipatin(4) + 1
       if (sigtot < div0Lim .and. sigtot > -div0Lim) then
          sigtot = div0Lim
          write(sDCM%io%message,1000) "119, 121"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       if (tin1 <= pionProdThresh(id) .and. ipatin(2) == 0) then
          betael = one - sig1pi/sigtot
       else
          betael = sigelex/sigtot
          betael = min (one, betael)
       endif
       b2 = sDCM%rang()
       if (b2 <= betael) then
          if (ipatin(2) > 0) then
             call sDCM%setPhotoChannelSigma (tin1, photoData)
          endif
!   Elastic and SCX (1 meson produced in gamma reaction):
!   Criterion for ELEX is 2-body final state.
          call sDCM%elex (v, u, tin1, partin, ipatin, partne, ipatne, mv, np, &
               & l, mb, ksi, me, sigex, sigelex, photoData, results)
       else
!   Inelastic reaction of some type (pion production, 2 pi for gamma-N):
!   n-body final state with n > 2.
          r1 = sDCM%rang()
          call sDCM%binel (partin, ipatin, ipatne, l, ms, mb, ksi, me, v, &
               & u, tin1, mv, np, nin, &
               & clientTarg%numProtons(), clientTarg%numBaryons(), &
               & sigdif, sig1pi, r1, results)
          if (nin.ne.0) go to 10
       endif
    endif


    return

! ======================================================================
1000 format("Divide by zero error prevented in 'typint.f90', line(s) ", A)
! ======================================================================
  end subroutine typint
