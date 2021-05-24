
  function poten (clientTarg, i, ipatin)

! ======================================================================
!
!     Calculation of particle potential in  nuclear zone number i.
!
!   Called by: CASCAD PINPN REFRAC
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================
!
!  Definition of ipatin:
!                       ipatin(1); charge of particle
!                       ipatin(2); non-zero for photon interaction
!                       ipatin(3); strangeness of particle
!                       ipatin(4); particle baryon number
!                       ipatin(5); zone number of nucleus where particle
!                                  is located.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one
    use standardDCMDataClass, only: StandardDCMData

    implicit none
    class(StandardDCMData), intent(inout) :: clientTarg
    integer(int32),         intent(in   ) :: i
    integer(int32),         intent(in   ) :: ipatin(5)
    real(real64)                          :: poten

    real(real64) :: temp1, temp2, temp3

! ======================================================================

    if (i-clientTarg%numZones() >= 2) then
       poten = zro
    else
       temp1 = dble(ipatin(1))
       if (i-clientTarg%numZones() == 1) then
!  i = n+1; nuclear part of potential energy vanishes
          poten = temp1*clientTarg%coulombPote( clientTarg%numZones()+1 )
       else
          if (ipatin(3).ne.0) then
!  Coulomb energy only
             poten = temp1*clientTarg%coulombPote(i)
          else
             if (ipatin(4) == 0) then
!  Meson (pion binding + Coulomb energy):
                temp2 = dble(ipatin(2) - 1)
                poten = temp1*clientTarg%coulombPote(i) - &
                     & clientTarg%pionPote() * temp2
             else
!  Baryon (Coulomb + particle sep. eng. + Fermi energy)
                temp3 = dble(ipatin(4))
                poten = temp1*clientTarg%coulombPote(i) &
                     & + clientTarg%protFermiMom(i) * temp1 &
                     & + (one - temp1) * clientTarg%neutFermiMom(i) &
                     & + temp3*clientTarg%getSepEnergy()
             endif
          endif
       endif
    endif
    return
! ======================================================================
  end function poten
