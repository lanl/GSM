
  subroutine chabs (sDCM, l, ine, ne1, ne2, a, z, r1)

! ======================================================================
!
!     Determines isospin (p or n) of nucleons after pion absorption.
!     This modified version keeps track of the isospin of the original
!     first nucleon partner.
!
!   Called by: ABSORP
!
!       ine; charge of incident pion or photon
!       ne1; upon entry: charge of first nucleon partner;
!       ne1, ne2; upon return: the charges of the two outgoing nucleons.
!       l = 1 for gammas, otherwise = 0
!
!   CEM95 written by S. G. Mashnik
!   Edited by A. J. Sierk, LANL T-2, February, 1996.
!   Modified  from old CHABS by AJS, July, 1997
!   Added gamma absorption by FCG August, 2000.
!   "Last" change: 12-AUG-2003 by NVMokhov
!   Modified by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: one

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    integer(int32), intent(in   ) :: l
    integer(int32), intent(in   ) :: ine
    integer(int32), intent(inout) :: ne1
    integer(int32), intent(inout) :: ne2
    real(real64),   intent(in   ) :: a
    real(real64),   intent(in   ) :: z
    real(real64),   intent(in   ) :: r1

    real(real64) :: t1, t2, temp, temp1

! ======================================================================

    if (l.ne.0) then
       if (ne1 == 0) then
          ne2 = 1
       else
          ne2 = 0
       endif
    else
       t1 = z - one
       t2 = a - one
       temp = t2
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message, 1000) "61, 69, 79, 93"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       if (ine < 0) then
!   pi-
          ne2 = 0
          if (ne1 == 0) then
!    neutron; must have proton partner; final state = nn
             continue
          elseif (ne1 == 1) then
!    proton; final state must be nn or np:
             temp1 = t1/temp
!   Initial partner a neutron: fs = nn
             if (r1 > temp1) ne1 = 0
          endif
       elseif (ine == 0) then
!    pi0
          if (ne1 == 0) then
!    first partner is a neutron
             temp1 = z/temp
             if (r1 > temp1) then
!    2nd neutron in initial pair; fs = nn
                ne2 = 0
             else
!    proton in initial pair; fs = np
                ne2 = 1
             endif
          elseif (ne1 == 1) then
!    first partner is a proton
             temp1 = t1/temp
             if (r1 > temp1) then
!    Initial 2nd partner a neutron; fs = np
                ne2 = 0
             else
!    Initial 2nd partner a 2nd proton; fs = pp
                ne2 = 1
             endif
          endif
       elseif (ine == 1) then
!   pi+
          ne2 = 1
          if (ne1 == 0) then
!    first partner is a neutron
             temp1 = t1/temp
             if (r1 > temp1) then
!    2nd nucleon is a neutron; fs = np
                continue
             else
!    2nd nucleon is a proton; fs = pp
                ne1 = 1
             endif
          elseif (ne1 == 1) then
!    first partner is a proton; initial = np; fs = pp
             continue
          endif
       endif
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'chabs.f90', line(s) ", A)
! ======================================================================
  end subroutine chabs
