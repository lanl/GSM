
  subroutine resdist (gsmObj, event, zbeg, abeg, fisevent)

! ======================================================================
!
!   Accumulate the distribution of residual nuclides in the appropriate
!   arrays for:
!      rdis(1,.,.):  after the cascade;
!      rdis(2,.,.):  after preeq. decay;  (= sum of (3,..) and (4,..))
!      rdis(3,.,.):  nuclei before evaporation which will not fission;
!      rdis(4,.,.):  nuclei before evaporation which will fission;
!      rdis(5,.,.):  nuclei just before fissioning.
!
!      rdis(.,1,.):  A distribution;
!      rdis(.,2,.):  Z distribution;
!      rdis(.,3,.):  E* distribution;
!      rdis(.,4,.):  |P| distribution;
!      rdis(.,5,.):  L distribution.
!
!   Written by K. K. Gudima, December, 2004.
!   Edited by A. J. Sierk, LANL T-16, January, 2005.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: one
    use EventDataMod, only: EventData

    implicit none
    class(GSM),        intent(inout) :: gsmObj
    class(EventData),  intent(in   ) :: event
    real(real64),      intent(in   ) :: zbeg
    real(real64),      intent(in   ) :: abeg
    logical,           intent(in   ) :: fisevent

    integer(int32) :: ia, iam, iex, ipm, iz
    real(real64)   :: da, dz

! ======================================================================

    real(real64) :: rdis, dex, dpm
    common /resdis/ rdis(5,5,250), dex, dpm

! ======================================================================
!  KKG  12/02/04
!   Accumulate residual nuclei information:

!   Nuclei after cascade:
    da = abeg - event%residual(1)%numBaryons
    dz = zbeg - event%residual(1)%numProtons
    ia = nint(da) + 1
    ia = min(nint(abeg), ia)
    ia = min(ia, 247)
    iz = nint(dz) + 1
    iz = min(nint(zbeg) + 1, iz)
    iz = min(iz, 247)
    if ( dex < div0Lim .and. dex > -div0Lim ) then
       dex = div0Lim
       write(gsmObj%io%message,1000) "62"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    iex = int(event%residual(1)%energy/dex) + 1
    if ( dpm < div0Lim .and. dpm > -div0Lim ) then
       dpm = div0Lim
       write(gsmObj%io%message,1000) "62"
       call gsmObj%io%print(4, 3, gsmObj%io%message) 
   end if
    ipm = int(event%residual(1)%linearMomentum/dpm) + 1
    iam = event%residual(1)%angMom     + 1
    if (ia <= 247) then
       rdis(1,1,ia) = rdis(1,1,ia)  + one
       rdis(1,1,248)= rdis(1,1,248) + event%residual(1)%numBaryons
       rdis(1,1,249)= rdis(1,1,249) + event%residual(1)%numBaryons**2
       rdis(1,1,250)= rdis(1,1,250) + one
    endif
    if (iz <= 247) then
       rdis(1,2,iz) = rdis(1,2,iz)  + one
       rdis(1,2,248)= rdis(1,2,248) + event%residual(1)%numProtons
       rdis(1,2,249)= rdis(1,2,249) + event%residual(1)%numProtons**2
       rdis(1,2,250)= rdis(1,2,250) + one
    endif
    if (iex > 0 .and. iex <= 247) then
       rdis(1,3,iex) = rdis(1,3,iex) + one
       rdis(1,3,248) = rdis(1,3,248) + event%residual(1)%energy
       rdis(1,3,249) = rdis(1,3,249) + event%residual(1)%energy**2
       rdis(1,3,250) = rdis(1,3,250) + one
    endif
    if (ipm > 0 .and. ipm <= 247) then
       rdis(1,4,ipm) = rdis(1,4,ipm) + one
       rdis(1,4,248) = rdis(1,4,248) + event%residual(1)%linearMomentum
       rdis(1,4,249) = rdis(1,4,249) + event%residual(1)%linearMomentum**2
       rdis(1,4,250) = rdis(1,4,250) + one
    endif
    if (iam > 0 .and. iam <= 247) then
       rdis(1,5,iam) = rdis(1,5,iam) + one
       rdis(1,5,248) = rdis(1,5,248) + dble(event%residual(1)%angMom)
       rdis(1,5,249) = rdis(1,5,249) + dble(event%residual(1)%angMom)**2
       rdis(1,5,250) = rdis(1,5,250) + one
    endif
!   Nuclei after preequilibrium:
    da = abeg - event%residual(2)%numBaryons
    dz = zbeg - event%residual(2)%numProtons
    ia = nint(da) + 1
    ia = min(nint(abeg), ia)
    ia = min(ia, 247)
    iz = nint(dz) + 1
    iz = min(nint(zbeg) + 1, iz)
    iz = min(iz, 247)
    iex = int(event%residual(2)%energy/dex) + 1
    ipm = int(event%residual(2)%linearMomentum/dpm) + 1
    iam = event%residual(2)%angMom     + 1
    if (ia <= 247) then
       rdis(2,1,ia) = rdis(2,1,ia)  + one
       rdis(2,1,248)= rdis(2,1,248) + event%residual(2)%numBaryons
       rdis(2,1,249)= rdis(2,1,249) + event%residual(2)%numBaryons**2
       rdis(2,1,250)= rdis(2,1,250) + one
    endif
    if (iz <= 247) then
       rdis(2,2,iz) = rdis(2,2,iz)  + one
       rdis(2,2,248)= rdis(2,2,248) + event%residual(2)%numProtons
       rdis(2,2,249)= rdis(2,2,249) + event%residual(2)%numProtons**2
       rdis(2,2,250)= rdis(2,2,250) + one
    endif
    if (iex > 0 .and. iex <= 247) then
       rdis(2,3,iex) = rdis(2,3,iex) + one
       rdis(2,3,248) = rdis(2,3,248) + event%residual(2)%energy
       rdis(2,3,249) = rdis(2,3,249) + event%residual(2)%energy**2
       rdis(2,3,250) = rdis(2,3,250) + one
    endif
    if (ipm > 0 .and. ipm <= 247) then
       rdis(2,4,ipm) = rdis(2,4,ipm) + one
       rdis(2,4,248) = rdis(2,4,248) + event%residual(2)%linearMomentum
       rdis(2,4,249) = rdis(2,4,249) + event%residual(2)%linearMomentum**2
       rdis(2,4,250) = rdis(2,4,250) + one
    endif
    if (iam > 0 .and. iam <= 247) then
       rdis(2,5,iam) = rdis(2,5,iam) + one
       rdis(2,5,248) = rdis(2,5,248) + dble(event%residual(2)%angMom)
       rdis(2,5,249) = rdis(2,5,249) + dble(event%residual(2)%angMom)**2
       rdis(2,5,250) = rdis(2,5,250) + one
    endif
    if (.not. fisevent) then
!   Nuclei after preeq. which will NOT fission:
       if (ia <= 247) then
          rdis(3,1,ia) = rdis(3,1,ia)  + one
          rdis(3,1,248)= rdis(3,1,248) + event%residual(2)%numBaryons
          rdis(3,1,249)= rdis(3,1,249) + event%residual(2)%numBaryons**2
          rdis(3,1,250)= rdis(3,1,250) + one
       endif
       if (iz <= 247) then
          rdis(3,2,iz) = rdis(3,2,iz)  + one
          rdis(3,2,248)= rdis(3,2,248) + event%residual(2)%numProtons
          rdis(3,2,249)= rdis(3,2,249) + event%residual(2)%numProtons**2
          rdis(3,2,250)= rdis(3,2,250) + one
       endif
       if (iex > 0 .and. iex <= 247) then
          rdis(3,3,iex) = rdis(3,3,iex) + one
          rdis(3,3,248) = rdis(3,3,248) + event%residual(2)%energy
          rdis(3,3,249) = rdis(3,3,249) + event%residual(2)%energy**2
          rdis(3,3,250) = rdis(3,3,250) + one
       endif
       if (ipm > 0 .and. ipm <= 247) then
          rdis(3,4,ipm) = rdis(3,4,ipm) + one
          rdis(3,4,248) = rdis(3,4,248) + event%residual(2)%linearMomentum
          rdis(3,4,249) = rdis(3,4,249) + event%residual(2)%linearMomentum**2
          rdis(3,4,250) = rdis(3,4,250) + one
       endif
       if (iam > 0 .and. iam <= 247) then
          rdis(3,5,iam) = rdis(3,5,iam) + one
          rdis(3,5,248) = rdis(3,5,248) + dble(event%residual(2)%angMom)
          rdis(3,5,249) = rdis(3,5,249) + dble(event%residual(2)%angMom)**2
          rdis(3,5,250) = rdis(3,5,250) + one
       endif
    else
!   Nuclei after preeq. which will fission:
       if (ia <= 247) then
          rdis(4,1,ia) = rdis(4,1,ia)  + one
          rdis(4,1,248)= rdis(4,1,248) + event%residual(2)%numBaryons
          rdis(4,1,249)= rdis(4,1,249) + event%residual(2)%numBaryons**2
          rdis(4,1,250)= rdis(4,1,250) + one
       endif
       if (iz <= 247) then
          rdis(4,2,iz) = rdis(4,2,iz)  + one
          rdis(4,2,248)= rdis(4,2,248) + event%residual(2)%numProtons
          rdis(4,2,249)= rdis(4,2,249) + event%residual(2)%numProtons**2
          rdis(4,2,250)= rdis(4,2,250) + one
       endif
       if (iex > 0 .and. iex <= 247) then
          rdis(4,3,iex) = rdis(4,3,iex) + one
          rdis(4,3,248) = rdis(4,3,248) + event%residual(2)%energy
          rdis(4,3,249) = rdis(4,3,249) + event%residual(2)%energy**2
          rdis(4,3,250) = rdis(4,3,250) + one
       endif
       if (ipm > 0 .and. ipm <= 247) then
          rdis(4,4,ipm) = rdis(4,4,ipm) + one
          rdis(4,4,248) = rdis(4,4,248) + event%residual(2)%linearMomentum
          rdis(4,4,249) = rdis(4,4,249) + event%residual(2)%linearMomentum**2
          rdis(4,4,250) = rdis(4,4,250) + one
       endif
       if (iam > 0 .and. iam <= 247) then
          rdis(4,5,iam) = rdis(4,5,iam) + one
          rdis(4,5,248) = rdis(4,5,248) + dble(event%residual(2)%angMom)
          rdis(4,5,249) = rdis(4,5,249) + dble(event%residual(2)%angMom)**2
          rdis(4,5,250) = rdis(4,5,250) + one
       endif
!  Distributions of nuclei just before they fission;
!  NOTE: in GEM2 angular momentum is not changed!
       da = abeg - event%compound%numBaryons
       dz = zbeg - event%compound%numProtons
       ia = nint(da) + 1
       ia = min(nint(abeg), ia)
       ia = min(ia, 247)
       iz = nint(dz) + 1
       iz = min(nint(zbeg) + 1, iz)
       iz = min(iz, 247)
       iex = int(event%compound%kinEnergy/dex) + 1
       ipm = int(event%compound%linearMomTot/dpm) + 1
       iam = event%residual(2)%angMom  + 1
       if (ia <= 247) then
          rdis(5,1,ia) = rdis(5,1,ia)  + one
          rdis(5,1,248)= rdis(5,1,248) + event%compound%numBaryons
          rdis(5,1,249)= rdis(5,1,249) + event%compound%numBaryons**2
          rdis(5,1,250)= rdis(5,1,250) + one
       endif
       if (iz <= 247) then
          rdis(5,2,iz) = rdis(5,2,iz)  + one
          rdis(5,2,248)= rdis(5,2,248) + event%compound%numProtons
          rdis(5,2,249)= rdis(5,2,249) + event%compound%numProtons
          rdis(5,2,250)= rdis(5,2,250) + one
       endif
       if (iex > 0 .and. iex <= 247) then
          rdis(5,3,iex) = rdis(5,3,iex) + one
          rdis(5,3,248) = rdis(5,3,248) + event%compound%kinEnergy
          rdis(5,3,249) = rdis(5,3,249) + event%compound%kinEnergy**2
          rdis(5,3,250) = rdis(5,3,250) + one
       endif
       if (ipm > 0 .and. ipm <= 247) then
          rdis(5,4,ipm) = rdis(5,4,ipm) + one
          rdis(5,4,248) = rdis(5,4,248) + event%compound%linearMomTot
          rdis(5,4,249) = rdis(5,4,249) + event%compound%linearMomTot**2
          rdis(5,4,250) = rdis(5,4,250) + one
       endif
       if (iam > 0 .and. iam <= 247) then
          rdis(5,5,iam) = rdis(5,5,iam) + one
          rdis(5,5,248) = rdis(5,5,248) + dble(event%residual(2)%angMom)
          rdis(5,5,249) = rdis(5,5,249) + dble(event%residual(2)%angMom)**2
          rdis(5,5,250) = rdis(5,5,250) + one
       endif
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'resdist.f90', line(s) ", A)
! ======================================================================
  end subroutine resdist
