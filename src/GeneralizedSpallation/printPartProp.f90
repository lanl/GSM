
  subroutine printPartProp (gsmObj, proj, targ, results, ncas, intel)

! ======================================================================
!
!   Prints out a table of emitted particle properties for the first nnnp
!   reactions.
!   Basically not used in CEM03; except for debugging.
!
!  Definition of spt:
!                       spt(1,k) = sin(theta.k)
!                       spt(2,k) = cos(theta.k)
!                       spt(3,k) = kinetic energy of particle k (GeV)
!                       spt(4,k) = charge of particle k
!                       spt(5,k) = rest mass of particle k
!
!  Definition of parz:
!                       parz(1,k) = particle type #; 1-9 for n-pi+;
!                                   nint(A + 999*Z) for residual nuclei.
!                       parz(2,k) = kinetic energy of particle k (GeV)
!                       parz(3,k) = theta of particle k
!                       parz(4,k) = phi of particle k
!                       parz(5,k) = index: <100 for cascade particle,
!                        negative for hole;
!                       = 100 for preq.; 200 for coalescence;
!                       = 1000 for thermal; 2000 for evap. from fission.
!                         fragments.
!                       parz(6,k) = electric charge of particle k
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk  LANL  T-2  February, 1996.
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, int64, real64
    use gsm_params, only: thousand, radianToDegree

    implicit none
    class(GSM),          intent(inout) :: gsmObj
    type(GSMProjectile), intent(in   ) :: proj
    type(GSMTarget),     intent(in   ) :: targ
    type(GSMResults),    intent(in   ) :: results
    integer(int64),      intent(in   ) :: ncas
    integer(int64),      intent(in   ) :: intel

    integer(int32) :: k
    real(real64)   :: aa, am, fi, pt, tet, tk, zz

! ======================================================================

    write (16, 100) ncas, intel
    write (16, 200) targ%numBaryons, targ%numProtons, proj%kinEnergy, results%numProgeny
    write (16, 300)
    do k = 1, results%numProgeny
       tet = results%progenyBnk(k)%theta*radianToDegree
       fi = results%progenyBnk(k)%phi*radianToDegree
       tk = results%progenyBnk(k)%kinEnergy*thousand
       am = results%progenyBnk(k)%restMass*thousand
       zz = results%progenyBnk(k)%numProtons
       pt = results%progenyBnk(k)%typeID
       aa = results%progenyBnk(k)%numBaryons
       if (aa > 4.d0 .or. pt > 10.d0) then
          aa = results%progenyBnk(k)%typeID - 999.d0*zz
          am = aa
          pt = aa
       endif
       write (16, 400) pt, zz, tk, tet, fi, results%progenyBnk(k)%prodMech, &
            & results%progenyBnk(k)%sinTheta, results%progenyBnk(k)%cosTheta, am
    end do
    return

! ======================================================================
100 format (/1x,'ncas =',i7,', intel =',i7)
200 format (1x,'After the cascade, ', &
         &       'At = ',f4.0,', Zt =',f4.0,', Ut =',f7.2,', ktot =',i3)
300 format (2x,'Par  Q',4x,'T(MeV)',1x,'theta',3x,'phi',3x, &
         &       'ngen',5x,'st',6x,'ct',5x,'M(MeV)')
400 format (1x,f4.0,1x,f3.0,1x,f7.2,2x,f5.1,2x,f5.1,2x,f5.0,2x,f6.4, &
         &       2x,f6.4,1x,f7.1)
! ======================================================================
  end subroutine printPartProp
