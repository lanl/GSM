
  subroutine simulateDecay (gsmObj, gsmRxn, residual, exciton, afMult, &
       & czMult, wam)

! ======================================================================
!
!    Calculates pre-equlibrium and equlibrium particle emission.
!    Various parts removed to separate subroutines to streamline
!    the structure.  10/20/03
!
!    Called by: CEM03
!
!    Calls: AUXL DELTAM EVAPINTERFACE FAM MOLNIX PEQEMT PREQAUX STSTCS kin_energy
!
!    CEM95 written by S. G. Mashnik
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Edited by A. J. Sierk,  LANL  T-2  November, 1997.
!    Modified by AJS, February-March, 1999.
!    Modified by SGM to include reduced masses, 1999
!    Modified by SGM at 06/16/2000
!    Last modification by SGM of 04/27/01, to get the module for CEM2k
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.  
!    Modified by K. K. Gudima, December, 2004.  
!    Edited by AJS, January, 2005.
!    Edited by SGM 07/09/06 to account for the KKG 06/23/06 changes 
!      to use Fermi break-up model in Preco and Evap when A < 13.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Modified by LMK, 2014-2015, for expanison of preeq to heavy ions
! 
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use gsm_params, only : thousand
    use gsm_derived_types, only: nucleus
    use preequilibriumClass, only: &
         & residualState, compoundState, stableState, fermiBreakupState

        ! Excited Residual Nucleus (ERN) Properties:
    implicit none
    class(GSM),         intent(inout) :: gsmObj
    class(GSMReaction), intent(inout) :: gsmRxn
    type(GSMResidual),  intent(in   ) :: residual
    type(excitonData),  intent(in   ) :: exciton
    real(real64),       intent(in   ) :: afMult
    real(real64),       intent(in   ) :: czMult
    real(real64),       intent(inout) :: wam

    type(nucleus) :: gsmNucleus
    type(GSMResults), pointer:: results => NULL()

! ======================================================================

    ! "Reset" fission flags
    results => gsmRxn%results
    results%info%fusion = -1

    ! Setup compound data
    gsmNucleus%ap    = residual%numBaryons
    gsmNucleus%zp    = residual%numProtons
    gsmNucleus%up    = residual%kinEnergy
    gsmNucleus%pnx   = residual%linearMom(1)
    gsmNucleus%pny   = residual%linearMom(2)
    gsmNucleus%pnz   = residual%linearMom(3)
    gsmNucleus%elx   = residual%angularMom(1)
    gsmNucleus%ely   = residual%angularMom(2)
    gsmNucleus%elz   = residual%angularMom(3)
    gsmNucleus%state = 0   ! Nucleus is an excited residual

    ! Simulate preequilbrium physics
    call gsmObj%preequilibriumInterface (gsmRxn, &
         & gsmNucleus, exciton)

    ! Determine "next step" based on the state of the resulting compound nucleus
    if ( gsmNucleus%state == fermiBreakupState ) then
       ! Underwent Fermi Break-up; do nothing
    else if ( gsmNucleus%state == residualState ) then
       ! Warning: The particle bank is full, however the residual nucleus is still undergoing preequilibrium emission (warn user)
       wam = -13.d0
    else if ( gsmNucleus%state == compoundState ) then
       ! Underwent successful equlibration; now evaporate compound

       !   Calculate statistics on average E*, Z, A, L when nucleus enters
       !   statistical decay phase.
 
       !   [KKG 12/01/04; for distribution of A, Z, E*, P, L after preeq. decay:]
       if(results%tallySim) then
          results%info%residual(2)%numBaryons = gsmNucleus%ap
          results%info%residual(2)%numProtons = gsmNucleus%zp
          results%info%residual(2)%energy     = gsmNucleus%ue
          results%info%residual(2)%linearMomentum = &
               & sqrt(gsmNucleus%pnx**2 + gsmNucleus%pny**2 + &
               & gsmNucleus%pnz**2) * thousand
          results%info%residual(2)%angMom = gsmNucleus%ln
          call gsmObj%ststcs (gsmNucleus%ap, gsmNucleus%zp, gsmNucleus%ue, gsmNucleus%ln, &
               & gsmNucleus%bf0, 3)
       end if

       ! Undergo evaporation/fission
       call gsmObj%evaporationInterface (gsmNucleus%ap, gsmNucleus%zp, gsmNucleus%ue, &
            & gsmNucleus%trec, gsmNucleus%pnx, gsmNucleus%pny, gsmNucleus%pnz, &
            & gsmNucleus%ln, gsmNucleus%bf0, &
            & afMult, czMult, gsmRxn)

    else if ( gsmNucleus%state == stableState ) then
       ! Nucleus has low thermal energy and is thus stable; do nothing
    end if

    return
! ======================================================================
  end subroutine simulateDecay
