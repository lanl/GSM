
  subroutine simulatePreequilibrium (preeqObj, clientResidual, &
       & clientExcitonData, results)

! ======================================================================
!
!    Establishes properties of excited residual nucleus and makes
!    a call to simulate the pre-equilibrium and equilibrium particle
!    emission.
!
!    Calls: auxl, molnix (defineEnergy), and equilibrate
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
    use preequilibriumParams, only : zro, one, two, thr, thousand, ato3rd, emnucb, &
         & emnuct
    
    ! Excited Residual Nucleus (ERN) Properties:
    implicit none
    class(Preequilibrium),       intent(inout) :: preeqObj
    type(residualNucleus),       intent(in   ) :: clientResidual
    type(preeqExcitonData),      intent(in   ) :: clientExcitonData
    type(preequilibriumResults), intent(inout) :: results

    integer(int32) :: ia
    real(real64)   :: delu, e, emx, pevap, un

    type(preequilibriumCalculation) :: calcVars
    type(preeqExcitonData) :: exciton

! ======================================================================

    ! Ensure no errors (class initialized)
    if ( .not. preeqObj%constructed ) then
       ! Class not constructed; print error and return
       write(preeqObj%io%message, 1000)
       call preeqObj%io%print(1, 1, preeqObj%io%message)
       write(preeqObj%io%message, 1999)
       call preeqObj%io%print(1, 1, preeqObj%io%message)
       results%simState = 10
       return
    end if


    ! Ensure data object was constructed:
    if ( .not. preeqObj%preeqData%properlyConstructed() ) then
       write(preeqObj%io%message, 1100)
       call preeqObj%io%print(1, 1, preeqObj%io%message)
       write(preeqObj%io%message, 1999)
       call preeqObj%io%print(1, 1, preeqObj%io%message)
       results%simState = 11
       return
    end if


    ! Ensure progeny array is associated w/ results object:
    if ( .not. results%constructed ) then
       write(preeqObj%io%message, 1200)
       call preeqObj%io%print(1, 1, preeqObj%io%message)
       write(preeqObj%io%message, 1999)
       call preeqObj%io%print(1, 1, preeqObj%io%message)
       results%simState = 12
       return
    end if


    ! Transfer residual nuclei's properties to local variables
    results%residual = clientResidual
    results%residual%state = residualState
    exciton = clientExcitonData

    ! Fill in missing properties of residual, adjust values
    ia = nint(results%residual%numBaryons)
    calcVars%athrd = ato3rd(ia)
    un = results%residual%numBaryons - results%residual%numProtons

    ! Obtain Mass excess
    calcVars%in = nint(un)
    calcVars%iz = nint(results%residual%numProtons)
    calcVars%in = max(1,calcVars%in)
    calcVars%iz = max(1,calcVars%iz)
    emx = preeqObj%molEnergy%defineEnergy &
         & ( calcVars%iz, calcVars%in, 2 )


    ! Obtain rest mass
    if (calcVars%iz > 7 .and. calcVars%in > 7) then
       results%residual%restMass = results%residual%numBaryons * emnucb + emx/thousand
    else
       results%residual%restMass = results%residual%numBaryons * emnuct + emx/thousand
    endif


    ! Total recoil energy and kinetic energy of entire nucleus
    e = sqrt(results%residual%linearMomX**2 + results%residual%linearMomY**2 + &
         & results%residual%linearMomZ**2 + results%residual%restMass**2)
    results%residual%recEnergy = (e - results%residual%restMass)*thousand
    if (e < div0Lim .and. e > -div0Lim) then
       e = div0Lim
       write(preeqObj%io%message,2000) "106-108"
       call preeqObj%io%print(4, 3, preeqObj%io%message)
    end if


    ! Veloctiy
    results%residual%normSpeed(1) = results%residual%linearMomx/e
    results%residual%normSpeed(2) = results%residual%linearMomy/e
    results%residual%normSpeed(3) = results%residual%linearMomz/e

    ! Kinetic Energy
    results%residual%kinEnergy = results%residual%kinEnergy - results%residual%recEnergy

    ! Obtain fission barrier height, angular momentum quantum number, rotational energy, and mass excess, respecitvely
    call preeqObj%preeqData%auxl (results%residual%numBaryons, results%residual%numProtons, &
         & results%residual%angularMom, results%residual%fissBarr, &
         & results%residual%angMomFlag, results%residual%rotEnergy, &
         & delu)

    ! Re-obtain mass excess
    results%residual%kinEnergy = results%residual%kinEnergy + delu
    ! Obtain pairing energy gap
    pevap = preeqObj%molEnergy%defineEnergy &
         & (calcVars%iz, calcVars%in, 3)
    ! Obtain thermal kinetic energy of residual
    results%residual%thermEnergy = results%residual%kinEnergy - pevap - &
         & results%residual%rotEnergy


    ! KKG 12/01/04; for distribution of A, Z, E*, P, L (of residual target nucleus) after the cascade
    results%initResidual = results%residual


    ! Simulate preequilbrium physics
    call preeqObj%equilibration ( calcVars, exciton, results )


    return
! ======================================================================
1000 format("The preequilibrium object has not been constructed.")
1100 format("Data for the preequilibrium simulation was not established.")
1200 format("The preequilibrium results object was not properly setup.")
1999 format("   The equilibration of the excited residual nucleus cannot be ", &
          & "simulated.")
2000 format("Divide by zero error prevented in ", &
          & "'simulatePreequilbrium.90', line ", A)
! ======================================================================

  end subroutine simulatePreequilibrium
