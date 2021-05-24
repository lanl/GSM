
  subroutine equilibration ( preeqObj, calcVars, exciton, results )

! ======================================================================
!
!    Calculates pre-equlibrium and equlibrium particle emission.
!    Various parts removed to separate subroutines to streamline
!    the structure.  10/20/03
!
!    Called by: simulatePreequilibrium
!
!    Calls: AUXL FAM MOLNIX PEQEMT PREQAUX STSTCS
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
    
    implicit none
    class(Preequilibrium),           intent(inout) :: preeqObj
    type(preequilibriumCalculation), intent(inout) :: calcVars
    type(preeqExcitonData),          intent(inout) :: exciton
    type(preequilibriumResults),     intent(inout) :: results

    integer(int32) :: ipflg, k, lm, nsp
    real(real64)   :: am, pevap, ppre, rnp, rr, temp
    logical preqgo


! ======================================================================

    ! Nucleus has low E* ( < 3 MeV) such that no preequilibrium emission will occur
    if (results%residual%thermEnergy < thr) then
       ! Keep track of nuclei returned with less than 3 MeV of E*:
       ! Residual is stable
       results%residual%state = stableState

       ! Check for photon emission:
       if ( preeqObj%usePhotonEmission ) then
          call preeqObj%photonEmission( results%progenyBnk( results%numProgeny ), &
               & results%residual)
       end if

       return
    endif

    ! ---------------------------------
    !   Beginning of main decay loop:
    ! ---------------------------------
    do k = 1, results%maxProgeny

       ! Undergo Fermi Break-up if E* > E_binding (i.e. small nuclei)
       if ( preeqObj%fbuObj%recommendFermiBreakUp(results%residual%numBaryons) ) then
          call  preeqObj%fermiPreeqInterface (results%residual, results )
          results%residual%state = fermiBreakupState
          return
       endif

       ! Obtain energy information for residual
       ! (mass excess)
       calcVars%dl = &
            & preeqObj%molEnergy%defineEnergy &
            & (calcVars%iz, calcVars%in, 2)
       ! (Pairing gap energies)
       pevap = preeqObj%molEnergy%defineEnergy &
            & (calcVars%iz, calcVars%in, 3)
       ! (Thermal energy of nucleus at rotating ground state)
       results%residual%thermEnergy = results%residual%kinEnergy - pevap - results%residual%rotEnergy


! Determine if excited nucleus has low E* (no preequilibrium/evaporation emissions will occur below this threshold)
       if (results%residual%thermEnergy <= 0.1d0) then
          ! Residual nucleus has reached a stable state
          ! NOTE: No compound nucleus remains at this point; it has reached a ground state (stable)
          results%residual%state = stableState

          ! Check for photon emission:
          if ( preeqObj%usePhotonEmission ) then
             call preeqObj%photonEmission( results%progenyBnk( results%numProgeny ), &
                  & results%residual)
          end if

          return
       endif


!   Set up auxiliary quantities for finding prequilibrium
!   emission widths:
       am = preeqObj%fam (results%residual%numBaryons, results%residual%numProtons, &
            & results%residual%thermEnergy, calcVars, 0)
       calcVars%ac = 0.595d0*am
       call preeqObj%preqaux (calcVars, rr, results)

       !  nsp is the number of excitons at which a compound nucleus is
       !  assumed to be formed.
       nsp = int(   sqrt(  abs( 1.19d0 * am * results%residual%numBaryons * &
            & results%residual%kinEnergy + 0.5d0 )  )   )
       preqgo = (preeqObj%options%excludePreeq == 0) .and. (rr > zro)
       lm = 0

10     continue
!  Restore old logic of sharp transition when n > nsp AJS 07/08/05.
       temp = dble(nsp)
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(preeqObj%io%message,2000) "172"
          call preeqObj%io%print(4, 3, preeqObj%io%message)
       end if
       rnp = dble(exciton%numTotal)/temp

       ! Probability of undergoing pre-equilibrium emission (Eqn. 38 in CEM Manual)
       ! NOTE: The equation used is NOT the same as in the Manual! Either...
       !       (a) the manual needs corrected, or...
       !       (b) the parameters used ought to be refit.
       ppre = one - exp(-0.5d0*((rnp - one)/preeqObj%options%emissionWidth)**2)   ! Gaussian Probability Density Equation
!       ppre = one - exp( -(rnp - one) / (2 * preeqObj%options%emissionWidth**2) )   ! Eqn. 38 of LA-UR-12-01364 (CEM Manual [2012])

       if (exciton%numTotal < nsp .and. preeqObj%rng() < ppre .and. preqgo) then

!  **************** Pre-equilibrium emission *************************
          ipflg = 0
          calcVars%exn = dble(exciton%numTotal)
          
          ! Undergo preequilibrium emission
          call preeqObj%peqemt (calcVars, exciton, lm, ipflg, results)

          !  ipflg = 1 if number of excitons has changed;
          if (ipflg == 1) go to 10
          !  Otherwise, prequilibrium particle of type lm was emitted (progeny stored into bank)

       else
!  ****************** equilibrium emission ***************************

!   Calculate statistics on average E*, Z, A, L when nucleus enters
!   statistical decay phase.
          results%residual%state = compoundState
          return
       endif
    end do

    ! Preequilibrium particle bank is completely filled during equilibration
    write(preeqObj%io%message,1000) k
    call preeqObj%io%print(1, 3, preeqObj%io%message)
    write(preeqObj%io%message,1100)
    call preeqObj%io%print(1, 3, preeqObj%io%message)
    results%simState = 1

    return
! ======================================================================
1000 format("Preequilibrium progeny bank is full (", i3, " progeny created).", &
          & "Cannot simulate further equilibration of excited residual nucleus.")
1100 format("   Results may be suspect.")
2000 format("Divide by zero error prevented in 'equlibration.90', line ", A)
! ======================================================================
  end subroutine equilibration
