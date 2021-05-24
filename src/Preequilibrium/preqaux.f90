
  subroutine preqaux (preeqObj, calcVars, rr, results)

! ======================================================================
!
!   This routine defines several quantities for the particle-
!   emission channels needed for the preequilibrium emission 
!   calculation.  Extracted from the old PRECOF routine.
!
!   CALLED BY: PRECOF
!   CALLS: DELTAM FAM MOLNIX CB_NASA
!
!    Written by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (inclued error protection)
!    Modified by LMK, XCP-3, July 2014, to include new NASA coul. bar. and x_sec
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use preequilibriumParams, only: zro, one, two, four, ato3rd, proton_mass
    use preequilibriumData, only: aj, zj, dlmn, emuredf

! LMK 07/2014
    use coulomb_barrier, only: cb_nasa

    implicit none
    class(Preequilibrium), intent(inout) :: preeqObj
    type(preequilibriumCalculation), intent(inout) :: calcVars
    real(real64),          intent(  out) :: rr  ! Total free energy at the coulomb barrier
    type(preequilibriumResults),     intent(inout) :: results

    integer(int32) :: i, ia, iaa, iaj, iunj, izj
    real(real64)   :: unj    ! Number of neutrons (in real form)
    real(real64)   :: w      ! Total energy in center-of-momentum frame (MeV); LMK 07/2014
    real(real64)   :: t_cm   ! Kinetic energy in center-of-momentum frame (MeV); LMK 07/2014
    real(real64)   :: redcom ! Reduced mass - NOTE: For redcom==0, the decay channel is NOT allowed

! ======================================================================

    rr = zro
    iaa = nint(results%residual%numBaryons)

    do i = 1,preeqObj%options%numPreeqType   ! LMK 07/2012
       calcVars%afj(i) = results%residual%numBaryons - aj(i)
       if (calcVars%afj(i) < 0.5) then
          cycle
       end if
       ia = nint(calcVars%afj(i))
       iaj = nint(aj(i))
       calcVars%afjthr(i) = ato3rd(ia)
       calcVars%zfj(i) = results%residual%numProtons - zj(i)
       redcom = emuredf(iaa,iaj)
! Transform from lab frame to cm frame, LMK 07/2014
       w = sqrt( (aj(i)*proton_mass)**2 + (calcVars%afj(i) * proton_mass)**2 + &
            & 2. * calcVars%afj(i) * proton_mass * (results%residual%kinEnergy + &
            & proton_mass*aj(i)) )
       t_cm = w - (proton_mass*aj(i) + proton_mass*calcVars%afj(i))
! LMK 12/2012, to fix emuredf=zero and redpre div/0 bug
       if (redcom <=  div0Lim) then
          redcom = div0Lim
          write(preeqObj%io%message,1000) &
               & results%residual%numBaryons, results%residual%numProtons
          call preeqObj%io%print(4, 3, preeqObj%io%message)
          write(preeqObj%io%message,1100)
          call preeqObj%io%print(4, 3, preeqObj%io%message)
       end if
       calcVars%redpre(i) = sqrt(abs(one/redcom))
       if (i <= 2) calcVars%redpre(i) = redcom

       ! Check for emission of light fragments ( < He-4 in A or Z)
       if (   calcVars%afj(i) < four .or. &
            & calcVars%zfj(i) < two  .or. &
            & calcVars%zfj(i) >= calcVars%afj(i) ) then
          calcVars%bj(i) = zro
          calcVars%rj(i) = zro
          calcVars%vj(i) = zro
          calcVars%pevapj(i) = zro
          go to 10
       endif

       unj = calcVars%afj(i) - calcVars%zfj(i)
       iunj = nint(unj)
       izj = nint(calcVars%zfj(i))
       calcVars%pevapj(i) = preeqObj%molEnergy%defineEnergy &
            & (izj, iunj, 3)

! "Free" energy in the i'th particle-emission channel:
       calcVars%uej(i) = results%residual%kinEnergy - &
            & calcVars%pevapj(i) - results%residual%rotEnergy

       ! Obtain free energy at coulomb barrier
       if (calcVars%uej(i) <= one) then
          ! No free energy for fragments with low thermal energy
          calcVars%rj(i) = zro
       else
          ! Obtain level density parameter
          calcVars%ami(i) = preeqObj%fam (calcVars%afj(i), &
               & calcVars%zfj(i), calcVars%uej(i), calcVars, i )

          ! Obtain coulomb barrier (use NASA systematics) - LMK [7/2014]
          call cb_nasa(calcVars%zfj(i), calcVars%afj(i), zj(i), &
               & aj(i), t_cm, calcVars%vj(i))

          ! Obtain Binding energy
          calcVars%bj(i) = preeqObj%deltam (calcVars%afj(i), &
               & calcVars%zfj(i)) - ( calcVars%dl - dlmn(i) )

          !   Free "thermal" energy at the Coulomb barrier (corrected for w/ binding energy) [= U_thermal - BE - V_coulomb]
          calcVars%rj(i) = calcVars%uej(i) - &
               & ( calcVars%bj(i) + calcVars%vj(i) )
          calcVars%rj(i) = max (zro, calcVars%rj(i)) ! Ensure valid
       endif

       ! Increase the total free energy at the coulomb barrier
       rr = rr + calcVars%rj(i)

10     continue
    end do
    calcVars%afj(67) = results%residual%numBaryons            !LMK 07/2012
    calcVars%afjthr(67) = calcVars%athrd        !LMK 07/2012

    return
! ======================================================================
1000 format("An unphysical decay channel is possibly being ", &
          & "allowed in the Equilibration (preequilibrium model) of the ", &
          & "excited nucleus (A=", f5.2, ", Z=", f5.2, ").")
1100 format("   Results may be suspect.")
! ======================================================================
  end subroutine preqaux
