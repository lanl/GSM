
  subroutine renorm (gsmObj, results, ibadren)

! ======================================================================
!
!    Checks the law of conservation for the baryon and charge numbers
!    after the cascade stage of a reaction and renormalizes the
!    excitation energy and momentum of the excited residual nucleus
!    imposed by momentum-energy conservation for each simulated
!    inelastic event using real binding energies and masses.
!
!   CALLS: DELTAM
!
!    Written by S. G. Mashnik, LANL, T-2, 1998
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by CMJ, LANL XCP-3 (06/2018) to include A/Z conservation w/
!         LAQGSM INC (yields incorrect results when used w/ LAQGSM INC)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro, two, ibarpar, parmass

    implicit none
    class(GSM),       intent(inout) :: gsmObj
    type(GSMResults), intent(inout), target :: results
    integer(int32),   intent(  out) :: ibadren

    integer(int32) :: ibaryon, icharge, indx, innbar, innch, k
    real(real64)   :: aat, bez, dlm, emmass, eres, etarget, pnuc2, &
         & pxem, pyem, pz0, pzem, rec, tem, uae

    real(real64), dimension(results%numProgeny) :: pxp, pyp, pzp, pp

    type(GSMProjectile), pointer :: proj  => NULL()
    type(GSMTarget),     pointer :: targ  => NULL()
    type(GSMResidual),   pointer :: resid => NULL()

! ======================================================================

    ! Initialize pointers and return type
    ibadren = 0
    proj  => results%initialProj
    targ  => results%initialTarg
    resid => results%targRes

    do k = 1, results%numProgeny
       ! Assigning particle momenta (total, in x, y, and z)
       indx = nint(results%progenyBnk(k)%typeID)
       pp(k)  = sqrt(abs( results%progenyBnk(k)%kinEnergy * &
            & (results%progenyBnk(k)%kinEnergy + two*parmass(indx)) ))
       pxp(k) = pp(k)*results%progenyBnk(k)%sinTheta*cos(results%progenyBnk(k)%phi)
       pyp(k) = pp(k)*results%progenyBnk(k)%sinTheta*sin(results%progenyBnk(k)%phi)
       pzp(k) = pp(k)*results%progenyBnk(k)%cosTheta
    end do
    icharge = 0
    ibaryon = 0
    emmass = zro
    pxem   = zro
    pyem   = zro
    pzem   = zro
    tem    = zro
    pz0    = sqrt(abs(proj%kinEnergy*(proj%kinEnergy + two*proj%restMass)))
    do k = 1, results%numProgeny
       indx = nint(results%progenyBnk(k)%typeID)
       icharge = icharge + nint(results%progenyBnk(k)%numProtons)
       ibaryon = ibaryon + ibarpar(indx)
       emmass = emmass + parmass(indx)
       pxem   = pxem + pxp(k) ! total x momentum
       pyem   = pyem + pyp(k) ! total y momentum
       pzem   = pzem + pzp(k) ! total z momentum
       tem    = tem + results%progenyBnk(k)%kinEnergy
    end do
    icharge = icharge + nint(resid%numProtons) ! total charge (residuals + target)
    ibaryon = ibaryon + nint(resid%numBaryons) ! total mass (residuals + target)
    innbar  = proj%numBaryons + nint(targ%numBaryons) ! initial mass (total nucleons)
    innch   = proj%numProtons + nint(targ%numProtons) ! initial charge (total protons)
    if (innbar.ne.ibaryon .or. innch.ne.icharge) then
! Conservation failure
       if (innbar.ne.ibaryon .and. innch.ne.icharge) then
          write(gsmObj%io%message, 1000)
       elseif (innbar.ne.ibaryon ) then
          write(gsmObj%io%message, 1030)
       elseif (innch.ne.icharge) then
          write(gsmObj%io%message, 1060)
       end if
       call gsmObj%io%print(3, 3, gsmObj%io%message)
       write(gsmObj%io%message, 1999)
       call gsmObj%io%print(3, 3, gsmObj%io%message)
       ibadren = 1
       return
    endif
!   Calculation of the target mass in GeV:
!   Parametrization used by Peter Moller (1997) for the total atomic
!   binding energy (noted as Be(Z) in formula (1) by Audi & Wapstra,
!   NP, A565 (1993) 1-65); We use here formula (4) and all constants
!   as in:  P. Moller, J.R. Nix, and K.-L. Kratz, At. Data Nuc. Tab.,
!           66 (1997) 131
!
    bez = 0.00001433d0*targ%numProtons**2.39d0
    dlm = gsmObj%genData%molnixEnergies%massExcess (targ%numBaryons, targ%numProtons)
    uae = 931.5014d0
    aat = uae*targ%numBaryons + dlm
    etarget = (aat - targ%numProtons*0.51099906d0 + bez)/1000.d0

    !   Calculation of the residual nucleus mass in GeV:
    bez = 0.00001433d0*resid%numProtons**2.39d0
    dlm = gsmObj%genData%molnixEnergies%massExcess (resid%numBaryons, resid%numProtons)
    aat = uae*resid%numBaryons + dlm
    eres = (aat - resid%numProtons*0.51099906d0 + bez)/1000.d0
    resid%linearMom(1) = -pxem
    resid%linearMom(2) = -pyem
    resid%linearMom(3) = pz0 - pzem
    pnuc2 = resid%linearMom(1) * resid%linearMom(1) + &
         &  resid%linearMom(2) * resid%linearMom(2) + &
         &  resid%linearMom(3) * resid%linearMom(3)
    rec = sqrt(pnuc2 + eres*eres) - eres
    resid%kinEnergy = proj%kinEnergy + proj%restMass + etarget - emmass - eres - tem - rec
    if (resid%kinEnergy <= zro) then
       write(gsmObj%io%message, 1200) resid%kinEnergy
       call gsmObj%io%print(3, 3, gsmObj%io%message)
       write(gsmObj%io%message, 1999)
       call gsmObj%io%print(3, 3, gsmObj%io%message)
       resid%kinEnergy = zro
       ibadren = 1
    endif

    return
! ======================================================================
1000 format("Mass and charge conservation failure from ", &
         & "INC process.")
1030 format("Mass conservation failure from ", &
         & "INC process.")
1060 format("Charge conservation failure from ", &
         & "INC process.")
1200 format("A negative particle energy from INC process was found", &
         & "(E_k=", ES15.4, ").")
1999 format("   Restarting event...")
! ======================================================================
  end subroutine renorm
