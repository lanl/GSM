
  subroutine gamma (evapObj, compound, uran, beta, calcVars)

! ======================================================================
!
! GAMMA
!  Calculate decay width for each particle emission channel.
!
!  Called by: STDCAY
!
!  Calls:
! ====================================================================
! <Subroutine>
!   eye    : Calculate I0,I1,I2,I3
! ====================================================================
! <Function>
!  dost    : Calculate kp, cp, or k_alpha
!  paire   : Calculate pairing energy (NOT USED)
! ====================================================================
! <variables>
!     a    :   residual mass before emission                     (IN)
!     z    :   charge number of res nucleus before emission      (IN)
!     u    :   excitation energy of nucleus before emission      (IN)
!   uran   :   random number                                     (IN)
!     r    :   decay width                                       (OUT)
!     rr   :   decay width enchancement factor by excited-state
!                            particle emission                   (OUT)
!   ifa    :   mass # of emitted particle                        (OUT)
!   ifz    :   charge # of emitted particle                      (OUT)
!  omega   :   spin of emittor                                   (OUT)
!   gj     :   gj in eq.(39)  (~ spin multiplicity factor*mu)    (OUT)
!   q      :   Q-value                                           (OUT)
!  delta   :   pairing energy                                    (OUT)
!  smalla  :   level density parameter                           (OUT)
!   v      :   Coulomb barrier                                   (OUT)
!   exm    :   mass of an excited state [MeV]                    (OUT)
!  spin    :   spin of an excited state                          (OUT)
!  width   :   lifetime of an excited state [MeV]                (OUT)
!   gamn   :   Decay width for neutron emission                  (OUT)
!   an     :   Level density parameter for neutron emission      (OUT)
! beta,alp :   Inverse cross section params for neutron emission (OUT)
!   s      : NOT USED!  Removed 1/24/05 AJS.
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Modified by A. J. Sierk, LANL T-16, January, 2005.
!    Modified by A. J. Sierk, LANL T-16, March, 2005.
!    Modified by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (Evap class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, hlf, one, two, thr, hbarc, amu, pi
    use evaporationFissionData, only: exm, spin, width, omega, ifa, &
         & ifz, pz, pn, paire0, maxSize, aMax, omfct, alpp, betap

    implicit none
    class(Evaporation),  intent(inout) :: evapObj
    type(evapCompound),  intent(in   ) :: compound
    real(real64),        intent(in   ) :: uran
    real(real64),        intent(  out) :: beta
    type(evaporationCalculation), intent(inout) :: calcVars

    integer(int32) :: i, ia, iaa, &
         & isdum, iz, izz, j, k, nn, nn0
    real(real64)   :: aa, alp, aq, bett, del, emred, gq, q1, q2, qq, &
         & rhop, rmass, rq, rrq, rrqg, temp, zz

! ======================================================================

    real(real64), parameter :: cst0 = amu / ( pi * (hbarc**2) )

! ======================================================================

!   Parent nucleus:
    ia = nint(compound%numBaryons)
    iz = nint(compound%numProtons)
    q1 = evapObj%energy (iz, ia)
    if (q1 >= 1.d10) then
       write(evapObj%io%message, 2000) ia, iz
       call evapObj%io%print(2, 2, evapObj%io%message)
       q1 = 1.d10
    endif

    if (evapObj%options%inverseParameter == zro) then
!   Set cj and kj for Dostrovsky's parameter set
       calcVars%couk(2) = dost (1, compound%numProtons-dble(ifz(2)))
       calcVars%couk(3) = calcVars%couk(2) + 6.d-2
       calcVars%couk(4) = calcVars%couk(2) + 1.2d-1
       calcVars%couk(6) = dost (2, compound%numProtons-dble(ifz(6)))
       calcVars%couk(5) = calcVars%couk(6) - 6.d-2
       calcVars%couc(2) = dost (3, compound%numProtons-dble(ifz(2)))
       calcVars%couc(3) = calcVars%couc(2)*hlf
       calcVars%couc(4) = calcVars%couc(2)/thr
!   ck is given by expression: ck = C1*log(Z) + C2*log(A) + C3
    elseif (evapObj%options%inverseParameter < zro) then                     
       do j = 2,6                                 
          calcVars%couk(j) = ckcal (j, compound%numProtons - dble(ifz(j)), compound%numBaryons - dble(ifa(j)))   
       end do
    endif
    nn0 = ia - iz
    calcVars%smalla0 = geta (evapObj%data%alev(), compound%kinEnergy, &
         & iz, nn0, isdum)
!  Pairing energy for the  parent nucleus:
    del = paire0(iz,ia-iz)
!  Rho_i (level density of parent))
    rhop = evapObj%rho (ia, iz, compound%kinEnergy, del, calcVars)
    if (rhop <= zro) then
       write(evapObj%io%message, 2100) rhop
       call evapObj%io%print(2, 3, evapObj%io%message)
    end if

!     Start calculating Gamma for each particle type:
    do j = 1, evapObj%options%numEvapType

!  Initialization:
       calcVars%r(j) = zro
       calcVars%rr(j) = one
       calcVars%gj(j) = zro
       if (ifa(j) == 0) go to 20
! ----------------------------------------------------------------------
!    Reduce calculation time:
! ----------------------------------------------------------------------
       if (evapObj%options%redCompTime.ne.0) then
          if (j > 6) then
             if (ia > 40 .and. uran < 0.95d0) go to 20
             if (ia > 30 .and. ia <= 40 .and. uran < 0.93d0) go to 20
             if (ia > 20 .and. ia <= 30 .and. uran < 0.7d0) go to 20
          endif
       endif
!  Daughter nucleus mass and charge:
       iaa = ia - ifa(j)
       aa = dble(iaa)
       izz = iz - ifz(j)
       zz = dble(izz)
       nn = iaa - izz
! ----------------------------------------------------------------------
!    Check the residual nuclei after the emission for reasonable Z, A
!    and avoid double counting:
! ----------------------------------------------------------------------
       if (iaa <= 0 .or. izz <= 0 .or. iaa < izz .or. nn <= 0) &
            & go to 20
       if (iaa < ifa(j) .or. izz < ifz(j)) then
          do k = 1,j-1
             if (iaa == ifa(k) .and. izz == ifz(k)) go to 20
          end do
       endif

!   Q-value
       q2 = evapObj%energy (izz, iaa)
       if (q2 >= 1.d10) then
          write(evapObj%io%message,2000) iaa, izz
          call evapObj%io%print(2, 2, evapObj%io%message)
          q2 = 1.d10
       endif
       calcVars%q(j) = q2 - q1 + evapObj%energy (ifz(j), ifa(j))

!   Pairing energy of daughter:
       calcVars%delta(j) = paire0(izz, nn)

!   ............................................................
!  Coulomb potential
       calcVars%v(j) = evapObj%vcoul (zz, aa, ifz(j), ifa(j), calcVars%couk(j), j)
!    ..Set Coulomb potential to zero if Q-value < 0, for light residual
!    ..... e.g. 8Be, 9B
       if (calcVars%q(j) <= zro .and. iaa <= 20) calcVars%v(j) = zro

!   Check whether the emission is energetically possible or not:
       if (compound%kinEnergy - calcVars%q(j) - calcVars%v(j) <= zro) go to 20

!   alpha and beta for neutrons: the precise parameter set:
       if (j == 1) then
          if (evapObj%options%inverseParameter <= zro) then
             alp = alpp(iaa)
             beta = betap(iaa)
          else
!   Simple model:
             alp = one
             beta = zro
          endif
          bett = beta
       else   
          bett = -calcVars%v(j)
       endif

       if (evapObj%options%inverseParameter == -one .and. j > 6) then
!    Exact calculation for subthreshold reaction
          if (calcVars%r(1) + calcVars%r(2) + calcVars%r(3) + calcVars%r(4) + &
               & calcVars%r(5) + calcVars%r(6).ne.zro) call evapObj%calr &
               & (aa, zz, compound%kinEnergy, calcVars%q(j), calcVars%v(j), calcVars%delta(j), &
               & calcVars%r(j), calcVars)
       else
!   Level density parameter of daughter:
          calcVars%smalla(j) = geta ( evapObj%data%alev(), &
               & (compound%kinEnergy-calcVars%q(j)-calcVars%v(j)), &
               & izz, nn, isdum)
          call evapObj%eye (iaa, izz, nn, compound%kinEnergy, calcVars%q(j), calcVars%v(j), &
               & calcVars%delta(j), calcVars%smalla(j), bett, calcVars%r(j), calcVars)
       endif
!   ............................................................

!   (2Sj + 1)*muj*alpha
       temp = compound%numBaryons
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '308'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       emred = dble(ifa(j))*aa/temp
       calcVars%gj(j) = omfct(j)*emred
       if (j == 1) then
          calcVars%gj(j) = calcVars%gj(j)*alp
       else
          calcVars%gj(j) = calcVars%gj(j)*(one + calcVars%couc(j))
       endif

!  Geometric cross section... sigma R
       rmass = evapObj%rb (aa, ifa(j), zz, ifz(j), j)
       calcVars%gj(j) = rmass*rmass*calcVars%gj(j)
       if (iz==7 .and. ia==13 .and. j==2) then
!          print *, 'Geom. radius: ', rmass
!          print *, 'cp = ', calcVars%couc(2)
!          print *, 'Fragment energy: ', q(2), 'GEM2 cross section: ', gj(2)
       end if
!  Decay width:
       calcVars%r(j) = calcVars%gj(j)*calcVars%r(j)
       calcVars%r(j) = max(calcVars%r(j), zro)

!  Decay width calculation from excited states:
       if (calcVars%r(j) == zro .or. j <= 6) go to 20

!   Move outside the loop over j; AJS (02/11/09)
!  Pairing energy for the  parent nucleus:
!  Rho_i (level density of parent))
!   End of j-independent code fragment.  (02/11/09)
       rrq = calcVars%r(j)
       do i = 1,200
          if (width(j-6,i) == zro) go to 10
!   Q-value for an excited state:
          qq = calcVars%q(j) + exm(j-6,i)
!   Level density for an excited state:
          aq = geta (evapObj%data%alev(), compound%kinEnergy - qq - calcVars%v(j), izz, nn, isdum)
!   Integral part of gamma for an excited state:
          call evapObj%eye (iaa, izz, nn, compound%kinEnergy, qq, calcVars%v(j), &
               & calcVars%delta(j), aq, bett, rq, calcVars)
!   (2Sj + 1)*mu_j*alpha*sigma_R for an excited state:
          gq = (two*spin(j-6,i) + one)*emred
          gq = rmass*rmass*gq
!   gamma for an excited state [MeV]:
          temp = rhop
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(evapObj%io%message,1000) '354'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          rrqg = gq*rq*cst0/temp
!   Reject an excited state if the decay width [MeV] < level width [MeV]:
          if (width(j-6,i) >= rrqg) go to 10
          rrq = rrq + gq*rq
       end do
10     temp = rrq
       if (temp < 1.0d-90 .and. temp > -1.0d-90) then
          print *, temp
          temp = 1.0d-90
          write(evapObj%io%message,1000) '370'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       calcVars%rr(j) = calcVars%r(j)/temp
       calcVars%r(j) = rrq
20     continue
    end do
!  End of j (ejectile index) DO loop ^
    calcVars%sigma = zro
    do j = 1, evapObj%options%numEvapType
       calcVars%sigma = calcVars%sigma + calcVars%r(j)
    end do

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'gamma.f90', line ", A)
2000 format("Mass calculation error occurred in function ", &
            & "'energy' for the nucleus with ia = ", i3, &
            & " and iz = ", i3, ".")
2100 format("Negative level density parameter (", f8.2, &
          & ") ocurred in function 'rho'.")
! ======================================================================
  end subroutine gamma
