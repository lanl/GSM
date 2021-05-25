
  function fprob (evapObj, compound, beta, calcVars)

! ======================================================================
!
! ***************This routine was originally in LAHET code**************
!  FPROB
!    Calculate fission probability
! ======================================================================
! <variables>
!     a    :   mass of nucleus before fission                     (IN)
!     z    :   charge  of nucleus before fission                  (IN)
!     e    :   excitation energy of nucleus before fission        (IN)
!  sigma   :   total decay width calculated in sub. gamma         (IN)
!   an     :   level density parameter for neutron emission
!                                calculated in sub. gamma         (IN)
!  frob    :   fission probability                               (OUT)
!
!    Function to compute the fission probability.
!     for z<90 uses ..........
!     uses statistical model fits.
!
!    Modified by SGM to fit a_f/a_n, 12/04/01
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Modified by KKG, 11/24/04
!    Edited by A. J. Sierk, January, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams,             only: zro, hlf, one, two, ato3rd
    use evaporationFissionData,        only: inn, iiz, paire0

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    type(evapCompound), intent(in   ) :: compound
    real(real64),       intent(in   ) :: beta
    type(evaporationCalculation), intent(inout) :: calcVars
    real(real64)                      :: fprob

    integer(int32) :: ia, in, iz
    real(real64)   :: a1thrd, af, agoes, an, an1, an2, e1, e2, ef, &
         & epsav, exps, eye0, eye1, eye10, eye2, gamnf, i0, i1, s2, &
         & se, ss, temp, u, x

! ======================================================================

    integer(int32),                 parameter :: lnfis = 18
    real(real64), dimension(lnfis), parameter :: slope = &
         & [ 0.23d0,    0.233d0,   0.12225d0, 0.14727d0, 0.13559d0, &
         &   0.15735d0, 0.16597d0, 0.17589d0, 0.18018d0, 0.19568d0, &
         &   0.16313d0, 0.17123d0, 0.1700d0,  0.1700d0,  0.1700d0,  &
         &   0.1700d0,  0.1700d0,  0.1700d0  ]
    real(real64), dimension(lnfis), parameter :: anort = &
         & [    219.4d0,  226.9d0,  229.75d0, 234.04d0, 238.88d0, &
         &      241.34d0, 243.04d0, 245.52d0, 246.84d0, 250.18d0, &
         &      254.0d0,  257.8d0,  261.3d0,  264.8d0,  268.3d0, &
         &      271.8d0,  275.3d0,  278.8d0   ]
    real(real64), parameter :: a1    =   0.2185024d0
    real(real64), parameter :: a2    =  16.70314d0
    real(real64), parameter :: a3    = 321.175d0
    real(real64), parameter :: const =   0.3518099d0
    real(real64), parameter :: c1    =   1.089257d0
    real(real64), parameter :: c2    =   0.01097896d0
    real(real64), parameter :: c3    =  31.08551d0

! ======================================================================

!   Initialization
    fprob = zro
    iz = nint(compound%numProtons)
    ia = nint(compound%numBaryons)
    in = nint(compound%numBaryons - compound%numProtons)
    u = compound%kinEnergy

    if (iz > 88) then
!   High-Z fission probability for Z>= 89 and <100 ........
!   Use the systematics of vandenbosch & Huizenga with a ball park
!    observation that fission probability drops for most nuclei at 6 MeV
       iz = iz - 88
       if (iz > lnfis .or. u < 6.d0) return
       gamnf = slope(iz)*(compound%numBaryons - anort(iz))
       gamnf = 1.d1**gamnf
       gamnf = gamnf*compound%czMultiplier

!  Here, compound%czMultiplier is an additional fitting parameter in Eq. (2)
!  (to not change for the moment Atchison's parameter C(Z) defined
!  in this function by the operator: data slope /.../,
!  for fissioning nuclei with  Z>88; compound%czMultiplier = C(Z)[CEM]/C(Z)[RAL]
!  SGM, 01/03/2002
       agoes = one
    else
!   Calculate -Qn + the pairing energy + the shell correction:
       se = evapObj%energy (iz,ia-1) + evapObj%energy (0,1) - evapObj%energy (iz,ia)
!  KKG  11/24/04
       if (iz > 0 .and. iz <= 98 .and. in > 0 .and. in <= 150) then
          se = se + paire0(iz, in)
       endif
       !   Calculate fission barrier:
       temp = compound%numBaryons
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message, 1000) '448'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       x = compound%numProtons*compound%numProtons/temp
       ef = x*(a1*x - a2) + a3 + se
!   Excited energy is below fission barrier:
       if (ef > u) return    
!   Calculate level density parameter for neutron emission:
       an = 0.125d0*(compound%numBaryons - one)
       temp = an
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '122, 123'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       an1 = hlf/temp
       an2 = hlf*hlf*an1/temp
       !   Calculate level density parameter for fission:
       af = x - c3
       af = an*(c1 + c2*af*af)
       af = af*compound%afMultiplier
       a1thrd = ato3rd(ia)
       temp = an*(u - se)
       if (temp < 0.0d0) then
          temp = 0.01d0
          write(evapObj%io%message,1100) '134'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       ss = two*sqrt(temp)
!   Calculate Gamma_n/Gamma_f:
       if (ss > 10.d0) then
!   I0 = J0:
          i0 = an1*(ss - one)
!   I1 = J1: ! Mistake coming from original Atchison; fixed (AJS 10/03)
          i1 = an2*(ss*(ss + ss - 6.d0) + 6.d0)
          gamnf = const*(((0.76d0*i1 - 5.d-2*i0)*a1thrd + 1.93d0*i1) &
               & *a1thrd +  1.66d0*i0)
          temp = af*(u - ef)
          if (temp < 0.0d0) then
             temp = 0.01d0
             write(evapObj%io%message,1100) '148'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          s2 = two*sqrt(temp)
          exps = zro
          if ((ss - s2) > -150.d0) exps = exp(ss - s2)
          temp = s2 - one
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(evapObj%io%message,1000) '156'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          gamnf = gamnf*exps*af/(temp)
       else
          e1 = exp(ss)
          i0 = ((ss - one)*e1 + one)*an1
          i1 = ((6.d0 + ss*(ss + ss - 6.d0))*e1 + ss*ss - 6.d0)*an2
          gamnf = const*(((0.76d0*i1 - 5.d-2*i0)*a1thrd + 1.93d0*i1) &
               & *a1thrd + 1.66d0*i0)
          temp = af*(u - ef)
          if (temp < 0.0d0) then
             temp = 0.01d0
             write(evapObj%io%message,1100) '168'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          s2 = two*sqrt(temp)
          temp = af
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(evapObj%io%message,1000) '174'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          e2 = ((s2 - one)*exp(s2) + one)/temp
          temp = e2
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(evapObj%io%message,1000) '180'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          gamnf = gamnf/temp
       endif
       call evapObj%drein1 (1, ss, an, eye1, eye0)
       call evapObj%drein2 (ss, an, eye2)
       eye10 = eye1 + beta*eye0
       if (eye10 < div0Lim)  then
          write (evapObj%io%message, 2000) compound%numBaryons, compound%numProtons, &
               & u, eye1, eye0, eye2, ss, se, &
               & calcVars%smalla(1), beta
          call evapObj%io%print(2, 2, evapObj%io%message)
          return
       endif
       epsav = (eye2 + beta*eye1)/eye10
       temp = epsav + 7.d0
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '200'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       agoes = (u - 7.d0)/(temp)
       agoes = max(agoes, one)
    endif
    temp = (one + gamnf)*agoes
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) '209'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    fprob = one/(temp)
    return
! ======================================================================
1000 format("Divide by zero error prevented in 'fprob.f90', line ", A)
1100 format("Square root error prevented in 'fprob.f90', line ", A)
2000 format("In 'fprob' routine: a=", f5.2, ", z=", f5.2, &
          & ", u=", f5.2, ", eye1=", f5.2, ", eye0=", f5.2, ", eye2=", &
          & f5.2, ", ss=", f5.2, ", se=", f5.2, ", am=", f5.2, ", beta=", &
          & f5.2)
! ======================================================================
  end function fprob
