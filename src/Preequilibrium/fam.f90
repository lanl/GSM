
  function fam (preeqObj, a, z, e, calcVars, j)

! ======================================================================
!
!   Calculates a level density parameter, a
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk  LANL  T-2  February, 1996.
!   Modified to add levelDenParam = 11, 12; May, 1996.
!   Modified by AJS, March, 1999.
!
!   "Last" change: 12-AUG-2003
!   Edited by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use preequilibriumParams, only: one

    implicit none
    class(Preequilibrium), intent(inout) :: preeqObj
    real(real64),          intent(in   ) :: a
    real(real64),          intent(in   ) :: z
    real(real64),          intent(in   ) :: e
    type(preequilibriumCalculation), intent(inout) :: calcVars
    integer(int32),        intent(in   ) :: j
    real(real64)                         :: fam

    integer(int32) :: ln
    real(real64)   :: afg, al, be, bs, de, fct, ga, sh, temp, x

! ======================================================================
!   Set of parameters derived by Sierk in May 1996, following method
!   of Iljinov, et al., using nearly the same data set; several
!   errors in J's, Bn values, and a few apparent misprints in the
!   Tables 1 & 2 were corrected. Using Moller-Nix-Myers-Swiatecki (1992)
!   shell corrections; pairing = 11/sqrt(A); without collective effects;
!   without A dependnce of ga. Gives nearly identical fit to levelDenParam = 7,
!   above, but without the time-consuming calculation of shell 
!   corrections following Myers & Swiatecki.

    if (preeqObj%options%levelDenParam == 11) then
       al =  0.1450d0
       be = -0.0694d0
       ga =  0.0511d0

!   Set of parameters derived by Sierk in May 1996, following method
!   of Iljinov, et al., using nearly the same data set; several
!   errors in J's, Bn values, and a few apparent misprints in the
!   Tables 1 & 2 were corrected. Using Moller-Nix-Myers-Swiatecki (1992)
!   shell corrections; Moller-Nix pairing gaps for the shifted Fermi
!   gas level density formula; without collective effects;
!   without A dependnce of ga. Gives not quite as good fit as 
!   levelDenParam=11, above. [2nd version uses rho = 0.010 x exp(2 x sqrt(aE*))]

    elseif (preeqObj%options%levelDenParam == 12) then
       al =  0.1463d0
       be = -0.0716d0
       ga =  0.0542d0
    endif

    if (e == -one) then
!   If finding saddle point level density, do not include GS
!   shell & pairing correction (e input as -1.0).
!   For saddle-point level density, find bs:
       if (z < 15) then
          bs = 1.25d0
       else
          ln = 0 
          call preeqObj%fissBarr%bsfit (z, a, ln, bs)
       endif
       fct = one
    else
       bs = one
       sh = preeqObj%molEnergy%shellEnergy(a, z)
       x = ga*e
       if (x <= 0.00001d0) then
          fct = ga
       else
          temp = e
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(preeqObj%io%message,1000) "91"
             call preeqObj%io%print(4, 3, preeqObj%io%message)
          end if
          fct = (one - exp(-x))/temp
       endif
       fct = one + fct*sh
    endif
    if (j == 0) then
       de = calcVars%athrd
    else
       de = calcVars%afjthr(j)
    endif
    temp = de
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(preeqObj%io%message,1000) "105"
       call preeqObj%io%print(4, 3, preeqObj%io%message)
    end if


    afg = al + be*bs/temp
    fam = afg*fct

    return
! ======================================================================
1000 format("Divide by zero error prevented in 'fam.f90', line ", A)
! ======================================================================
  end function fam
