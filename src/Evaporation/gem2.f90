
  function efms (evapObj, z, a)

! ======================================================================
!
!  EFMS
!  Fission barrier given by Myers & Swiatecki (PRC 60, 014606, 1999)
!  Note that this fit is only good for Z >= 70; Z <70 results may be
!  suspect!!!  AJS  10/22/03
!
! ====================================================================
! <variables>
!     a   :   the mass of a fissioning nucleus      (IN)
!     z   :   the charge of a fissioning nucleus    (IN)
!   efms  :   fission barrier  [MeV]                (OUT)
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, two, ato3rd
    use evaporationFissionData, only: shellc

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64),       intent(in   ) :: a
    real(real64),       intent(in   ) :: z
    real(real64)                      :: efms

    integer(int32) :: ia, in, iz
    real(real64)   :: ai, c, sh, ss, temp, x, xx

! ======================================================================

!    8/15/1999
    real(real64), parameter :: x0 = 48.5428d0
    real(real64), parameter :: x1 = 34.15d0

! ======================================================================

    efms = zro
    iz = nint(z)
    ia = nint(a)

    c = 1.9d0 + (z - 80.d0)/75.d0
    temp = a
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) '783'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    ai = one - two*(z/temp)
    xx = one - c*ai**2
    ss = xx*ato3rd(ia)**2
    temp = a*xx
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) '791'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    x = z**2/(temp)
    in = ia - iz
    if (in <= 0 .or. in > 250 .or. iz > 150 .or. iz < 1) then
       sh = zro
    else
       sh = shellc (iz, in)
    endif

    if (x >= x1 .and. x <= x0) then
       efms = ss*f1(x0, x) - sh
    elseif (x >= 20.d0 .and. x < x0) then
       efms = ss*f2(x1, x) - sh
    else
       efms = -one
    endif

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'gem2.f90', line ", A)
! ======================================================================
  end function efms


  function f1 ( x0, t )

! ======================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    real(real64), intent(in   ) :: x0
    real(real64), intent(in   ) :: t
    real(real64)                :: f1

! ======================================================================

    f1 = 1.99749d-4*(x0 - t)**3

    return
! ======================================================================
  end function f1

  function f2 ( x1, t )

! ======================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    real(real64), intent(in   ) :: x1
    real(real64), intent(in   ) :: t
    real(real64)                :: f2

! ======================================================================

    f2 = 5.95553d-1 - 0.124136d0*(t - x1)

    return
! ======================================================================
  end function f2



! ======================================================================

  function radgem(ia)

! ======================================================================
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
!  =====================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, ato3rd

    implicit none
    integer(int32), intent(in   ) :: ia
    real(real64)                  :: radgem

! ======================================================================

    real(real64), parameter :: rho0 = 1.414d0

! ======================================================================

    if (ia < 10) then
       if (ia == 9) radgem = 3.25d0
       if (ia == 8) radgem = 2.83d0
       if (ia == 7) radgem = 2.42d0
       if (ia == 6) radgem = 2.02d0
       if (ia == 5) radgem = 2.02d0
       if (ia <= 4) radgem = 1.2d0
       if (ia == 1) radgem = zro
    else
       radgem = rho0*ato3rd(ia) + one
    endif
    return

! ======================================================================
  end function radgem

  function dost (i, z)

! ======================================================================
!
! <variables>
!   t   : set in subroutine setup
!       This routine was originally in the HETC code.
!
! ======================================================================
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: one, ten
    use evaporationFissionData, only: t

    implicit none
    integer(int32), intent(in   ) :: i
    real(real64),   intent(in   ) :: z
    real(real64)                  :: dost

    integer(int32) :: n
    real(real64)   :: x

! ======================================================================

!      i = 1 -> Calculate kp
!      i = 2 -> Calculate k_alpha
!      i = 3 -> Calculate cp

! Linear interpolation, flat extrapolation, LMK
    if (z >= 50) then
       dost = t(i,4)
    elseif (z <= 20) then
       dost = t(i,1)
    else
       n = int(z/ten)
       x = ten*(dble(n) + one)
       x = (x - z)/ten
       dost = x*t(i,n-1) + (one - x)*t(i,n)
    endif
    return

! ======================================================================
  end function dost

! ======================================================================

  function rb (evapObj, a1, ia2, z1, iz2, j)

! ======================================================================
!
!  RB
!    Calculate nuclear radius for a geometric cross section
! ====================================================================
!
! <variables>
!     a1  :   mass of nucleus #1                    (IN)
!     z1  :   charge  of nucleus #1                 (IN)
!    ia2  :   mass of nucleus #2                    (IN)
!    iz2  :   charge  of nucleus #2                 (IN)
!     j   :   type of nucleus #2 (index of ejectile)(IN)
!    ck   :   transmission probability              (IN)
!    rb   :   Nuclear radius [fm]                   (OUT)
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, ato3rd
    use evaporationFissionData, only: aMax

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64),       intent(in   ) :: a1
    integer(int32),     intent(in   ) :: ia2
    real(real64),       intent(in   ) :: z1
    integer(int32),     intent(in   ) :: iz2
    integer(int32),     intent(in   ) :: j
    real(real64)                      :: rb

    integer(int32) :: ia1, l
    real(real64)   :: a1thrd, a2thrd, r1, r2, rr, rzero, temp, ziz

! ======================================================================

    real(real64), parameter                  :: r0     = 1.5d0
    real(real64), parameter, dimension(aMax) :: rmat = &
         & [ ( 1.12d0*ato3rd(l) - 0.86d0/ato3rd(l), l=1, aMax ) ]

! ======================================================================

    ia1 = nint(a1)
    a1thrd = ato3rd(ia1)
    a2thrd = ato3rd(ia2)

    if (evapObj%options%inverseParameter > zro) then
!   Simple form
       if (evapObj%options%inverseParameter == 10.d0) then
          rr = r0
       else
          rr = evapObj%options%inverseParameter
       endif
       r1 = rr*a1thrd
       r2 = rr*a2thrd
       rb = r1 + r2
    elseif (evapObj%options%inverseParameter == -one .and. j > 6) then
!   Expression in NP A475 (1987) 663 (for heavy ions)
       ziz = z1*dble(iz2)
       temp = one + 9.443d-3*ziz
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '1161' 
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       rzero = 2.173d0*(one + 6.103d-3*ziz)/(temp)
       r1 = a1thrd
       r2 = a2thrd
       rb = rzero*(r1 + r2)
    elseif (evapObj%options%inverseParameter < zro .and. j <= 6) then
!   FB form  ...2001/3/15
       r1 = radgem (ia1)
       r2 = radgem (ia2)
       rb =  r1 + r2
    else
!   Dostrovsky et. al.
       if (j <= 6) then
          if (j > 2) then
             rb  =  a1thrd + a2thrd
          else
             rb  =  a1thrd
          endif
          rb  =  rb*r0
       else
!    Matsuse et al. PRC 26 (1982) 2338
          rb = rmat(ia1) + rmat(ia2) + 2.85d0
       endif
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'gem2.f90', line ", A)
! ======================================================================
  end function rb

  function rho (evapObj, ia, iz, u, delta, calcVars)

! ======================================================================
!
!    Calculates unnormalized level density as a function of energy.
!
!     variable       IN/OUT
!     ia             IN        mass before emission (=A)
!     iz             IN        charge before emission (=Z)
!     u              IN        excitation energy before emission
!     delta          IN        pairing energy
!     smalla         IN        level density parameter at U
!
! -------------------------------------------------------------------
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Modified by AJS, January, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: two
    use evaporationFissionData, only: ux0

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    integer(int32),     intent(in   ) :: ia
    integer(int32),     intent(in   ) :: iz
    real(real64),       intent(in   ) :: u
    real(real64),       intent(in   ) :: delta
    type(evaporationCalculation), intent(in   ) :: calcVars
    real(real64)                      :: rho

    integer(int32) :: nn
    real(real64)   :: ex0, s0, smalla14, temp, umd, umd14, ut

! ======================================================================

    ex0 = ux0(ia) + delta
    s0 = two*sqrt(abs(calcVars%smalla0*u))
    if (u >= ex0) then
!  Eq. (44) of CEM03 Manual: [modified by mysterious exp(-s0) factor]
       smalla14 = sqrt(sqrt(abs(calcVars%smalla0)))
       umd = u - delta
       umd14 = sqrt(sqrt(abs(umd)))
       temp = smalla14*umd*umd14
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '1287'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       rho = exp(two*sqrt(abs(calcVars%smalla0*umd)) - s0) &
            & /(temp)
    else
!  Eq. (46) of CEM03 Manual: [modified by mysterious exp(-s0) factor]
       nn = ia - iz
       ut = u*evapObj%data%taux0(iz,nn)
       rho = evapObj%data%fact10(iz,nn) * evapObj%data%taux0(iz,nn) * &
            & exp(ut + evapObj%data%sx0(iz,nn) - s0)
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'gem2.f90', line ", A)
! ======================================================================
  end function rho

  function pe (evapObj, e, q, delta, bet, tau, e0, ex, smalla, x, calcVars)

! ======================================================================
!
!   The unnormalized emission probability for evaporation of a particle
!   as a function of energy
!
!   variable :   IN/OUT
!     e      :     IN      excitation energy of nucleus before emission
!     q      :     IN      Q-value
!     delta  :     IN      pairing energy
!     bet    :     IN      beta in the Dostrovsky inverse cross section
!                          approximation; beta for n, -V for chgd. par.
!     tau    :     IN      auxiliary quantity for low excitation energy
!                          defined in EYE subroutine 
!     e0     :     IN      auxiliary quantity; defined after Eq. (46)
!     ex     :     IN      auxiliary quantity; point of joining of the
!                          low e and high e approximations to l. d. parm.
!                          = Ux0 + delta; Ux0 defined in EYE
!    smalla  :     IN      level-density parameter
!     x      :     IN      independent variable; KE of ejectile
!     pe     :    OUT      unnormalized emission probability as a funct
!                          of outgoing KE; see Eq. (39) of CEM03 man'l.
!
!    Called by: CALR, SELECTE
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: real64
    use evaporationParams, only: zro, two

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64),       intent(in   ) :: e
    real(real64),       intent(in   ) :: q
    real(real64),       intent(in   ) :: delta
    real(real64),       intent(in   ) :: bet
    real(real64),       intent(in   ) :: tau
    real(real64),       intent(in   ) :: e0
    real(real64),       intent(in   ) :: ex
    real(real64),       intent(in   ) :: smalla
    real(real64),       intent(in   ) :: x
    type(evaporationCalculation), intent(in   ) :: calcVars
    real(real64)                      :: pe

    real(real64)   :: eqex, eqx, s0, sqax, temp, xxx

! ======================================================================

    if (e - q <= zro) then
       pe = zro
       return
    endif
    s0 = two*sqrt(abs(calcVars%smalla0*e))
    eqx = e - q - x
    eqex = e - q - ex
    if (x > eqex .and. eqx >= zro) then
       temp = tau
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '1542'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       pe = exp((eqx - e0)/temp - s0)/temp
!  Modified 2001/4/13:
    elseif (x >= zro .and. x <= eqex) then  
       xxx = eqx - delta
       sqax = sqrt(abs(smalla*xxx))
       pe = exp(two*sqax - s0)
       temp = sqax
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '1553'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       pe = pe/sqrt(abs(temp))
       temp = xxx
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '1559'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       pe = pe/temp
    else
       pe = zro
    endif
    pe = pe*(x + bet)
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'gem2.f90', line ", A)
! ======================================================================
  end function pe

  subroutine calr (evapObj, aa, zz, u, q, v, delta, r, calcVars)

! ======================================================================
!
!    Called by: GAMMA
!
!    Calls: GETA, PE
!
!    Decay width
!    Related to the probability for emission under the barrier of an
!    evaporated complex particle (AJS)
!
!     variable       IN/OUT
!     aa             IN        residual mass after emission (=A-Aj)
!     zz             IN        charge number of res nuclei after emission
!     u              IN        excitation energy of nuclei before emission
!     q              IN        Q-value
!     v              IN        kv in the equation
!     delta          IN        pairing energy
!     r              OUT       r in the equation
! -------------------------------------------------------------------
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one ,two

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64),       intent(in   ) :: aa
    real(real64),       intent(in   ) :: zz
    real(real64),       intent(in   ) :: u
    real(real64),       intent(in   ) :: q
    real(real64),       intent(in   ) :: v
    real(real64),       intent(in   ) :: delta
    real(real64),       intent(  out) :: r
    type(evaporationCalculation), intent(inout) :: calcVars

    integer(int32) :: ifuri, imax, isdum, izz, nn
    real(real64)   :: aaa, alpha, ax, dt, e0, emin, ex, sx, tau, &
         & temp, tempOld, temp2, ux, x, xold, xxx

! ======================================================================

    izz = nint (zz)
    nn = nint (aa - zz)
    temp = zz
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) '1625'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    alpha = 0.869d0 + 9.91d0/temp
    temp = aa
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) '1631'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    ux = 2.5d0 + 150.d0/temp
    ex = ux + delta
    ax = geta (evapObj%data%alev(), ex, izz, nn, isdum)
    temp = ux
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) '1639'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    tau = sqrt(abs(ax/temp)) - 1.5d0/temp
    temp = tau
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) '1645'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    tau = one/temp
    sx = two*sqrt(abs(ax*ux))
    e0 = ex - tau*(log(tau) - 0.25d0*log(ax) - 1.25d0*log(ux) + sx)
    
!   ..Calculate exact decay width
    imax = 1000
    emin = v*0.1d0
    xold = zro
    tempold = zro
    r = zro
    do ifuri = 0,imax
       xxx = u - q - v
       temp = dble(imax)
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '1662'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       x = emin + xxx*dble(ifuri)/temp
       aaa = geta (evapObj%data%alev(), xxx, izz, nn, isdum)
       temp = evapObj%pe (u, q, delta, -v, tau, e0, ex, aaa, x, calcVars)
       if (x < v+one) then
          temp2 = x - v
          if (temp2 < div0Lim .and. temp2 > -div0Lim) then
             temp2 = div0Lim
             write(evapObj%io%message,1000) '1671'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          temp = temp/(temp2)
          temp2 = one + v
          if (temp2 < div0Lim .and. temp2 > -div0Lim) then
             temp2 = div0Lim
             write(evapObj%io%message,1000) '1677'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          temp = temp/(temp)*exp(alpha*(x - v - one))*x
       endif
       dt = (tempold + temp)*(x - xold)/two  
       r = r + dt
       temp2 = r
       if (temp2 < div0Lim .and. temp2 > -div0Lim) then
          temp2 = div0Lim
          write(evapObj%io%message,1000) '1686'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       if (x > v .and. (dt/temp2) < 1.d-4) return
       xold = x
       tempold = temp
    end do
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'gem2.f90', line ", A)
! ======================================================================
  end subroutine calr

  function ckcal (j, z, a)

! ====================================================================
!
! <variables>
!     j   :   ejectile  ID                      (IN)
!     a   :   mass of daughter nucleus          (IN)
!     z   :   charge  of daughter nucleus       (IN)
!
! ====================================================================
! 
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one

    implicit none
    integer(int32), intent(in   ) :: j
    real(real64),   intent(in   ) :: z
    real(real64),   intent(in   ) :: a
    real(real64)                  :: ckcal

    integer(int32) :: i

! ======================================================================

    integer(int32), parameter :: jul = 65
    real(real64), dimension(jul), parameter :: c1 = &
         & [     0.0615d0, 0.0556d0, 0.0530d0, 0.0484d0, 0.0468d0, &
         &       0.0440d0, 0.0423d0, 0.0419d0, 0.0409d0, 0.0404d0, 0.0395d0, &
         &       0.0397d0, 0.0382d0, 0.0378d0, 0.0374d0, 0.0370d0, 0.0380d0, &
         &       0.0368d0, 0.0362d0, 0.0358d0, 0.0354d0, 0.0358d0, 0.0352d0, &
         &       0.0349d0, 0.0351d0, 0.0353d0, 0.0354d0, 0.0353d0, 0.0352d0, &
         &       0.0353d0, 0.0352d0, 0.0349d0, 0.0348d0, 0.0346d0, 0.0347d0, &
         &       0.0345d0, 0.0345d0, 0.0344d0, 0.0344d0, 0.0343d0, 0.0345d0, &
         &       0.0342d0, 0.0344d0, 0.0346d0, 0.0346d0, 0.0347d0, 0.0346d0, &
         &       0.0345d0, 0.0345d0, 0.0346d0, 0.0347d0, 0.0347d0, 0.0347d0, &
         &       0.0346d0, 0.0347d0, 0.0347d0, 0.0347d0, 0.0346d0, 0.0345d0, &
         &       0.0346d0, 0.0347d0, 0.0347d0, 0.0347d0, 0.0347d0, 0.0347d0  ]

    real(real64), dimension(jul), parameter :: c2 = &
         & [     0.0167d0, 0.0135d0, 0.0134d0, 0.0122d0, 0.0122d0, &
         &       0.0116d0, 0.0111d0, 0.0110d0, 0.0107d0, 0.0105d0, 0.0106d0, &
         &       0.0105d0, 0.0101d0, 0.0097d0, 0.0100d0, 0.0104d0, 0.0098d0, &
         &       0.0099d0, 0.0101d0, 0.0103d0, 0.0104d0, 0.0098d0, 0.0099d0, &
         &       0.0099d0, 0.0096d0, 0.0094d0, 0.0096d0, 0.0101d0, 0.0090d0, &
         &       0.0092d0, 0.0097d0, 0.0103d0, 0.0109d0, 0.0114d0, 0.0104d0, &
         &       0.0109d0, 0.0114d0, 0.0119d0, 0.0123d0, 0.0128d0, 0.0130d0, &
         &       0.0122d0, 0.0125d0, 0.0128d0, 0.0131d0, 0.0134d0, 0.0147d0, &
         &       0.0151d0, 0.0155d0, 0.0159d0, 0.0162d0, 0.0166d0, 0.0169d0, &
         &       0.0140d0, 0.0143d0, 0.0148d0, 0.0152d0, 0.0156d0, 0.0126d0, &
         &       0.0129d0, 0.0132d0, 0.0137d0, 0.0141d0, 0.0145d0, 0.0149d0 ]

    real(real64), dimension(jul), parameter :: c3 = &
         & [     0.3227d0, 0.4067d0, 0.4374d0, 0.4938d0, 0.5120d0, &
         &       0.5460d0, 0.5691d0, 0.5735d0, 0.5856d0, 0.5953d0, 0.6035d0, &
         &       0.6029d0, 0.6216d0, 0.6310d0, 0.6341d0, 0.6370d0, 0.6262d0, &
         &       0.6421d0, 0.6460d0, 0.6499d0, 0.6537d0, 0.6527d0, 0.6573d0, &
         &       0.6616d0, 0.6648d0, 0.6668d0, 0.6679d0, 0.6680d0, 0.6696d0, &
         &       0.6709d0, 0.6714d0, 0.6717d0, 0.6718d0, 0.6718d0, 0.6747d0, &
         &       0.6750d0, 0.6751d0, 0.6751d0, 0.6750d0, 0.6748d0, 0.6746d0, &
         &       0.6781d0, 0.6781d0, 0.6779d0, 0.6778d0, 0.6778d0, 0.6844d0, &
         &       0.6839d0, 0.6833d0, 0.6829d0, 0.6823d0, 0.6818d0, 0.6812d0, &
         &       0.6827d0, 0.6823d0, 0.6818d0, 0.6813d0, 0.6808d0, 0.6810d0, &
         &       0.6809d0, 0.6808d0, 0.6804d0, 0.6801d0, 0.6796d0, 0.6792d0 ]

! ======================================================================

    if (j <= 1 .or. z <= zro .or. a < one) then
       ckcal = zro 
    else
       i = min(j-1,jul)
       ckcal = c1(i)*log(z) + c2(i)*log(a) + c3(i)
    endif

    return

! ======================================================================
  end function ckcal

