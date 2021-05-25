
  subroutine eye (evapObj, iaa, izz, in, u, q, v, delta, smalla, &
       & beta, r, calcVars)

! ======================================================================
!
!    Calculates the integral in the decay width calculation; the 
!    quantity in the square brackets of Eq. (48) from the CEM03 manual.
!
!    Called by: GAMMA
!
!    Calls: EY2 EY3
!
!     variable       IN/OUT
!     iaa            In      residual mass after emission ( = A-Aj)
!     izz            In      charge number of res nucleus after emission
!     in             In      neutron number of res nucleus after emiss.   
!     u              In      excitation energy of nucleus before emiss.   
!     q              In      Q-value
!     v              In      v in the equation
!     delta          In      pairing energy
!     smalla         In      level density parameter at U-Q-delta-V
!     beta           In      beta in the equation [= -v for j > 1]
!     r              OUT     r in the equation
!     r              OUT     (decay width with factors removed)
!  [r is the quantity in brackets in Eq. (48) in the CEM03 manual.] AJS
!                                      _______________
!     s              OUT     s  =  2 \/a(u-q-delta-v) (NOT USED!)
!
!     ux0            Ux      in CEM03 Manual
!     ux140         (Ux)^(1/4)
!
! -------------------------------------------------------------------
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Modified by A. J. Sierk, LANL T-16, November, 2004.
!    Modified by AJS, January, 2005.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (Evap class creation)
!
! ====================================================================== 

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams,             only: zro, one, two
    use evaporationFissionData,        only: paire0, ux0

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    integer(int32),     intent(in   ) :: iaa
    integer(int32),     intent(in   ) :: izz
    integer(int32),     intent(in   ) :: in
    real(real64),       intent(in   ) :: u
    real(real64),       intent(in   ) :: q
    real(real64),       intent(in   ) :: v
    real(real64),       intent(in   ) :: delta
    real(real64),       intent(in   ) :: smalla
    real(real64),       intent(in   ) :: beta
    real(real64),       intent(out)   :: r
    type(evaporationCalculation), intent(inout) :: calcVars

    real(real64)   :: ax, bv, ex, exets, expssx, exss0, ext, extx, eye0, &
         & eye1, eye2, eye3, s, s0, sx, t, tau, temp, tp1, tx, uqv, &
         & uqvd, ux

! ======================================================================

    uqv = u - q - v
    if (uqv <= zro) then
       r = zro
       return
    endif

    ux = ux0(iaa)
    ex = ux + paire0(izz,in)
    ax = evapObj%data%ax0(izz,in)
    if (ax <= zro .or. ux <= zro)  then
       r = zro
       return
    endif

    tau = evapObj%data%taux0(izz,in)
!  This tau is exactly T, see above:
    temp = tau
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) '166'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    tau = one/temp
!   AJS 11/02/04:
    if (tau <= zro .or. u <= zro)  then
       r = zro
    else     
       tx  = evapObj%data%tx0(izz,in)
       t   = uqv*evapObj%data%taux0(izz,in)
       sx  = evapObj%data%sx0(izz,in)
       temp = calcVars%smalla0*u
       if (temp < 0.0d0) then
          temp = 0.01d0
          write(evapObj%io%message,1100) '183'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       s0 = two*sqrt(temp)
       bv = beta + v
       exets = evapObj%data%fact10(izz,in)*exp(sx - s0)
       tp1 = t + one
       if (uqv <= ex) then
          ext = exp(t)
          eye1 = tau*(ext - tp1)
          eye0 = ext - one
!   Eq(47) in CEM03 manual:
          r = exets*(eye1 + bv*eye0)
       else
          uqvd = uqv - delta
          if (uqvd <= zro) then
             s = zro
          else
             temp = smalla*uqvd
             if (temp < 0.0d0) then
                temp = 0.01d0
                write(evapObj%io%message,1100) '203'
                call evapObj%io%print(4, 3, evapObj%io%message)
             end if
             s = two*sqrt(temp)
          endif
          extx = evapObj%data%extx0(izz, in)
          exss0 = exp(s - s0)
          eye1 = exets*tau*(extx*(tp1 - tx) - tp1)
          eye0 = exets*(extx - one)
          expssx = exp(sx - s)
          eye3 = exss0*evapObj%ey3 (s, sx, smalla, expssx)
          eye2 = exss0*evapObj%ey2 (s, sx, expssx)
!   Eq(48) in CEM03 manual:
          r =  eye3 + eye1 + bv*(eye0 + eye2)
       endif
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'eye.f90', line ", A)
1100 format("Square root error prevented in 'eye.f90', line ", A)
! ======================================================================
  end subroutine eye

  function ey2 (evapObj, s, sx, expssx)

! ======================================================================
!
!   Called by: EYE
!
!   Integral I_2 defined on pg. 18 of the CEM03 manual:
!
! -------------------------------------------------------------------
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Modified by AJS, January, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (Evap class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: real64
    use evaporationParams,             only: zro, one, cst0

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64),       intent(in   ) :: s
    real(real64),       intent(in   ) :: sx
    real(real64),       intent(in   ) :: expssx
    real(real64)                     :: ey2

    real(real64) :: s2, s3, sqs, sqx, sx2, sx3, temps, tempsx, temp, &
         & tempx

! ======================================================================

    temps = s
    if (temps <= zro) then
       temp = zro
    else
       if (temps < div0Lim .and. temps > -div0Lim) then
          temps = div0Lim
          write(evapObj%io%message,1000) '268'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       s2 = temps * temps
       s3 = s2 * temps
       sqs = sqrt(temps)
       temp = (one + 1.5d0/temps + 3.75d0/s2 + 13.125d0/s3)/(temps*sqs)
    endif

    if (sx <= zro) then
       tempx = zro
    else
       tempsx = sx
       if (tempsx < div0Lim .and. tempsx > -div0Lim) then
          tempsx = div0Lim
          write(evapObj%io%message,1000) '280'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       sqx = sqrt(tempsx)
       sx2 = tempsx * tempsx
       sx3 = sx2 * tempsx
       tempx = (one + 1.5d0/tempsx + 3.75d0/sx2 + 13.125d0/sx3)/(tempsx*sqx)
    endif

    ey2 = cst0*(temp - expssx*tempx)

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'eye.f90', line ", A)
! ======================================================================
  end function ey2

  function ey3 (evapObj, s, sx, a, expssx)

! ======================================================================
!
!   Called by: EYE
!
!   Integral I_3 defined on pg. 18 of the CEM03 manual:
!
! -------------------------------------------------------------------
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Modified by AJS, January, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (Evap class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: real64
    use evaporationParams, only: two, four, sqrtTwo

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64),       intent(in   ) :: s
    real(real64),       intent(in   ) :: sx
    real(real64),       intent(in   ) :: a
    real(real64),       intent(in   ) :: expssx
    real(real64)                      :: ey3

    real(real64) :: esxs, s2, s3,s4, sqs, sqsx, sx2, sx3, sx4, sx5, &
         & sx6, temps, tempsx, temp, temp2

! ======================================================================

    temps = s
    if (temps < div0Lim) then
       temps = 0.01d0
       write(evapObj%io%message,1200) '334-357'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    tempsx = sx
    if (tempsx < div0Lim) then
       tempsx = 0.01d0
       write(evapObj%io%message,1200) '334-357'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    sqs = sqrt(temps)
    sqsx = sqrt(tempsx)
    esxs = expssx/sqsx
    s2 = temps * temps
    s3 = s2 * temps
    s4 = s2*s2
    sx2 = tempsx * tempsx
    sx3 = sx2 * tempsx
    sx4 = sx2*sx2
    sx5 = sx2*sx3
    sx6 = sx3*sx3
    temp = 325.125d0/(sqs*s4) - &
         & esxs*(324.8d0*s2   +   3.28d0*sx2)/sx6
    temp = temp + 60.d0/(sqs*s3) - &
         & esxs*(59.0625d0*s2 + 0.9375d0*sx2)/sx5
    temp = temp + 13.5d0/(sqs*s2) - &
         & esxs*(12.875d0*s2  +  0.625d0*sx2)/sx4
    temp = temp + four/(sqs*temps) - &
         & esxs*(3.75d0*s2    +   0.25d0*sx2)/sx3
    temp = temp + two/sqs      - &
         & esxs*(1.5d0*s2/sx2 +    0.5d0    )
    temp = temp - esxs*(s2 - sx2)/tempsx
    temp2 = a*sqrtTwo
    if (temp2 < div0Lim .and. temp2 > -div0Lim) then
       temp2 = div0Lim
       write(evapObj%io%message,1000) '365'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    ey3 = temp/(temp2)

    return

! ======================================================================
1000 format("Divide by zero error prevented in ", &
          & "'eye.f90', line ", A)
1200 format("Square root/divide by zero error prevented in ", &
          & "'eye.f90', line ", A)
! ======================================================================
  end function ey3
