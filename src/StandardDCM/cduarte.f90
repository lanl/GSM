
  function cduarte (sDCM, tin1, ie1, ie2, results)

! ======================================================================
!
!     Angular distribution simulation for n + p and p + p
!     for energies < 2GeV using Duarte's approximations:
!     n + p : ds/dom ~ exp(b*t) + a*exp(b*u) + c*exp(alc*u)
!     n + n OR p + p:  ds/dom ~ exp(b*t) + a*exp(b1*t)
!     with Mandelstam's variables t, u
!
!   Called by: ELEX
!     H.Duarte http://www.fjfi.cvut.cz/con_adtt99/papers/Mo-o-c17.pdf
!
!   INPUT:  ie1 = charge of particle 1;
!           ie2 = charge of particle 2;
!          tin1 = CM kinetic energy in units of GeV.
!
!     Written by K. K. Gudima at LANL, Fall 2003.
!     Rewritten by A. J. Sierk, LANL T-16, March, 2004.
!     Edited by AJS, LANL T-2, December, 2011.
!     Edited by LMK, XCP-3, July 2013 (inclued error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, hlf, one, two, four, thsn, pi, twpi, &
         & emnucg

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    real(real64),   intent(in   ) :: tin1
    integer(int32), intent(in   ) :: ie1
    integer(int32), intent(in   ) :: ie2
    type(StandardDCMResults), intent(inout) :: results
    real(real64)                  :: cduarte

    real(real64) :: a, alc, b, b1, c, c1, c2, cn, cn1, cn2, cn3, cts, &
         & el3h, el67, el6e, elab, ps2, ps2b, ps2c, r1, r2, t, temp, &
         & tm, twpia, u, um

! ======================================================================

    real(real64), parameter :: eps = 1.0e-04

! ======================================================================

    cduarte = zro
    elab = tin1*thsn
    ps2 = hlf*tin1*emnucg
    el3h = elab - 300.d0
    el67 = elab - 670.d0
!   p + p or n + n:
    if (ie1 == ie2) then
       a  =  0.2d0
       b1 =  zro
       if (el3h <= zro) then
          b  =  zro
       elseif (elab <= 670.d0) then
          b  =  9.87d-8*el3h**3
       elseif (elab <= 1100.d0) then
          b  =  4.56d-3*el67 + 4.76d0
          a  =  97.02d3*exp(-2.0d-2*elab) + 0.053d0
          b1 =  9.72d-8*exp(-5.0d-3*el67)*el67**3
       else
          temp = el3h**2.23d0
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "63"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          b  =  7.4d0/(one + 3.0d5/temp)
          a  =  0.28d0*exp(-1.5d-3*elab)
          b1 =  1.94d0*exp(-7.0d-4*elab)
       endif
       if (b <= eps) then
          cn1  =  two*twpi
       else
          temp = ps2*b
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "75"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          cn1  =  twpi*(one - exp(-two*ps2*b))/(temp)
       endif
       twpia = twpi*a
       if (b1 <= eps) then
          cn2  =  two*twpia
       else
          temp = ps2*b1
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "86"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          cn2 =   twpia*(one - exp(-two*ps2*b1))/(temp)
       endif
       cn  =  cn1 + cn2
       temp = cn
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message,1000) "94 and 95"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       c1  =  cn1/temp
       c2  =  cn2/temp
       r1 = sDCM%rang()
       r2 = sDCM%rang()
       tm = -two*ps2
       if (r1 <= c1) then
          if (b <= eps) then
             t =  tm*r2
          else
             temp = b
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(sDCM%io%message,1000) "108"
                call sDCM%io%print(4, 3, sDCM%io%message)
             end if
             t = log(one - r2*(one - exp(tm*b)))/temp
          endif
       else
          if (b1 <= eps) then
             t =  tm*r2
          else
             temp = b1
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(sDCM%io%message,1000) "119"
                call sDCM%io%print(4, 3, sDCM%io%message)
             end if
             t = log(one - r2*(one - exp(tm*b1)))/temp
          endif
       endif
       temp = tm
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message,1000) "127"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       cts = one - t/temp
       if (sDCM%rang() <= hlf) then
          cduarte  =  cts
       else
          cduarte  = -cts
       endif
!  n + p  OR  p + n:
    elseif ((ie1 == 0 .and. ie2 == 1) .or. (ie1 == 1 .and. ie2 == 0)) then
       el6e = exp(-6.d-3*elab)
       b  =  25.0d0*el6e + 2.0d-3*elab + 2.8d0
       a  =  (3.0d-3*elab + 0.1d0)*exp(1.55d0 - 4.9d-3*elab) + eps*elab
       c  =  6.0d-5*elab*elab*el6e + 15.0d-5*elab
       if (elab < 100.d0) then
          alc = 330.0d0 - elab
       else
          alc = 80.0d0 + 15.0d3/elab
       endif
       ps2b = ps2*b
       temp = ps2b
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message,1000) "150"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       cn1  =  pi*(one - exp(-four*ps2b))/temp
       cn2  =  cn1*a
       ps2c = ps2*alc
       temp = ps2c
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message,1000) "158"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       cn3  =  pi*c*(one - exp(-four*ps2c))/temp
       cn  =  cn1 + cn2 + cn3
       temp = cn
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message,1000) "165 and 166"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       c1 = cn1/temp
       c2 = cn2/temp
       tm = -four*ps2
       um = tm
       r1 = sDCM%rang()
       r2 = sDCM%rang()
       temp = b
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message,1000) "180 and 192"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       if (r1 <= c1) then
          if (b <= eps) then
             t =  tm*r2
          else
             t = log(one - r2*(one - exp(tm*b)))/temp
          endif
          temp = tm
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "187"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          cts  =  one - two*t/temp
       elseif (r1 <= (c1 + c2)) then
          if (b <= eps) then
             u =  um*r2
          else
             u = log(one - r2*(one - exp(um*b)))/temp
          endif
          temp = um
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "199"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          cts  = -one + two*u/temp
       else
          if (alc <= eps) then
             u =  um*r2
          else
             temp = alc
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(sDCM%io%message,1000) "209"
                call sDCM%io%print(4, 3, sDCM%io%message)
             end if
             u = log(one - r2*(one - exp(um*alc)))/temp
          endif
          temp = um
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(sDCM%io%message,1000) "216"
             call sDCM%io%print(4, 3, sDCM%io%message)
          end if
          cts = -one + two*u/temp
       endif
       if (ie1 == 1 .and. ie2 == 0) then
!  parametrization is for n + p case; change sign for p + n
          cduarte = -cts
       else
          cduarte =  cts
       endif
    else
       write(sDCM%io%message, 2000) ie1, ie2
       call sDCM%io%print(2, 3, sDCM%io%message)
       results%simState = 2   ! Flag that an error occurred
       return
    endif


    if (abs(cduarte) > one) then
       cduarte = sign(one, cduarte)
    elseif (abs(cduarte) < 1.0d-10) then
       cduarte = zro
    endif

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'cduarte.f90' line(s) ", A)
2000 format("ie1 (", i3, ") and ie2 (", i3, ") are wrong in ", &
          & "'cduarte' - unknown error.")
! ======================================================================
  end function cduarte
