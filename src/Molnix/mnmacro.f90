
  function mnmacro (molObj, iz, in)

! ======================================================================
!
!   This function calculates the Moller-Nix macroscopic mass excess
!   for (spherical) nuclei lying outside the tabulated range of Z and N.
!   The inputs are the proton and neutron numbers.
!   Uses formula (62) from At. Data Nucl. data Tables 59, 185 (1995).
!   Uses FRLDM, not FRDM.
!
!   Called by: MOLNIX
!
!    Written by A. J. Sierk  LANL  T-2  November, 1997.
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use molnixParams, only: zro, hlf, one, two, thr, fiv, &
         & thrd, twthrd, pi, ato3rd

    implicit none
    class(Molnix),  intent(inout) :: molObj
    integer(int32), intent(in   ) :: in
    integer(int32), intent(in   ) :: iz
    real(real64)                  :: mnmacro

    real(real64) :: a, ai, akf, athrd, b3, bi, deln, delnp, delp, &
         & fkprp, temp, un, wt, xx, yc, yc2, ys, z
    logical      :: nodd, podd

! ======================================================================

    real(real64), parameter :: &
         &  emp = 7.289034d0,  emn = 8.071431d0, &
         &  esq = 1.4399764d0, ael = 1.433d-5, &
         &   rp = 0.80d0,       r0 = 1.16d0, &
         & aype = 0.68d0,       ay = 0.70d0, &
         & rmac = 4.80d0,        h = 6.6d0, &
         &    w = 30.d0,        av = 16.00126d0, &
         &  akv = 1.92240d0,    as = 21.18466d0, &
         &  aks = 2.345d0,      a0 = 2.615d0, &
         &   ca = 0.10289d0

    real(real64), parameter :: d1 = 3.020833333d0
    real(real64), parameter :: d2 = 0.113541666667d0
    real(real64), parameter :: d3 = 0.00126240079d0

    real(real64), parameter :: c1 = thr*esq/(fiv*r0)
    real(real64), parameter :: c4 = 1.25d0*c1*(1.5d0/pi)**twthrd
    real(real64), parameter :: c5 = (2.25d0*pi)**thrd/r0

! ======================================================================

    un = dble(in)
    z = dble(iz)
    a = z + un
    athrd = ato3rd(in + iz)
    if (a < div0Lim .and. a > -div0Lim) then
       a = div0Lim
       write(molObj%io%message,1000) "292, 323"
       call molObj%io%print(4, 3, molObj%io%message)
    end if
    ai = (un - z)/a
    deln = zro
    delnp = zro
    nodd = .false.
    podd = .false.
    if (2*(in/2).ne.in) then
       temp = ato3rd(in)
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(molObj%io%message,1000) "304"
          call molObj%io%print(4, 3, molObj%io%message)
       end if
       deln = rmac/temp
       nodd = .true.
    endif
    delp = zro
    if (2*(iz/2).ne.iz) then
       temp = ato3rd(iz)
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(molObj%io%message,1000) "314"
          call molObj%io%print(4, 3, molObj%io%message)
       end if
       delp = rmac/temp
       podd = .true.
    endif
    wt = zro
    if (nodd .and. podd) then
       if (athrd < div0Lim .and. athrd > -div0Lim) then
          athrd = div0Lim
          write(molObj%io%message,1000) "322, 325"
          call molObj%io%print(4, 3, molObj%io%message)
       end if
       delnp = h/athrd**2
       if (in == iz) wt = one/a
    endif
    akf = c5*ato3rd(iz)/athrd
    xx = akf*rp
    fkprp = -0.125d0*esq*rp**2*(d1 - d2*xx**2 + d3*xx**4)/r0**3
    ys = r0*athrd/aype
    yc = ys*aype/ay
    if (ys < div0Lim .and. ys > -div0Lim) then
       ys = div0Lim
    end if
    if (yc < div0Lim .and. yc > -div0Lim) then
       yc = div0Lim
       write(molObj%io%message,1000) "344"
       call molObj%io%print(4, 3, molObj%io%message)
    end if
    yc2 = yc*yc

    bi = one - thr/ys**2 + (one + ys)*( two + thr * (one + one/ys)/ys) &
         & * exp(-two*ys)
    b3 = one - fiv*(one - 1.875d0*(one - 1.4d0/yc2)/yc - 0.75d0 * &
         & (one + 4.5d0/yc + 7.0d0*(one + hlf/yc)/yc2)*exp(-two*yc)) / &
         & yc2

    mnmacro = emp*z + emn*un - av*a*(one - akv*ai**2) + as*bi*  &
         (one - aks*ai**2)*athrd**2+ a0 + c1*b3*z**2/athrd -  &
         c4*ato3rd(iz)**4/athrd + fkprp*z**2/a - ca*ai*a +  &
         w*(abs(ai) + wt) + delp + deln - delnp - ael*z**2.39d0

    return
! ======================================================================
1000 format("Divide by zero error prevented in 'mnmacro.f90', line(s) ", &
          & A, ".")
! ======================================================================
  end function mnmacro
