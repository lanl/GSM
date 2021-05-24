
  function wopt (sDCM, clientTarg, e1, r, i1)

! ======================================================================
!
!    Negative of the imaginary part of optical potential.
!    Wopt is the negative of the W(r) defined in Becchetti & Greenlees.
!
!   Called by CASCAD
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February-March, 1996.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two, four

    implicit none
    class(StandardDCM),     intent(inout) :: sDCM
    class(StandardDCMData), intent(inout) :: clientTarg
    real(real64),           intent(in   ) :: e1
    real(real64),           intent(in   ) :: r
    integer(int32),         intent(in   ) :: i1
    real(real64)                          :: wopt

    real(real64) :: am, e, f, rm, temp, wsf, wv, x

! ======================================================================

    x = (clientTarg%numBaryons() - two*clientTarg%numProtons()) / &
         & clientTarg%numBaryons()
    e = e1*1000.d0
    if (i1 <= 0) then
       if (e <= 25.d0) then
!  Imaginary part of Becchetti & Greenlees neutron optical potential.
!  [Phys. Rev. 182, 1190 (1969)] 0 < E < 40 MeV.
          rm = 1.26d0*clientTarg%aTargThrd()
          am = 0.58d0
          wv = 0.22d0*e - 1.56d0
          wsf = 13.d0 - 0.25d0*e - 12.d0*x
       else
!  Imaginary part of Marshak, Langford, Tamura, & Wong neutron optical
!  potential [Phys. Rev. C 2, 1862 (1970)]. 2 < E < 129 MeV.
          rm = 1.21d0*clientTarg%aTargThrd()
          am = 0.6448d0
          wv = 0.459d0 + 0.111d0*e
          wsf = 4.28d0 - 0.0414d0*e
       endif
    else
       if (e <= 25.d0) then
!  Imaginary part of Becchetti & Greenlees proton optical potential.
          rm = 1.32d0*clientTarg%aTargThrd()
          am = 0.51d0 + 0.7d0*x
          wv = 0.22d0*e - 2.7d0
          wsf = 11.8d0 - 0.25d0*e + 12.d0*x
       else
!  Imaginary part of Menet, Gross, Malanify, & Zucker proton optical
!  potential [Phys. Rev. C 4, 1114 (1971)]  30 < E < 60 MeV.
          rm = 1.37d0*clientTarg%aTargThrd()
          am = 0.74d0 - 0.008d0*e + x
          wv = 1.2d0 + 0.09d0*e
          wsf = 4.2d0 - 0.05d0*e + 15.5d0*x
       endif
    endif
    wv = max (wv, zro)
    wsf = max (wsf, zro)
    if (am < div0Lim .and. am > -div0Lim) then
       am = div0Lim
       write(sDCM%io%message,1000) "76"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    temp = one + exp((r - rm)/am)
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(sDCM%io%message,1000) "81"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    f = one/(temp)
!  df/dr = (f - 1)/(am*f)
    wopt = f*(wv + four*wsf*(one - f))
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'wopt.f90' line(s) ", A)
! ======================================================================
  end function wopt
