
  function cosel (sDCM, l, mb, ksi, t, cm, photoData)

! ======================================================================
!
!     Cosine of a random angle, weighted by the angular distribution
!     for elastic scattering.
!     mb is baryon number; l is hardwired to be 0.
!     Reinstated photon logic (l /=/ 0) AJS  (10/10/03)
!     Call to COSTAN, which removes unphysical features near 0 and pi
!     was added by K. Gudima, 2/04
!     Exact kinematic limit tmax inserted in place of approximate one
!     KKG, 2/04.
!
!   Called by: ELEX
!
!   Calls: COSGAMN COSTA COSTAN
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk  LANL  T-2  February, 1996.
!   Edited by AJS, July, 1997.
!   Modified by AJS, January, 1999.
!   "Last" change: 12-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by K. K. Gudima, Feb., 2004.
!    Modified by A. J. Sierk, LANL T-16, March, 2004.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: hlf, one, two, thr, emnucg

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    integer(int32), intent(in   ) :: l
    integer(int32), intent(in   ) :: mb
    integer(int32), intent(in   ) :: ksi
    real(real64),   intent(in   ) :: t
    real(real64),   intent(in   ) :: cm
    type(sDCMPhotonCrossSections), intent(in   ) :: photoData
    real(real64)                  :: cosel

    integer(int32) :: ik, jdel
    real(real64)   :: r1, tcm, temp, tm7, tm8, tmax

! ======================================================================

    r1 = sDCM%rang()
    if (l.ne.0) then
!  Photon reactions:
       if (t <= 0.45d0) then
!  KKG 02/12/04
          cosel = sDCM%cosgamn (12, r1, photoData)
       else
          cosel = sDCM%cosgamn (13, r1, photoData)
       endif
    else
       if (mb >= 2) then
!  N-N scattering:
          if (t <= 2.8d0) then
             cosel = (one + sDCM%costan (1, t, r1))/two
          elseif (t <= 10.d0) then
             cosel = 0.25d0*(thr - sDCM%costan (2, t, r1))
          else
!  KKG 02/16/04:
             tmax = two*t*cm
             tm8 = 8.7d0*tmax
             temp = tm8
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(sDCM%io%message,1000) "73"
                call sDCM%io%print(4, 3, sDCM%io%message)
             end if
             cosel = one + (two*log(one + r1*(exp(-tm8) - one)))/temp
          endif
       else
!  pi-N scattering:
          if (ksi < 2) then
!  pi+ p or pi- n scattering:
             jdel = 0
          elseif (ksi == 2) then
!  pi+ n or pi- p scattering:
             jdel = 4
          else
!  pi0 p or pi0 n scattering:
             if (sDCM%rang() <= hlf) then
                jdel = 0
             else
                jdel = 4
             endif
          endif
          if (t <= 0.08d0) then
             ik = 4
          elseif (t <= 0.3d0) then
             ik = 5
          elseif (t <= one) then
             ik = 6
          elseif (t <= 2.4d0) then
             ik = 7
          endif
          if (t <= 2.4d0) then
             cosel = sDCM%costan (ik+jdel, t, r1)
          else
!  KKG 02/16/04 (exact definition ot tmax):
             tcm = two*emnucg
             temp = tcm*t + (emnucg + cm)**2
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(sDCM%io%message,1000) "110"
                call sDCM%io%print(4, 3, sDCM%io%message)
             end if
             tmax = tcm*t*tcm*(t + two*cm)/(temp)
             tm7 = 7.5d0*tmax
             temp = tm7
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(sDCM%io%message,1000) "117"
                call sDCM%io%print(4, 3, sDCM%io%message)
             end if
             cosel = one + (two*log(one + r1*(exp(-tm7) - one)))/temp
          endif
       endif
    endif

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'cosel.f90' line(s) ", A)
! ======================================================================
  end function cosel
