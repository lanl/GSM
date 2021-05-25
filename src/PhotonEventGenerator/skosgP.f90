
  subroutine skosgp (e, sigkosgp)
!  subroutine skosgp (photonEG, e, sigkosgp)

! ======================================================================
!
!    gamma + H cross section. Only called with a0 = 1!
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64 

    implicit none
!    class(PhotonEventGenerator), intent(inout) :: photonEG   ! UNUSED!
    real(real64),   intent(in   ) :: e
    real(real64),   intent(  out) :: sigkosgp

    real(real64) :: ar, fp, fr, g4, g8, hpp, rdel, rh, s1, s2, s3, s4, &
         & temp, temp2, z, zud, zz1

! ======================================================================

    real(real64), parameter :: udelz = 5.813_real64
    real(real64), parameter :: wdelz = 0.056_real64
    real(real64), parameter :: whaz  = 0.045_real64

! ======================================================================

    z = log(e)
    zud = z - udelz
    temp2 = wdelz
    temp = one + zud*zud/temp2
    rdel = 0.55d0/(temp)
    zz1 = z - 6.57d0
    temp2 = whaz
    temp = one + zz1*zz1/temp2
    rh = 0.223d0/(temp)

    temp = one + exp(12.d0*(7.25d0 - z))
    g4 = exp(four*(6.27d0 - z))/(temp)

    temp = one + exp(24.d0*(6.9d0 - z))
    g8 = exp(eight*(6.66d0 - z))/(temp)

    temp = one + exp(four*(7.d0 - z))
    fp = one/(temp)
    hpp = 0.0375d0*(z - 16.5d0) + 1.07d0*exp(-0.11d0*z)

    ar = 25.d0*(5.24d0 - z)
    temp = one + exp(ar)
    fr = one/(temp)
    s1 = fr*(rdel + rh)
    s2 = g4
    s3 = g8
    s4 = fp*hpp

    sigkosgp = s1 + s2 + s3 + s4
    return
! ======================================================================
  end subroutine skosgp
