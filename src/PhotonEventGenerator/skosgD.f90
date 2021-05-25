
  subroutine skosgd (photonEG, e, sigkosgd)

! ======================================================================
!
!    gamma + d cross section.
!    Only called when a0 = 2 !!
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64 

    implicit none
    class(PhotonEventGenerator), intent(inout) :: photonEG
    real(real64),   intent(in   ) :: e
    real(real64),   intent(  out) :: sigkosgd

    real(real64) :: aln, ar, fp, fr, g1, g2, g4, g8, rdel, rh, s1, s2, &
         & s3, s4, sh, sp, sr, temp, temp2, tra, z, zs1, zud

! ======================================================================

    z = log(e)

    temp = one + exp(thr*(1.2d0 - z))
    g1 = exp(one*(1.86d0 - z))/(temp)
    g2 = exp(two*(2.11d0 - z))/(temp)

    temp = one + exp(thr*four*(7.1d0 - z))
    g4 = exp(four*(6.2d0 - z))/(temp)

    temp = one + exp(24.d0*(6.91d0 - z))
    g8 = exp(8.d0*(6.62d0 - z))/(temp)

    tra = 5.13d0 - 0.00150d0
    zud = z - udel(two)
    temp2 = wdel(two)
    if (temp2 < divZerLim .and. temp2 > -divZerLim) then
       temp2 = divZerLim
       write(photonEG%io%message,1000) "284"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    temp = one + zud*zud/temp2
    if (temp < divZerLim .and. temp > -divZerLim) then
       temp = divZerLim
       write(photonEG%io%message,1000) "289"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    rdel = 0.88d0/(temp)
    zs1 = z - 6.575d0
    temp2 = wha(two)
    if (temp2 < divZerLim .and. temp2 > -divZerLim) then
       temp2 = divZerLim
       write(photonEG%io%message,1000) "295"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    temp = one + zs1*zs1/temp2
    if (temp < divZerLim .and. temp > -divZerLim) then
       temp = divZerLim
       write(photonEG%io%message,1000) "300"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    rh = 0.348d0/(temp)

    ar = 25.d0*(tra - z)
    temp = one + exp(ar)
    fr = one/(temp)
    temp = one + exp(four *(7.d0 - z))

    fp = one/(temp)
    sr = fr*(rdel + rh)
    s1 = g1
    s2 = g2
    s3 = g4
    s4 = g8
    aln = log(two)
    sp = two*(1.d0 - 0.072d0*aln)
    sh = fp*sp*hpa (two, e)

    sigkosgd = sr + s1 + s2 + s3 + s4 + sh

    return

! ======================================================================
1000 format("Divide by zero error prevented in ", &
          & "'photoEventGenerator.f90' line(s) ", A)
! ======================================================================
  end subroutine skosgd
