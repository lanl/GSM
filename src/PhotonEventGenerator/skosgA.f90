
  subroutine skosga (photonEG, e, a0, sigkosga)

! ======================================================================
!
!    gamma + A cross section for Z > 2.
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
    real(real64),   intent(in   ) :: a0
    real(real64),   intent(  out) :: sigkosga

    real(real64) :: aln, fp, fr, rdel, rh, s0, sh, sigkosgdr, sp, &
         & temp, temp2, uha, z, zud, zuh

! ======================================================================

    z = log(e)
    aln = log(a0)
    zud = z - udel(a0)
    temp2 = wdel(a0)
    if (temp2 < divZerLim .and. temp2 > -divZerLim) then
       temp2 = divZerLim
       write(photonEG%io%message,1000) "185"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    temp = one + zud*zud/temp2
    if (temp < divZerLim .and. temp > -divZerLim) then
       temp = divZerLim
       write(photonEG%io%message,1000) "190"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    rdel = 0.39d0*a0/(temp)
    uha = 6.496d0 + 0.042d0*aln
    zuh = z - uha
    temp2 = wha(a0)
    if (temp2 < divZerLim .and. temp2 > -divZerLim) then
       temp2 = divZerLim
       write(photonEG%io%message,1000) "197"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    temp = one + zuh*zuh/temp2
    if (temp < divZerLim .and. temp > -divZerLim) then
       temp = divZerLim
       write(photonEG%io%message,1000) "205"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    if (aln < divZerLim .and. aln > -divZerLim) then
       aln = divZerLim
       write(photonEG%io%message,1000) "206"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    rh = (0.16d0*a0/sqrt(abs(aln)))/(temp)
    temp = one + exp(11.d0*((5.13d0 - 0.00075d0*a0) - z))
    fr = one/(temp)
    temp = one + exp(four*(7.d0 - z))

    fp = one/(temp)
    s0 = fr*(rdel + rh)
    sp = a0*(one  - 0.072d0*aln)
    sh = fp*sp*hpa (a0, e)

    call photonEG%skosgdr (e, a0, sigkosgdr)
    sigkosga = sigkosgdr + s0 + sh

    return

! ======================================================================
1000 format("Divide by zero error prevented in ", &
          & "'photoEventGenerator.f90' line(s) ", A)
! ======================================================================
  end subroutine skosga
