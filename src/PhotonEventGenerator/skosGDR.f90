
  subroutine skosgdr (photonEG, e, a0, sigkosgdr)

! ======================================================================
!
!   Giant dipole resonance cross section.
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
    real(real64),   intent(  out) :: sigkosgdr

    real(real64) :: aln, g1, g2, g4,g8, r1, r2, r4, r8, s2, t1, t2, &
         & t4, t8, temp, z

! ======================================================================

    z = log(e)
    aln = log(a0)
    temp = a0
    if (temp < divZerLim .and. temp > -divZerLim) then
       temp = divZerLim
       write(photonEG%io%message,1000) "112"
       call photonEG%io%print(4, 3, photonEG%io%message)
    end if
    s2 = one + (two/temp)**4
    r1 = (3.2d0 + 0.75d0*aln)/s2
    r2 = (four   + 0.125d0*aln)/s2
    r4 = 3.8d0  + 0.05d0*aln
    r8 = 3.65d0 - 0.05d0*aln
    t1 = (6.6d0 - 0.5d0*aln)/s2
    t2 = 3.4d0/s2
    t4 = 3.8d0  - 0.25d0*aln
    t8 = 3.5d0  - 0.16d0*aln

    temp = one + exp(thr*one*(t1 - z))
    g1 = exp(one*(r1 - z))/(temp)

    temp = one + exp(thr*two*(t2 - z))
    g2 = exp(two*(r2 - z))/(temp)

    temp = one + exp(thr*four*(t4 - z))
    g4 = exp(four*(r4 - z))/(temp)

    temp = one + exp(thr*two*four*(t8 - z))
    g8 = exp(two*four*(r8 - z))/(temp)

    sigkosgdr = g1 + g2 + g4 + g8
    return

! ======================================================================
1000 format("Divide by zero error prevented in ", &
          & "'photoEventGenerator.f90' line(s) ", A)
! ======================================================================
  end subroutine skosgdr
