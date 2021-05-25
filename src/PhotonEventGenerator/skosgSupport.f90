
  function hpa (a, e)

! ======================================================================
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    real(real64),   intent(in   ) :: a
    real(real64),   intent(in   ) :: e
    real(real64)                  :: hpa

    real(real64) :: aln, sha, z

! ======================================================================

    aln = log(a)
    sha = 1.0663d0 - 0.0023d0*aln
    z = log(e)
    hpa = 0.0375d0*(z - 16.5d0) + sha*exp(-0.11d0*z)
    return

! ======================================================================
  end function hpa


  function udel (a)

! ======================================================================
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    real(real64),   intent(in   ) :: a
    real(real64)                  :: udel

    real(real64) :: temp

! ======================================================================

    temp = 1.d0 + 0.003d0*a*a
    udel = 5.82d0 - 0.07d0/(temp)
    return

! ======================================================================
  end function udel


  function wdel (a)

! ======================================================================
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    real(real64),   intent(in   ) :: a
    real(real64)                  :: wdel

    real(real64) :: aln

! ======================================================================

    aln = log(a)
    wdel = 0.056d0 + aln*(0.03d0 - 0.001d0*aln)
    return

! ======================================================================
  end function wdel


  function wha (a)

! ======================================================================
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    real(real64),   intent(in   ) :: a
    real(real64)                  :: wha

    real(real64) :: ala

! ======================================================================

    ala = log(a)
    wha = 0.045d0 + 0.04d0*ala*sqrt(abs(ala))
    return

! ======================================================================
  end function wha
