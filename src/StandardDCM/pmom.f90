
  function pmom (j, t, r1)

! ======================================================================
!
!     Calculation of a secondary particle's random 3-momentum, weighted
!     by the momentum distribution.  See SR BD1 for significance of
!     index j.
!
!   Called by: VMNSP
!
!   CEM95 written by S. G. Mashnik
!
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Edited by AJS, July, 1997.
!    Modified by AJS, February, 1999.
!   "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by KKG, February, 2004.
!    Modified by A. J. Sierk, LANL T-16, March, 2004.
!    Modified by A. J. Sierk, LANL T-16, January, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one
    use standardDCMData,  only: bnkj, ckj

    implicit none
    integer(int32), intent(in   ) :: j
    real(real64),   intent(in   ) :: t
    real(real64),   intent(in   ) :: r1
    real(real64)                  :: pmom

    integer(int32) :: k, n
    real(real64)   :: pmax, s1, s2, term
    real(real64)   :: rn(5), tk(4)

! ======================================================================

    s1 = zro
    s2 = zro
    pmax = zro
    rn(1) = one
    rn(2) = r1
    rn(3) = r1*r1
    rn(4) = r1*r1*r1
    rn(5) = rn(3)*rn(3)
    tk(1) = one
    tk(2) = t
    tk(3) = t*t
    tk(4) = t*t*t
    do n = 1,4
       do k = 1,4
          term = bnkj(n,k,j)*tk(k)
          s1 = s1 + term*rn(n)
          s2 = s2 + term
       end do
    end do
    do k = 1,3
       pmax = pmax + ckj(k,j)*tk(k)
    end do
    pmom = pmax*sqrt(abs(r1))*(s1 + (one - s2)*rn(5))
    return

! ====================================================================
  end function pmom
