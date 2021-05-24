
  function totalLinearMomentum(residual) result(totMomentum)

! ====================================================================
!
! This procedure returns the total linear momentum of a residual nucleus
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(GSMResidual), intent(in   ) :: residual
    real(real64)                      :: totMomentum

! ====================================================================

    totMomentum = sqrt( residual%linearMom(1)**2 + &
         & residual%linearMom(2)**2 + residual%linearMom(3)**2 )

    return
! ====================================================================
  end function totalLinearMomentum


  function totalAngularMomentum(residual) result(totMomentum)

! ====================================================================
!
! This procedure returns the total angular momentum of a residual nucleus
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(GSMResidual), intent(in   ) :: residual
    real(real64)                      :: totMomentum

! ====================================================================

    totMomentum = sqrt( residual%angularMom(1)**2 + &
         & residual%angularMom(2)**2 + residual%angularMom(3)**2 )

    return
! ====================================================================
  end function totalAngularMomentum
