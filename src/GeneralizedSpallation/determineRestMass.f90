
  function determineRestMass(nuclA, nuclZ ) result(restMass)

! ====================================================================
!
! Determines the rest mass of a nucleus given its A/Z combination
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    integer(int32), intent(in   ) :: nuclA
    integer(int32), intent(in   ) :: nuclZ
    real(real64) :: restMass

! ====================================================================

    restMass = dble(nuclZ) * 0.9382723d0 + &
         & dble(nuclA - nuclZ) * 0.9395656

    return
! ====================================================================
  end function determineRestMass
