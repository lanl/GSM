
  function validResidual(numBaryons, numProtons) &
       & result(isValid)

! ====================================================================
!
! This procedure verifies that a given residual is valid (based on its
! size)
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: real64
    use gsm_params, only: one, four

    implicit none
    real(real64),   intent(in   ) :: numBaryons
    real(real64),   intent(in   ) :: numProtons
    logical                       :: isValid

! ====================================================================

    isValid = .TRUE.

    if ( numBaryons < four .or. numProtons < one .or. &
         & numBaryons < numProtons ) then
       isValid = .FALSE.
    end if

    return
! ====================================================================
  end function validResidual
