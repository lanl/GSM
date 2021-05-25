
  function sampleEnergy( gsmObj, centralValue, stdDev ) result(energy)

! ====================================================================
!
! Samples an energy around a central energy value (assumed normal distrbution)
! NOTE: Direct sampling is used!
!
! Formula from slide 93 of LA-UR-16-29043 (see normal dist.)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: real64
    use gsm_params, only: two, twpi

    implicit none
    class(GSM),   intent(in   ) :: gsmObj
    real(real64), intent(in   ) :: centralValue
    real(real64), intent(in   ) :: stdDev
    real(real64) :: energy

! ====================================================================

    real(real64), parameter :: sqrtTwPi = sqrt( twpi )

! ====================================================================

    ! Sample around some width, if desired:
    if( gsmObj%options%smoothTransition ) then
       energy = centralValue + stdDev * &
            & sqrt( -two * log( gsmObj%rang() ) ) * &
            & cos( twpi * gsmObj%rang() )
    else
       energy = centralValue
    end if

    return
! ====================================================================
  end function sampleEnergy
