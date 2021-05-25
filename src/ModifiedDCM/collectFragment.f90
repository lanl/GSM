
  subroutine collectFragment(kCurrent)

! ==============================================================================
!
! Subroutine meant to:
!    (1) Check if CEM/GSM array is exceeded
!    (2) Check if the fragment is allowed (basic check)
!    (3) If all tests pass, collect fragment data into CEM/GSM arrays
! 
! Written by CMJ, XCP-3 (5/2018)
! 
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    integer(int32), intent(in   ) :: kCurrent

    integer(int32) :: indx
    integer(int32) :: im(9)
    real(real64)   :: pm(9)
    logical        :: validPart

    real(real64)   :: pmemoLAQ
    integer(int32) :: imemoLAQ
    common /memorylaq/ pmemoLAQ(9,5999), imemoLAQ(5,5999)

! ==============================================================================

    ! Obtain the current particle to be tallied:
    do indx = 1, 9
       pm(indx) = pmemoLAQ(indx, kCurrent)
    enddo
    do indx = 1, 5
       im(indx) = imemoLAQ(indx, kCurrent)
    enddo


    ! Check if fragment is valid, if so collect/tally into arrays
    call isValidParticle(im(4), im(1), pm(8), pm(9), validPart)


    ! Collect fragment if valid
    if ( validPart ) then
       call tallyMDCMProgeny(im, pm)
    end if

    return
! ==============================================================================
  end subroutine collectFragment
