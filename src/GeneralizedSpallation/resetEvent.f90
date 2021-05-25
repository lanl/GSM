
! ======================================================================
!
! Resets the event-specific members of the results object
!
! ======================================================================

  subroutine resetEvent(results)

    use, intrinsic:: iso_fortran_env, only: int32
    implicit none
    class(GSMResults), intent(inout) :: results

    integer(int32) :: iter

! ======================================================================

    ! Reset progeny:
    results%numProgeny    = 0_int32
    results%progenyBnk(:) = GSMProgeny()

    ! Reset residuals and exciton data
    results%projRes = GSMResidual()
    results%targRes = GSMResidual()
    results%projExc = ExcitonData()
    results%targExc = ExcitonData()

    ! Reset misc.
    results%numElasticEvents = 0_int32
    results%modelUsage = GSMModelUsage()

    return
! ======================================================================
  end subroutine resetEvent
