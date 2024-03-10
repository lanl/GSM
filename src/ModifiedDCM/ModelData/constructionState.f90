
  function constructionState(mDCMDataObj) result(errorFlag)

! ====================================================================
!
! This function returns to the client an integer flag that indicates
! if the Modified DCM data class was properly constructed (if not,
! the value of the integer should indicate what error(s) occurred).
!
!
! Written by CMJ, XCP-3, 02/2024
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    class(ModifiedDCMData), intent(in) :: mDCMDataObj
    integer(int32) :: errorFlag

! ====================================================================

    errorFlag = mDCMDataObj%constructed

    return
! ====================================================================
  end function constructionState
