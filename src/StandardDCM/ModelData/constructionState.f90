
  function constructionState(sDCMDataObj) result(errorFlag)

! ====================================================================
!
! This function returns to the client an integer flag that indicates
! if the Standard DCM data class was properly constructed (if not,
! the value of the integer should indicate what error(s) occurred).
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    class(StandardDCMData), intent(in) :: sDCMDataObj
    integer(int32) :: errorFlag

! ====================================================================

    errorFlag = sDCMDataObj%constructed

    return
! ====================================================================
  end function constructionState
