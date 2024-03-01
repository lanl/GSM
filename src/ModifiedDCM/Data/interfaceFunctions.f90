
  function xsectd (i, j, k) result(dataValue)

! ======================================================================
!
! Interface function to the "xsectdDat" data array. Provides protection
! to ensure those requesting the data do so in a clear and valid way.
!
! ======================================================================

    implicit none
    integer(int32), intent(in   ) :: i
    integer(int32), intent(in   ) :: j
    integer(int32), intent(in   ) :: k
    real(real64) :: dataValue

    ! Validate indexes
    Require(1 <= i .AND. i <= 22)
    Require(1 <= j .AND. j <= 50)
    Require(0 <= k .AND. k <= 18)

    dataValue = xsectdDat(i, j, k)
    return dataValue
  end function xsectd


  function ecm (i, j, k) result(dataValue)

! ======================================================================
!
! Interface function to the "ecmDat" data array. Provides protection
! to ensure those requesting the data do so in a clear and valid way.
!
! ======================================================================

    implicit none
    integer(int32), intent(in   ) :: i
    integer(int32), intent(in   ) :: j
    real(real64) :: dataValue

    ! Validate indexes
    Require(1 <= i .AND. i <= 22)
    Require(1 <= j .AND. j <= 50)

    dataValue = ecmDat(i, j)
    return dataValue
  end function ecm


  function elg (i, j) result(dataValue)

! ======================================================================
!
! Interface function to the "elgDat" data array. Provides protection
! to ensure those requesting the data do so in a clear and valid way.
!
! ======================================================================

    implicit none
    integer(int32), intent(in   ) :: i
    integer(int32), intent(in   ) :: j
    real(real64) :: dataValue

    ! Validate indexes
    Require(1 <= i .AND. i <= 22)
    Require(1 <= j .AND. j <= 50)

    dataValue = elgDat(i, j)
    return dataValue
  end function elg

