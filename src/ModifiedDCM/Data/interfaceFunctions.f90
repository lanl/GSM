
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
    return
  end function xsectd


  function ecm (i, j) result(dataValue)

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
    return
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
    return
  end function elg


  function jamin (i) result(dataValue)

! ======================================================================
!
! Interface function to the "jaminDat" data array. Provides protection
! to ensure those requesting the data do so in a clear and valid way.
!
! ======================================================================

    implicit none
    integer(int32), intent(in   ) :: i
    integer(int32) :: dataValue

    ! Validate indexes
    Require(1 <= i .AND. i <= numNuclideLimits)

    dataValue = jaminDat(i)
    return
  end function jamin


  function jamax (i) result(dataValue)

! ======================================================================
!
! Interface function to the "jamaxDat" data array. Provides protection
! to ensure those requesting the data do so in a clear and valid way.
!
! ======================================================================

    implicit none
    integer(int32), intent(in   ) :: i
    integer(int32) :: dataValue

    ! Validate indexes
    Require(1 <= i .AND. i <= numNuclideLimits)

    dataValue = jamaxDat(i)
    return
  end function jamax


  function rms (i) result(dataValue)

! ======================================================================
!
! Interface function to the "rmsDat" data array. Provides protection
! to ensure those requesting the data do so in a clear and valid way.
!
! ======================================================================

    implicit none
    integer(int32), intent(in   ) :: i
    real(real64) :: dataValue

    ! Validate indexes
    Require(1 <= i .AND. i <= 10)

    dataValue = rmsDat(i)
    return
  end function rms


  function theta (i) result(dataValue)

! ======================================================================
!
! Interface function to the "thetaDat" data array. Provides protection
! to ensure those requesting the data do so in a clear and valid way.
!
! ======================================================================

    implicit none
    integer(int32), intent(in   ) :: i
    real(real64) :: dataValue

    ! Validate indexes
    Require(1 <= i .AND. i <= numThetaBins)

    dataValue = thetaDat(i)
    return
  end function theta


  function ctheta (i) result(dataValue)

! ======================================================================
!
! Interface function to the "cthetaDat" data array. Provides protection
! to ensure those requesting the data do so in a clear and valid way.
!
! ======================================================================

    implicit none
    integer(int32), intent(in   ) :: i
    real(real64) :: dataValue

    ! Validate indexes
    Require(1 <= i .AND. i <= numThetaBins)

    dataValue = cthetaDat(i)
    return
  end function ctheta



