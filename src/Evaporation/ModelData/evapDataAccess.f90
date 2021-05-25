
  subroutine checkIndex( evapData, iz, in )

! ====================================================================
!
! This procedure makes the iz/in index valid if it is not
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32
    use evaporationFissionData,       only: iiz, inn

    implicit none
    class(EvaporationData), intent(inout) :: evapData
    integer(int32),         intent(inout) :: iz
    integer(int32),         intent(inout) :: in

! ====================================================================

    ! Verify Z index:
    if ( iz <   1 ) then
       write(evapData%io%message, 1000) iz, 1
       call evapData%io%print(3, 3, evapData%io%message)
       iz = 1
    else if ( iz > iiz ) then
       write(evapData%io%message, 1000) iz, iiz
       call evapData%io%print(3, 3, evapData%io%message)
       iz = iiz
    end if


    ! Verify N index:
    if ( in <   1 ) then
       write(evapData%io%message, 2000) in, 1
       call evapData%io%print(3, 3, evapData%io%message)
       in = 1
    else if ( in > inn ) then
       write(evapData%io%message, 2000) in, inn
       call evapData%io%print(3, 3, evapData%io%message)
       in = inn
    end if

    return
! ====================================================================
1000 format("Invalid Z index for evaporation data (", i5, &
          & "). Approximating to ", i5, ".")
2000 format("Invalid N index for evaporation data (", i5, &
          & "). Approximating to ", i5, ".")
! ====================================================================
  end subroutine checkIndex


  function alev ( evapData )

! ====================================================================
!
! This function returns the value of 'alev' to the client
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(EvaporationData), intent(in   ) :: evapData
    real(real64)                          :: alev

! ====================================================================

    alev = evapData%options%alev
    return
! ====================================================================
  end function alev


  function fact10 ( evapData, iz, in )

! ====================================================================
!
! This function returns the value of 'fact10' to the client
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(EvaporationData), intent(inout) :: evapData
    integer(int32),         intent(in   ) :: iz
    integer(int32),         intent(in   ) :: in
    real(real64)                          :: fact10

    integer(int32) :: thisZ, thisN

! ====================================================================

    thisZ = iz
    thisN = in
    call evapData%checkIndex( thisZ, thisN )

    fact10 = evapData%data%fact10(thisZ, thisN)
    return
! ====================================================================
  end function fact10



  function sx0 ( evapData, iz, in )

! ====================================================================
!
! This function returns the value of 'sx0' to the client
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(EvaporationData), intent(inout) :: evapData
    integer(int32),         intent(in   ) :: iz
    integer(int32),         intent(in   ) :: in
    real(real64)                          :: sx0

    integer(int32) :: thisZ, thisN

! ====================================================================

    thisZ = iz
    thisN = in
    call evapData%checkIndex( thisZ, thisN )

    sx0 = evapData%data%sx0( thisZ, thisN )
    return
! ====================================================================
  end function sx0



  function taux0 ( evapData, iz, in )

! ====================================================================
!
! This function returns the value of 'taux0' to the client
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(EvaporationData), intent(inout) :: evapData
    integer(int32),         intent(in   ) :: iz
    integer(int32),         intent(in   ) :: in
    real(real64)                          :: taux0

    integer(int32) :: thisZ, thisN

! ====================================================================

    thisZ = iz
    thisN = in
    call evapData%checkIndex( thisZ, thisN )

    taux0 = evapData%data%taux0( thisZ, thisN )
    return
! ====================================================================
  end function taux0



  function ax0 ( evapData, iz, in )

! ====================================================================
!
! This function returns the value of 'ax0' to the client
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(EvaporationData), intent(inout) :: evapData
    integer(int32),         intent(in   ) :: iz
    integer(int32),         intent(in   ) :: in
    real(real64)                          :: ax0

    integer(int32) :: thisZ, thisN

! ====================================================================

    thisZ = iz
    thisN = in
    call evapData%checkIndex( thisZ, thisN )

    ax0 = evapData%data%ax0( thisZ, thisN )
    return
! ====================================================================
  end function ax0



  function tx0 ( evapData, iz, in )

! ====================================================================
!
! This function returns the value of 'tx0' to the client
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(EvaporationData), intent(inout) :: evapData
    integer(int32),         intent(in   ) :: iz
    integer(int32),         intent(in   ) :: in
    real(real64)                          :: tx0

    integer(int32) :: thisZ, thisN

! ====================================================================

    thisZ = iz
    thisN = in
    call evapData%checkIndex( thisZ, thisN )

    tx0 = evapData%data%tx0( thisZ, thisN )
    return
! ====================================================================
  end function tx0



  function extx0 ( evapData, iz, in )

! ====================================================================
!
! This function returns the value of 'extx0' to the client
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(EvaporationData), intent(inout) :: evapData
    integer(int32),         intent(in   ) :: iz
    integer(int32),         intent(in   ) :: in
    real(real64)                          :: extx0

    integer(int32) :: thisZ, thisN

! ====================================================================

    thisZ = iz
    thisN = in
    call evapData%checkIndex( thisZ, thisN )

    extx0 = evapData%data%extx0( thisZ, thisN )
    return
! ====================================================================
  end function extx0


  function queryOptions( evapData ) result(options)

! ====================================================================
!
! Returns the options contained by the object to the client
!
! ====================================================================

    implicit none
    class(EvaporationData), intent(in   ) :: evapData
    type(evaporationDataOptions) :: options

! ====================================================================

    options = evapData%options
    return
! ====================================================================
  end function queryOptions
