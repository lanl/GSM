
  function nmina ( molObj, i )

! ====================================================================
!
! This function returns the value of 'nmina' from the data module
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(Molnix),  intent(inout) :: molObj
    integer(int32), intent(in   ) :: i
    real(real64)                  :: nmina

    integer(int32) :: usedI

! ====================================================================

    ! Verify index:
    if ( i < 1 ) then
       usedI = 1
       write(molObj%io%message, 1000) i, usedI
       call molObj%io%print(3, 3, molObj%io%message)
    else if ( i > 7 ) then
       usedI = 7
       write(molObj%io%message, 1000) i, usedI
       call molObj%io%print(3, 3, molObj%io%message)
    else
       usedI = i
    end if

    nmina = nminaData( usedI )
    return
! ====================================================================
1000 format("Element ", i4, " does not exist. Using element ", i4, &
          & " instead.")
! ====================================================================
  end function nmina


  function nmaxa ( molObj, i )

! ====================================================================
!
! This function returns the value of 'nmaxa' from the data module
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(Molnix),  intent(inout) :: molObj
    integer(int32), intent(in   ) :: i
    real(real64)                  :: nmaxa

    integer(int32) :: usedI

! ====================================================================

    ! Verify index:
    if ( i < 1 ) then
       usedI = 1
       write(molObj%io%message, 1000) i, usedI
       call molObj%io%print(3, 3, molObj%io%message)
    else if ( i > 7 ) then
       usedI = 7
       write(molObj%io%message, 1000) i, usedI
       call molObj%io%print(3, 3, molObj%io%message)
    else
       usedI = i
    end if

    nmaxa = nmaxaData( usedI )
    return
! ====================================================================
1000 format("Element ", i4, " does not exist. Using element ", i4, &
          & " instead.")
! ====================================================================
  end function nmaxa


  function nmin ( molObj, i )

! ====================================================================
!
! This function returns the value of 'nmin' from the data module
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(Molnix),  intent(inout) :: molObj
    integer(int32), intent(in   ) :: i
    real(real64)                  :: nmin

    integer(int32) :: usedI

! ====================================================================

    ! Verify index:
    if ( i < 1 ) then
       usedI = 1
       write(molObj%io%message, 1000) i, usedI
       call molObj%io%print(3, 3, molObj%io%message)
    else if ( i > 93 ) then
       usedI = 93
       write(molObj%io%message, 1000) i, usedI
       call molObj%io%print(3, 3, molObj%io%message)
    else
       usedI = i
    end if

    nmin = nminData( usedI )
    return
! ====================================================================
1000 format("Element ", i4, " does not exist. Using element ", i4, &
          & " instead.")
! ====================================================================
  end function nmin


  function nmax ( molObj, i )

! ====================================================================
!
! This function returns the value of 'nmax' from the data module
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(Molnix),  intent(inout) :: molObj
    integer(int32), intent(in   ) :: i
    real(real64)                  :: nmax

    integer(int32) :: usedI

! ====================================================================

    ! Verify index:
    if ( i < 1 ) then
       usedI = 1
       write(molObj%io%message, 1000) i, usedI
       call molObj%io%print(3, 3, molObj%io%message)
    else if ( i > 93 ) then
       usedI = 93
       write(molObj%io%message, 1000) i, usedI
       call molObj%io%print(3, 3, molObj%io%message)
    else
       usedI = i
    end if

    nmax = nmaxData( usedI )
    return
! ====================================================================
1000 format("Element ", i4, " does not exist. Using element ", i4, &
          & " instead.")
! ====================================================================
  end function nmax


  function properlyConstructed(molObj) result(constructed)

! ====================================================================
!
! Returns whether or not the object was constructed
!
! ====================================================================

    implicit none
    class(Molnix), intent(in   ) :: molObj
    logical :: constructed

! ====================================================================

    constructed = molObj%constructed

    return
! ====================================================================
  end function properlyConstructed


  function queryOptions(molObj) result(options)

! ====================================================================
!
! Returns the options used by the Molnix object
!
! ====================================================================

    implicit none
    class(Molnix), intent(in   ) :: molObj
    type(molnixOptions) :: options

! ====================================================================

    options = molObj%options

    return
! ====================================================================
  end function queryOptions
