
  function getMolnix( preeqData ) result(molObj)

! ====================================================================
!
! This function returns the molnix object used by the preequlibrium
! data object
!
! ====================================================================

    use molnixClass, only: Molnix

    implicit none
    class(PreequilibriumData), intent(in   ) :: preeqData
    class(Molnix),             pointer       :: molObj

! ====================================================================

    molObj => preeqData%molEnergy
    return
! ====================================================================
  end function getMolnix


  function getFissionBarrier( preeqData ) result(fbObj)

! ====================================================================
!
! This function returns the molnix object used by the preequlibrium
! data object
!
! ====================================================================

    use fissionBarrierClass, only: FissionBarrier

    implicit none
    class(PreequilibriumData), intent(in   ) :: preeqData
    class(FissionBarrier),     pointer       :: fbObj

! ====================================================================

    fbObj => preeqData%fissBarr
    return
! ====================================================================
  end function getFissionBarrier


  function properlyConstructed ( preeqData ) result(constructed)

! ====================================================================
!
! Returns a logical flag indicating if the class was constructed or not
!
! ====================================================================

    implicit none
    class(PreequilibriumData), intent(in   ) :: preeqData
    logical                                  :: constructed

! ====================================================================

    constructed = preeqData%constructed
    return
! ====================================================================
  end function properlyConstructed


  subroutine checkIndex(preeqData, i, k, j) 

! ====================================================================
!
! This subroutine checks that the passed in indices exist for the
! eb/egs arrays
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    class(PreequilibriumData), intent(inout) :: preeqData
    integer(int32),            intent(inout) :: i
    integer(int32),            intent(inout) :: j
    integer(int32),            intent(inout) :: k

! ====================================================================

    ! Check index "i"
    if ( i < 1 ) then
       write(preeqData%io%message, 1000) i, 1
       call preeqData%io%print(3, 3, preeqData%io%message)
       i = 1
    else if ( i > dataEgsEbDim1 ) then
       write(preeqData%io%message, 1000) i, dataEgsEbDim1
       call preeqData%io%print(3, 3, preeqData%io%message)
       i = dataEgsEbDim1
    end if


    ! Check index "j"
    if ( j < 1 ) then
       write(preeqData%io%message, 1000) j, 1
       call preeqData%io%print(3, 3, preeqData%io%message)
       j = 1
    else if ( j > dataEgsEbDim2 ) then
       write(preeqData%io%message, 1000) j, dataEgsEbDim2
       call preeqData%io%print(3, 3, preeqData%io%message)
       j = dataEgsEbDim1
    end if


    ! Check index "k"
    if ( k < 1 ) then
       write(preeqData%io%message, 1000) k, 1
       call preeqData%io%print(3, 3, preeqData%io%message)
       k = 1
    else if ( k > dataEgsEbDim3 ) then
       write(preeqData%io%message, 1000) k, dataEgsEbDim1
       call preeqData%io%print(3, 3, preeqData%io%message)
       k = dataEgsEbDim1
    end if

    return
! ====================================================================
1000 format("Invalid element (", i3, ") in preeq. data array.", &
          & " Using element ", i3, ".")
! ====================================================================
  end subroutine checkIndex


  function eb(preeqData, i, j, k)

! ====================================================================
!
! Returns the value of "eb(i,j,k)" in the case of valid indices
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(PreequilibriumData), intent(inout) :: preeqData
    integer(int32),            intent(in   ) :: i
    integer(int32),            intent(in   ) :: j
    integer(int32),            intent(in   ) :: k
    real(real64)                             :: eb

    integer(int32) :: usedI, usedJ, usedK

! ====================================================================

    ! Verify a valid element of the array was chosen:
    usedI = i
    usedJ = j
    usedK = k
    call preeqData%checkIndex(usedI, usedJ, usedK)

    ! Obtain element:
    eb = preeqData%compEnergy%eb(usedI, usedJ, usedK)
    return
! ====================================================================
  end function eb


  function egs(preeqData, i, j, k)

! ====================================================================
!
! Returns the value of "egs(i,j,k)" in the case of valid indices
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(PreequilibriumData), intent(inout) :: preeqData
    integer(int32),            intent(in   ) :: i
    integer(int32),            intent(in   ) :: j
    integer(int32),            intent(in   ) :: k
    real(real64)                             :: egs

    integer(int32) :: usedI, usedJ, usedK

! ====================================================================

    ! Verify a valid element of the array was chosen:
    usedI = i
    usedJ = j
    usedK = k
    call preeqData%checkIndex(usedI, usedJ, usedK)

    ! Obtain element:
    egs = preeqData%compEnergy%egs(usedI, usedJ, usedK)
    return
! ====================================================================
  end function egs


  function numBaryons ( preeqData )

! ====================================================================
!
! Returns the number of baryons in the compound nucleus at its creation
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(PreequilibriumData), intent(in   ) :: preeqData
    real(real64)                             :: numBaryons

! ====================================================================

    numBaryons = preeqData%compound%numBaryons
    return
! ====================================================================
  end function numBaryons


  function numProtons ( preeqData )

! ====================================================================
!
! Returns the number of baryons in the compound nucleus at its creation
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(PreequilibriumData), intent(in   ) :: preeqData
    real(real64)                             :: numProtons

! ====================================================================

    numProtons = preeqData%compound%numProtons
    return
! ====================================================================
  end function numProtons


  function kinEnergy ( preeqData )

! ====================================================================
!
! Returns the number of baryons in the compound nucleus at its creation
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: real64

    implicit none
    class(PreequilibriumData), intent(in   ) :: preeqData
    real(real64)                             :: kinEnergy

! ====================================================================

    kinEnergy = preeqData%compound%kinEnergy
    return
! ====================================================================
  end function kinEnergy
