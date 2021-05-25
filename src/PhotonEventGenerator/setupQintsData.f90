
  subroutine setupQintsData ( photonEG, thetaValues, numElements )

! ======================================================================
!
! Sets up data for the "qintxs" function (i.e. angle bins for interpolation)
!
!
! Written by CMJ, XCP-3, 12/2018
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(PhotonEventGenerator), intent(inout) :: photonEG
    real(real64),                intent(in   ) :: thetaValues(:)
    integer(int32),              intent(in   ) :: numElements

    integer(int32) :: k
    integer(int32) :: usedElements = 0_int32

! ======================================================================
 
    ! Ensure array bounds are not exceeded based on user input
    if ( numElements > 0 ) then
       usedElements = min( numElements, size(thetaValues) )
    else
       usedElements = size( thetaValues )
    end if


    ! Allocate memory for the data
10  continue
    if ( .not. allocated(photonEG%data%x11) ) then
       allocate( photonEG%data%theta(usedElements) )
       allocate( photonEG%data%x11(usedElements) )
       allocate( photonEG%data%x22(usedElements) )
       allocate( photonEG%data%x33(usedElements) )
       allocate( photonEG%data%x12(usedElements) )
       allocate( photonEG%data%x23(usedElements) )
       allocate( photonEG%data%x31(usedElements) )
       allocate( photonEG%data%d(usedElements)   )
       photonEG%data%numElements = usedElements
    else
       ! Deallocate and then re-allocate (ensures correct size)
       deallocate( photonEG%data%theta )
       deallocate( photonEG%data%x11 )
       deallocate( photonEG%data%x22 )
       deallocate( photonEG%data%x33 )
       deallocate( photonEG%data%x12 )
       deallocate( photonEG%data%x23 )
       deallocate( photonEG%data%x31 )
       deallocate( photonEG%data%d )
       photonEG%data%numElements = 0
       go to 10
    end if


    ! Initialize data
    photonEG%data%theta(:) = thetaValues(:)
    photonEG%data%x11(:) = zro
    photonEG%data%x22(:) = zro
    photonEG%data%x33(:) = zro
    photonEG%data%x12(:) = zro
    photonEG%data%x23(:) = zro
    photonEG%data%x31(:) = zro
    do k = 2, (photonEG%data%numElements-1)
       photonEG%data%x11(k) = photonEG%data%theta(k-1)**2
       photonEG%data%x22(k) = photonEG%data%theta(k)**2
       photonEG%data%x33(k) = photonEG%data%theta(k+1)**2
       photonEG%data%x12(k) = photonEG%data%theta(k-1) - photonEG%data%theta(k)
       photonEG%data%x31(k) = photonEG%data%theta(k+1) - photonEG%data%theta(k-1)
       photonEG%data%x23(k) = photonEG%data%theta(k)   - photonEG%data%theta(k+1)
       photonEG%data%d(k)   = photonEG%data%x23(k)*photonEG%data%x11(k) + &
            & photonEG%data%x31(k)*photonEG%data%x22(k) + &
            & photonEG%data%x12(k)*photonEG%data%x33(k)
    end do

    return
! ======================================================================
  end subroutine setupQintsData
