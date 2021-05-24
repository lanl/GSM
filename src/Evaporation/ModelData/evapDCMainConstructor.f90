
  function evapDCMainConstructor(clientOptions, clientIO) result(evapData)

! ======================================================================
!
! Returns to user an "EvaporationData" object to be used by the Evaporation
! class based on calculation options.
!
!
! Written by CMJ, XCP-3, 8/2018 (Evap Class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, two
    use evaporationFissionData, only: aMax, maxSize, iiz, inn, &
         & omega, paire0, ux140, ux0

    implicit none
    type(evaporationDataOptions), intent(in   ), optional  :: clientOptions
    procedure(IOHANDLER), intent(in   ), optional, pointer :: clientIO
    type(EvaporationData)                         :: evapData

    integer(int32) :: iz0, in0, ia0, isdum
    real(real64)   :: ax140, ex0

! ======================================================================

    ! Set where messages get printed to:
    if ( present(clientIO) ) evapData%io%print => clientIO


    ! Set level density parameter option
    if ( present(clientOptions) ) then
       evapData%options = clientOptions
       call evapData%validateOptions()
    end if


    ! Establish arrays based on (Z,N)
    do iz0 = 1,iiz
       do in0 = 1,inn

          ia0 = iz0 + in0
          ex0 = ux0(ia0) + paire0(iz0,in0)

          ! Set ax0
          evapData%data%ax0(iz0,in0) = geta (evapData%options%alev, ex0, iz0, in0, isdum)
          evapData%data%ax0(iz0,in0) = max(zro,evapData%data%ax0(iz0,in0))

          ! Set taux0
          evapData%data%taux0(iz0,in0) = sqrt(evapData%data%ax0(iz0,in0)/ux0(ia0)) - &
               & 1.5d0/ux0(ia0)

          ! Set tx0, extx0
          evapData%data%tx0(iz0,in0) = ex0*evapData%data%taux0(iz0,in0)
          evapData%data%extx0(iz0,in0) = exp(evapData%data%tx0(iz0,in0))

          ! Set sx0
          ax140 = sqrt(sqrt(evapData%data%ax0(iz0,in0)))
          evapData%data%sx0(iz0,in0) = two*(ax140*ux140(ia0))**2

          ! Set fact10
          evapData%data%fact10(iz0,in0) = zro
          if (ax140 > zro .and. evapData%data%taux0(iz0,in0) > zro) then
             evapData%data%fact10(iz0,in0) = one/( evapData%data%taux0(iz0,in0)*ux0(ia0) * &
                  & ux140(ia0)*ax140 )
          end if
          evapData%data%fact10(iz0,in0) = evapData%data%fact10(iz0,in0) / &
               & ( evapData%data%extx0(iz0,in0) )

       end do
    end do

    ! Flag that data object was successfully created/established
    evapData%constructed = .TRUE.


    return
! ======================================================================
  end function evapDCMainConstructor
