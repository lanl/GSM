
  function qintxs (photonEG, x, sig, l, m, n)

! ======================================================================
!
!    Interpolation of gamma+N differential cross sections.
!
!    Written by K. K. Gudima, Fall 2003?
!    Modified by AJS, May, 2004.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none
    class(PhotonEventGenerator), intent(inout) :: photonEG
    integer(int32),              intent(in   ) :: m
    integer(int32),              intent(in   ) :: n
    real(real64),                intent(in   ) :: x
    real(real64),                intent(in   ) :: sig(m, n)
    integer(int32),              intent(in   ) :: l
    real(real64)                               :: qintxs

    integer(int32) :: k, k1, k2, k3, numElements
    real(real64)   :: a, b, c, d1, da, db, dc, y1, y2, y3

! ======================================================================

    ! Provide default value
    qintxs = 0


    ! Ensure course angle mesh data was created:
    if ( photonEG%data%numElements <= 0 ) then
       ! Wasn't constructed OR no course angle mesh exists
       ! Warn and return to client
       write(photonEG%io%message, 2000)
       call photonEG%io%print(2, 2, photonEG%io%message)
       write(photonEG%io%message, 2010)
       call photonEG%io%print(2, 2, photonEG%io%message)
       return
    end if

    ! Look for points very close to the desired value
    numElements = min(n, photonEG%data%numElements)
    do k = 1, numElements
       if (abs(x - photonEG%data%theta(k)) <= 1.0d-3) then
          qintxs = sig(l,k)
          qintxs = max(qintxs, zro)
          return
       endif
    enddo

    ! Perform interpolation
    do k = 2, (numElements-1)
       ! Look for a point nearby
       if (x < photonEG%data%theta(k)) then
          k1 = k - 1
          k2 = k
          k3 = k + 1
          go to 10
       endif
    enddo
    ! No points found - use maximum
    k1 = numElements - 2
    k2 = numElements - 1
    k3 = numElements


10  continue
    y1 = sig(l,k1)
    y2 = sig(l,k2)
    y3 = sig(l,k3)
    da = photonEG%data%x23(k2)*y1 + photonEG%data%x31(k2)*y2 + &
         & photonEG%data%x12(k2)*y3
    db =   (y2 - y3)*photonEG%data%x11(k2) + &
         & (y3 - y1)*photonEG%data%x22(k2) + &
         & (y1 - y2)*photonEG%data%x33(k2)
    dc =   (photonEG%data%theta(k2)*y3 - photonEG%data%theta(k3)*y2)*photonEG%data%x11(k2) + &
         & (photonEG%data%theta(k3)*y1 - photonEG%data%theta(k1)*y3)*photonEG%data%x22(k2) + &
         & (photonEG%data%theta(k1)*y2 - photonEG%data%theta(k2)*y1)*photonEG%data%x33(k2)
    d1 = zro
    if ( abs(photonEG%data%d(k2)-zro)>=divZerLim ) d1 = one/photonEG%data%d(k2)
    a = da*d1
    b = db*d1
    c = dc*d1


    qintxs = a*x*x + b*x + c
    qintxs = max(qintxs,zro)

    return

! ======================================================================
2000 format("No photon cross section interpolation data was created.")
2010 format("   Cannot interpolate cross section values in 'qintxs.f90'.")
! ======================================================================
  end function qintxs
