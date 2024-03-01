
  subroutine inigam0( egamma )

! ======================================================================
!
!    Main routine to extract ds/do for channel 1-22:
!
!    Written by K. K. Gudima, Fall 2003?
!    Modified by AJS, May, 2004.
!    Modified by KKG, Nov., 2004
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMParams, only: zro, two, twpi, degreeToRad
    use modifiedDCMData, only: ctheta, xsectd, elg

    implicit none
    real(real64), intent(in   ) :: egamma

    real(real64), external:: qintxsq

    integer(int32) :: ie, ieg1, ieg2, j, jch
    real(real64)   :: eg1, eg2
    real(real64)   :: dom, s1, s2, sint, temp1
!    real(real64)   :: dom, eg1, eg2, s1, s2, sint, temp1

! ======================================================================

    real(real64), parameter :: &
         & dtheti =  1.0_real64

    real(real64), dimension(22, 19) :: s, r
    real(real64), dimension(22    ) :: st

! ======================================================================

    real(real64) :: thetai, cthetai, si, ri
    common /ixsgpn/ thetai(181),cthetai(181),si(22,182),ri(22,181)

! ======================================================================

    ! Establish theta bins:
    do j = 1,181
       thetai(j) = dble(j-1)*dtheti
       cthetai(j) = cos(thetai(j)*degreeToRad)
    enddo

    ! Load data now:
    do jch = 1,22
       if (egamma <= elg(jch,2))      then
          ieg1 = 2
          ieg2 = 3
       elseif (egamma >= elg(jch,50))  then
          ieg1 = 49
          ieg2 = 50
       else
          do ie = 3,50
             if (egamma >= elg(jch,ie-1) .and. egamma <= elg(jch,ie)) &
                  & then
                ieg1 = ie-1
                ieg2 = ie
                go to 10
             endif
          enddo
       endif
10     if (ieg1 < 2 .or. ieg1 > 49 .or. ieg2 < 3 .or. ieg2 > 50) &
            & then
          write (*,*) '  stop in inigam: ieg1, ieg2 = ', ieg1, ieg2
          stop
       endif

       eg1 = elg(jch,ieg1)
       eg2 = elg(jch,ieg2)
       temp1 = (eg2 - eg1)
       sint = zro
       do j = 1,19
          s1 = xsectd(jch,ieg1,j-1)
          s2 = xsectd(jch,ieg2,j-1)

          s(jch,j) = s1 + ((s2 - s1)/temp1)*(egamma - eg1)

          if (j >= 2) then
             dom = twpi*(ctheta(j-1) - ctheta(j))
             sint = sint + dom*(s(jch,j-1) + s(jch,j))/two
          endif
          r(jch,j) = sint
       enddo
       if(sint < zro)  sint = zro
       st(jch) = sint
       temp1 = sint
       do j = 1,19
          if (sint > zro)  r(jch,j) = r(jch,j)/temp1
       enddo
    end do
    do jch = 1,22
       sint = zro
       do j = 1,181
          si(jch,j) = qintxsq(thetai(j), s, jch, 22, 19)
          if (j >= 2) then
             dom = twpi*(cthetai(j-1) - cthetai(j))
             sint = sint + dom*(si(jch,j-1) + si(jch,j))/two
          endif
          ri(jch,j) = sint
       end do
       if(sint < zro)  sint = zro
       si(jch,182) = sint
       temp1 = sint
       do j = 1,181
          if (sint > zro)  ri(jch,j) = ri(jch,j)/temp1
       enddo
    end do

    return
! ======================================================================
  end subroutine inigam0
