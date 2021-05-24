
  subroutine readPhotonData(photonFile)

! ======================================================================
!
!    Main routine to extract ds/do for channel 1-22:
!
!    Written by K. K. Gudima, Fall 2003?
!    Modified by AJS, May, 2004.
!    Modified by KKG, Nov., 2004
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use modifiedDCMParams, only: zro, two, twpi, degreeToRad, emnucg

    implicit none
    character(len=*), intent(in   ) :: photonFile

    integer(int32):: photonUnit
    integer(int32) :: inth, inw, jch

! ======================================================================

    ! Allocate memory for mDCM data
    if (.not.allocated(xsectd)) allocate(xsectd(22, 50, 0:18))
    if (.not.allocated(ecm)) allocate(ecm(22, 50))
    if (.not.allocated(elg)) allocate(elg(22, 50))

    ! ***      Read differential cross section data file
    open(newunit = photonUnit, file=photonFile, status="old", action="read")
    do jch = 1,22
       read (photonUnit, 40)
       do inw = 1,50
          read (photonUnit, 30) ecm(jch,inw), (xsectd(jch,inw,inth), inth=0,4)
          read (photonUnit, 20) (xsectd(jch,inw,inth), inth=5,18)
          elg(jch,inw) = (ecm(jch,inw)**2 - emnucg**2)/(two*emnucg)
       end do
    end do
    close (photonUnit)

    return
! ======================================================================
20  format (7e10.3)
30  format (7x,f5.2,8x,5e10.3)
40  format (1x)
! ======================================================================
  end subroutine readPhotonData
