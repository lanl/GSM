
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
    use modifiedDCMParams, only: two, emnucg

    implicit none
    character(len=*), intent(in   ) :: photonFile

    integer(int32):: photonUnit
    integer(int32) :: inth, inw, jch, rc

! ======================================================================

    ! Allocate memory for mDCM data
    if (.not.allocated(xsectdDat)) allocate(xsectdDat(22, 50, 0:18))
    if (.not.allocated(ecmDat)) allocate(ecmDat(22, 50))
    if (.not.allocated(elgDat)) allocate(elgDat(22, 50))

    ! ***      Read differential cross section data file
    open(newunit = photonUnit, &
        & file=photonFile, &
        & status="old", &
        & action="read", &
        iostat = rc)
    Insist (rc == 0, "Failed to read file: " // photonFile)
    do jch = 1,22
       read (photonUnit, 40, iostat = rc)
       call insist(rc == 0, &
           & "Failed during read of channel label in file: " // photonFile, &
           & __FILE__, &
           & __LINE__)
       do inw = 1,50
          read (photonUnit, 30, iostat = rc) ecmDat(jch,inw), (xsectdDat(jch,inw,inth), inth=0,4)
          call insist(rc == 0, &
              & "Failed to read ECM and first half of cross section data in file: " // photonFile, &
              & __FILE__, &
              & __LINE__)
          read (photonUnit, 20, iostat = rc) (xsectdDat(jch,inw,inth), inth=5,18)
          call insist(rc == 0, &
              & "Failed to read last half of cross section data in file: " // photonFile, &
              & __FILE__, &
              & __LINE__)

          ! Perform calc with data now
          elgDat(jch,inw) = (ecmDat(jch,inw)**2 - emnucg**2)/(two*emnucg)
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
