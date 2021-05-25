
  function initializeGammaData(desiredGammaFile) result(errorFlag)

! ======================================================================
!    Main routine to extract ds/do for channel 1-4 (single pion
!    production):
!
!    Written by K. K. Gudima, Fall 2003?
!    Modified by AJS, May, 2004.
!    Modified by KKG, Sep., 2004
!    Edited by AJS, January, 2005.
!    Modified by AJS, February, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64, output_unit
    use photonEventGeneratorClass, only: newPhotonEventGenerator
    use standardDCMParams, only: zro, one, two, pi, radianToDegree, emnucg

    implicit none
    character(LEN=*), intent(in   ) :: desiredGammaFile
    integer(int32)                  :: errorFlag

    integer(int32) :: i, inth, inw, ith, jch
    real(real64)   :: temp

    ! Set default file name and unit
    integer(int32)    :: gammaUnit
    character(LEN=64) :: gammaFile = defaultGammaFile

    ! Assume desired file unit is used and that the file does NOT exist
    logical        :: fileExists   = .false.

    ! (for when an error occurs during reading of the file)
    integer(int32), dimension(3) :: readError = 0_int32

! ======================================================================

    ! Data variables (transformed based on the angle bins established)
    ! Note: The below require ~31 kB of memory when allocated

    ! Cross section for gamma + d channel
    real(real64), dimension(:, :, :), allocatable :: xsectd

    ! Energy in CM system for channel?
    real(real64), dimension(:, :), allocatable :: ecm

! ======================================================================

    ! Assume no error occurs
    errorFlag = sDCMDataSetup


    ! Establish angle bins
    ! (Course angle bins)
    do j = 1,19
       theta(j) = dble(j-1) * dtheta
       ctheta(j) = cos( theta(j) * radianToDegree )
       if (j > 1) domo2(j) = pi * ( ctheta(j-1) - ctheta(j) )
    enddo

    ! (Fine angle bins)
    do j = 1,181
       thetai(j) = dble(j-1) * dthetai
       cthetai(j) = cos( thetai(j) * radianToDegree )
       if (j > 1) domo2i(j) = pi * ( cthetai(j-1) - cthetai(j) )
    enddo


    ! Establish photon event generator data (based on course angle mesh)
    photonEG = newPhotonEventGenerator( theta, 19, dataIO )
    if ( .not.photonEG%properlyConstructed() ) then
       write(message, 2400)
       call dataIO(0, 2, message)
       write(message, 2500)
       call dataIO(0, 2, message)
       errorFlag = errorFlag + 1
    end if



    ! Verify file exists and can be opened with the desired file unit
    gammaFile = desiredGammaFile
100 continue

    ! File is NOT open; verify file exists.
    inquire(file = gammaFile, exist = fileExists )
    if ( fileExists ) then
       ! The file exists and an available unit exists (file can be opened)
    else
       ! File doesn't exist; look for the default (if not being used)
       write(message, 2000) gammaFile   ! Warn client that file doesn't exist; try default if not used
       call dataIO(0, 2, message)
       if ( gammaFile /= defaultGammaFile ) then
          gammaFile = defaultGammaFile
          write(message, 2100) gammaFile
          call dataIO(0, 3, message)
          go to 100
       else
          ! Could not find the file (client defined or default); return with error code
          errorFlag = sDCMGammaError
       end if
    end if
    if ( errorFlag == sDCMGammaError ) return

    ! Create photon cross section data arrays
    if (.not.allocated(gppipn)) allocate(gppipn(50, 20))
    if (.not.allocated(gppi0p)) allocate(gppi0p(50, 20))
    if (.not.allocated(elg))    allocate(   elg( 4, 50))

    ! Allocate temporary data now (deallocate later)
    if (.not.allocated(xsectd)) allocate(xsectd(4, 50, 0:18))
    if (.not.allocated(ecm)) allocate(ecm(4, 50))

    ! Read differential cross section data file:
    ! NOTE: ONLY THE FIRST 2 CHANNELS (OF 50) ARE UTILIZED HERE (see "jch" values below)
    open( newunit = gammaUnit, file = gammaFile, status = "old", action = "read" )
    do jch = 1, 2
       readError(1) = 0
       read (gammaUnit, 10, iostat = readError(1) )
       do inw = 1,50
          readError(2:3) = 0

          ! Reading data
          read (gammaUnit, 20, iostat = readError(2) ) &
               & ecm(jch,inw), (xsectd(jch,inw,inth), inth=0,4)
          read (gammaUnit, 30, iostat = readError(3) ) &
               & (xsectd(jch,inw,inth), inth=5,18)


          ! Verify no errors or End Of File (EOF) reached during read
          do i = 1, 3
             if ( readError(i) /= 0 ) then
                if ( readError(i) > 0 ) then
                   ! Unknown error during read (not all variables initialized)
                   errorFlag = sDCMGammaError
                   write(message, 2300) trim(gammaFile)
                   call dataIO(0, 2, message)
                else if ( readError(i) < 0 ) then
                   ! EOF encountered (not all variables initialized)
                   errorFlag = sDCMGammaError
                   write(message, 2350) trim(gammaFile)
                   call dataIO(0, 2, message)
                end if
                
                write(message, 2500)
                call dataIO(0, 2, message)
                go to 200

             end if
          end do

          ! Manipulating data and storing result
          temp = two*emnucg
          elg(jch,inw) = (ecm(jch,inw)**2 - emnucg**2)/(temp)


          ! Store channel data for the reaction:
          do ith = 1,19
             if (jch == 1) then
                !  gamma + p --> pi+ + n
                gppipn(inw,ith) = xsectd(jch,inw,ith-1)
             elseif (jch == 2) then
                !  gamma + p --> pi0 + p
                gppi0p(inw,ith) = xsectd(jch,inw,ith-1)
             endif
          end do

       end do
    end do
200 continue

    close (gammaUnit)
    return

! ======================================================================
! (for reading data file)
10  format (1x)
20  format (7x,f5.2,8x,5e10.3)
30  format (7e10.3)
! (for messages to client)
2000 format("File '", A, "' could not be found for reading in photon ", &
          & "cross section data from.")
2100 format("Attempting to use file '", A, "' for photon cross section data.")
2300 format("An unkown error was encountered when reading the photon ", &
          & "cross section data file '", A, "'.")
2350 format("Unexpected end of file reached when reading the photon ", &
          & "cross section data file '", A, "'.")
2400 format("The data for the photon event generator failed to setup.")
2500 format("   Be wary of simulation results.")
! ======================================================================
  end function initializeGammaData
