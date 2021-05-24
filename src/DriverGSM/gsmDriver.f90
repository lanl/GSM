
  subroutine gsmDriver()

! ====================================================================
!
! This is the driver procedure for end-users to call to simulate GSM
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32
    use OpenMPAbstraction, only: &
         & multiThreadingAvailable, &
         & numAvailableThreads, &
         & setNumThreads

    implicit none

    integer(int32) :: numCmdErr = 0_int32
    character(LEN=128) :: gsmInput = defaultGSMInputName

    integer(int32),     dimension(8) :: timeValues
    character(LEN= 10), dimension(3) :: strTimeVals

! ====================================================================

    ! Print GSM version info to user:
    call date_and_time( strTimeVals(1), strTimeVals(2), strTimeVals(3), &
         & timeValues )
    write(*, 1000)
    write(*, 1100) "started", timeValues(5), timeValues(6), timeValues(7), timeValues(8), &
         & timeValues(1), timeValues(2), timeValues(3)


    ! Obtain, set, and validate all commands:
    call parseCommands( gsmInput, numCmdErr )

    ! HACK: Use all available OpenMP resources if GSM is built with it
    if (multiThreadingAvailable()) then
       call setNumThreads(numAvailableThreads())
    end if


    ! Read Input and perform simulation if no fatal command errors found:
    if ( numCmdErr == 0 ) then
       write(*, 1999)
       call readInput( trim(gsmInput) )
    end if


    ! Print end of simulation:
    call date_and_time( strTimeVals(1), strTimeVals(2), strTimeVals(3), &
         & timeValues )
    write(*, 1100) "ended", timeValues(5), timeValues(6), timeValues(7), timeValues(8), &
         & timeValues(1), timeValues(2), timeValues(3)

    return
! ====================================================================
1000 format(                               2/, &
          & " ==========================", 1/, &
          & " This is GSM, Version 1.0.0", 1/, &
          & " =========================="      )
1100 format(1/, &
          & "The simulation ", A, " at ", &
          &      i2, ":", i2, ":", i2, ".", i3, ",",    /, &   ! Time
          & "    on ", &
          &      i4, "-", i2, "-", i2, ".")   ! Date
1999 format("")
! ====================================================================
  end subroutine gsmDriver
