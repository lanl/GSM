
  subroutine readDecayData(decayFileIN, printStatusIN)

! ====================================================================
!
! This subroutine reads in the decay table and establishes decay
! data and mappings for the quark-based indent code.
!
! ====================================================================

    implicit none
    character(*), intent(in   ), optional :: decayFileIN
    logical,      intent(in   ), optional :: printStatusIN

    character(:), allocatable :: decayFile
    logical :: printStatus = .false.

    ! Values read from the data file
    integer(int32) :: ires   ! CUrrent particle ID
    real(real64)   :: br
    integer(int32), dimension(6) :: imode

    ! Interim variables
    logical :: &
        &      exitLoop = .false., &
        &      loopExceeded = .false.
    integer(int32) :: &
        & decayUnit, & ! File unit for the decay table
        & index, &     ! Particle index in lookup tables
        & iold, &      ! Prior particle ID
        & loop, &      ! Lookup index (also the loop number)
        & rc           ! File I/O status

! ====================================================================

    ! Parameters and unimportant temp values
    integer(int32) :: i, ifl1, ifl2, ifl3, itype, jspin, k
    character(len=8) :: lread(10), lmode(6), lres
    integer(int32), parameter :: loopLimit = 600_int32
    character(len=*), parameter :: &
         & iquit  = " end", &
         & iblank = "     "

! ====================================================================

    ! LOOK MUST BE DIMENSIONED TO THE MAXIMUM VALUE OF INDEX
    ! We should set these more dynamically; perhaps bigger than 600 and make
    ! these allocatable? Further, we can restrict their use to ensure
    ! indexes are in range of what's allocated
    integer(int32) :: look, mode
    real(real64)   :: cbr
    common/dkytab/look(400),cbr(600),mode(5,600)

! ====================================================================

    ! Parameter validation
    decayFile = defaultDecayFile
    if(present(decayFileIN)) then
       decayFile = decayFileIN
    end if
    if(present(printStatusIN)) then
       printStatus = printStatusIN
    end if
    if(printStatus) write(output_unit,10)

    ! Establish internal decay channel data
    defaultChannels%modelDecayOpt = .true.
    call defaultChannels%allocateData()

    ! Setup preliminary values
    loop=0
    iold=0
    exitLoop = .false.
    loopExceeded = .false.

    ! Check if file is already opened somewhere
    ! If so, rewind it / start at the front
    inquire(file = decayFile, &
        & number = decayUnit)
    if (decayUnit /= -1_int32) rewind decayUnit

    ! Start main reading loop
    primary: do while (.not. exitLoop)
       ! Increment loop and check if limit exceeded
       loop = loop + 1_int32
       exceeded: if (loop > loopLimit) then
          exitLoop = .true.
          loopExceeded = .true.
       else
! 220       continue
          ! Reset temporary variables
          imode(:) = 0_int32
          lmode(:) = iblank

          ! Open file if it's not already opened
          inquire(file = decayFile, &
             & number = decayUnit)
          if (decayUnit == -1_int32) then
             open(newunit = decayUnit, &
                & file = decayFile, &
                & status = "old", action = "read", iostat = rc)
             call insist(rc == 0, &
                & "Failed to open file: " // decayFile, &
                & __FILE__, &
                & __LINE__)
          end if

          ! Read file
          read(decayUnit, *, iostat = rc) &
             & ires, itype, br, (imode(i), i = 1, 5)
          call insist(rc <= 0, &
             & "A failure occurred when reading from decay data file: " // decayFile, &
             & __FILE__, &
             & __LINE__)

          processdata: if (ires == 0 .or. rc < 0) then
             ! Reached end of file; exit loop
             exitLoop = .true.
          else
             ! "ires" here is the particle ID associated to the decay channel
             ! IF(NOPI0.AND.IRES.EQ.110) GO TO 220
             ! IF(NOETA.AND.IRES.EQ.220) GO TO 220

             ! Assign a lookup number for a unique particle
             if(ires /= iold) then
                ! Get particle properties and associated index
                ! Assign into the lookup array [look]
                ! NOTE: we can map the "cbr" or "mode" values using the particle
                ! ID by obtaining the lookup value (to get the associated loop #),
                ! then can obtain the cbr or mode lookup using the loop number.
                call flavor(ires, ifl1, ifl2, ifl3, jspin, index)
                defaultChannels%lookDat(index) = loop
             end if
             iold = ires

             ! Assign specific data values associated to the particle
             defaultChannels%cbrDat(loop) = br
             do i = 1, 5
                defaultChannels%modeDat(i, loop) = imode(i)
                if(imode(i) /= 0) call label(lmode(i), imode(i))
             end do
             call label(lres, ires)

             ! Print status, if desired
             if (printStatus) then
                write(output_unit,20) &
                    & lres, (lmode(k), k = 1, 5), &
                    & br, ires, (imode(k), k = 1, 5)
             end if
          end if processdata
       end if exceeded
    end do primary
    defaultChannels%numEntriesDat = loop

    ! TODO: Remove common block
    look(:) = defaultChannels%lookDat(:)
    cbr(:) = defaultChannels%cbrDat(:)
    mode(:, :) = defaultChannels%modeDat(:, :)

    ! Read and print notes from the decay table file
    if (.not.loopExceeded .and. printStatus) then
        exitLoop = .false.
        readNotes: do while (.not.exitLoop)
           ! Read
           read(decayUnit, 1002, iostat=rc) lread

           ! If EOF, error, or last line, then exit
           ! Otherwise, print the line back to the user
           if (rc /= 0 .or. lread(1) == iquit) then
               exitLoop = .true.
           else
              write(output_unit, 1003) lread
           end if
       end do readNotes
    end if
    close(decayUnit)


    ! If loop limited was exceeded, then print error and close file
    if (loopExceeded) then
        write(output_unit, 30) __FILE__, &
            & "array size exceeded at ", &
            & loopLimit
    end if

    return
! ====================================================================
10  format(1h1,30(1h*)/2h *,28x,1h*/ &
         & 2h *,5x,17hcolli decay table,5x,1h*/ &
         & 2h *,28x,1h*/1x,30(1h*)// &
         & 6x,4hpart,18x,10hdecay mode,19x,6hcum br,15x,5hident,17x, &
         & 11hdecay ident/)
20  format(6x,a5,6x,5(a5,2x),3x,f8.5,15x,i5,4x,5(i5,2x))
30  format(//44h, "***** Error in ", A, ": ", A, i0, "*****")
1002   format(10a8)
1003   format(1x,10a8)
! ====================================================================
  end subroutine readDecayData
