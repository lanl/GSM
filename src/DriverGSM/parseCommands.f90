
  subroutine parseCommands( gsmInput, numCmdErr )

! ====================================================================
!
! This is the main procedure by which GSM can setup and simulated:
!
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsmClass, only: gsmVerbose

    implicit none
    character(LEN=*), intent(inout) :: gsmInput
    integer(int32),   intent(  out) :: numCmdErr

    integer(int32) :: convertError, newVerbose
    real(real64)   :: newRealVerbose

! ====================================================================

    ! Accepted commands:
    integer(int32),  parameter :: cmdArgLen = 3_int32
    character(LEN=cmdArgLen), parameter :: &
         & inpCmd  = "-i=", &
         & verbCmd = "-v="
    character(LEN=*), parameter :: helpCmd = "-help"
    integer(int32),   parameter :: helpLen = len( helpCmd )

! ====================================================================

    ! For commands:
    integer(int32)     :: cmdNum, cmdLen, cmdStatus
    character(LEN=512) :: totalCmd = ""   ! Contains the entire command line
    character(LEN=128) :: cmdStr   = ""   ! Contains the a single command (flag and body)
    character(LEN=cmdArgLen) :: cmdFlag = ""   ! The flag of the command
    character(LEN=125) :: cmdBody = ""         ! The body of the command

! ====================================================================

    ! Print out commands to user:
    call get_command( totalCmd )
    write(*,1000) trim(totalCmd)


    ! Obtain command line arguments at runtime:
    parseCommandLoop: do cmdNum = 1, command_argument_count()
       call get_command_argument( cmdNum, cmdStr, cmdLen, cmdStatus )
       if ( cmdStatus /= 0 ) then
          write(*,1100) "Command argument ", cmdNum, " failed."
          cycle parseCommandLoop
       end if

       ! Provide user with some 'help' if desired:
       if ( cmdStr( :helpLen) == helpCmd ) then
          call printHelp( inpCmd, verbCmd, helpCmd )
          cycle parseCommandLoop
       end if

       if ( cmdLen < (cmdArgLen+1) ) then
          write(*, 2000) cmdNum, trim(cmdStr)
          write(*, 2110)
          cycle parseCommandLoop
       end if

       ! Obtain command flags and bodies:
       cmdStr  = trim(adjustl(  cmdStr ))   ! Remove all leading and trailing spaces
       cmdFlag = trim(adjustl(  cmdStr(1:cmdArgLen)   ))   ! Obtain command flag
       cmdBody = trim(adjustl(  cmdStr(cmdArgLen+1:)  ))   ! Obtain command body

       ! Determine the command type:
       select case( cmdFlag )
          case( inpCmd )
             gsmInput = cmdBody
!             write(*,1200) trim(gsmInput)

          case( verbCmd )
             read(cmdBody, *, iostat=convertError) newRealVerbose
             if ( convertError == 0 ) then
                newVerbose = nint( newRealVerbose, int32 )
                gsmVerbose = newVerbose
!                write(*, 1300) gsmVerbose
             else
                write(*, 2200) cmdNum, trim(cmdBody)
             end if

          case default
             write(*, 2100) cmdNum, trim(cmdStr)
             write(*, 2110)

       end select
    end do parseCommandLoop

    ! Validate all commands:
    call validateCommands( gsmInput, numCmdErr )

    return
! ====================================================================
1000 format("Invoking GSM with the following command line:", 1/, &
          & "   ", A, 1/ )
1100 format(A, i2, A)
! 1200 format("   '", A, "' will be utilized as the input file for the simulation.")
! 1300 format("   The simulation will utilize a verbosity of ", i5, ".")
2000 format("   Command ", i2, ": The command does not have a body (", A, ").")
2100 format("   Command ", i2, ": Command not recognized: ", A)
2110 format("      The command will be ignored.")
2200 format("   Command ", i2, ": Unable to set simulation verbosity given the command's ", &
          & "body: ", A)
! ====================================================================
  end subroutine parseCommands
