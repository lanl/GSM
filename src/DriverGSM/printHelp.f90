
  subroutine printHelp( inpCmd, verbCmd, helpCmd )

! ====================================================================
!
! Prints out all command line to options to the user
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    character(LEN=*), intent(in   ) :: inpCmd
    character(LEN=*), intent(in   ) :: verbCmd
    character(LEN=*), intent(in   ) :: helpCmd

! ====================================================================

    write(*, 1000)
    write(*, 1100)   ! Invocation:
    write(*, 1200) & ! Command arguments
         & trim(inpCmd), defaultGSMInputName, &
         & trim(verbCmd), &
         & trim(helpCmd)

    return
! ====================================================================
1000 format(                                                       1/, &
          & " ***************   GSM Usage   ***************",      1/, &
          & " *********************************************"           )
1100 format(" GSM can be invoked by typing './{execName}', where", 1/, &
          & "    {execName} is the name of the GSM executible,",   1/, &
          & "    such as 'xgsm1', for example.",                   1/  )
1200 format(" Command line options for GSM execution include:",    1/, &
          & " ===============================================",    1/, &
          & " ", A8, ": Flags the input file name to be read ",        &
          & "(default is '", A, "')",                              1/, &
          & " ", A8, ": Sets the verbosity of the simulation",     1/, &
          & " ", A8, ": Provides information on GSM usage and ",       &
          & "command arguments",                                   1/, &
          &                                                        1/, &
          & "    NOTE: File names, paths, or both, must be",       1/, &
          & "          enclosed by quotation markers to be",       1/, &
          & "          properly accepted.",                        1/, &
          & " ===============================================",    1/, &
          &                                                        1/  )
    ! ====================================================================
  end subroutine printHelp
