! =============================================================================
!
!> \file
!> \brief  Contains the main Logger testing procedure
!> \author CMJ (XCP-3; LANL)
!
!> Contains the tests for the Loggers module.
!
! =============================================================================
!
!> \fn     tstLogger
!> \brief  Makes calls to all the Logger tests for verification of usage
!> \author CMJ (XCP-3; LANL)
!
!> Implements the various test procedures for the Loggers module
!
! =============================================================================
program tstLogger
    use tstLogger
    implicit none

    ! Call test procedures
    call tstLoggerLevel()
    !call tstLoggerStream()
    !call tstLoggerFilters()
    !call tstLoggerPrinters()

end program tstLogger

