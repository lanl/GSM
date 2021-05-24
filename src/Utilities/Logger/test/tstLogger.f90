! =============================================================================
!
!> \file
!> \brief  Contains the tstLogger module
!> \author CMJ (XCP-3; LANL)
!
!> Contains the tests for the Loggers module.
!
! =============================================================================
module tstLogger

    ! Import the Logger object and its constructor
    use Loggers, only: Logger, newLogger

    implicit none
    private

    public:: tstLoggerLevel
    !public: tstLoggerStream
    !public: tstLoggerFilters
    !public: tstLoggerPrinters

 contains

     include "tstLoggerLevel.f90"
     !include "tstLoggerStream.f90"
     !include "tstLoggerFilters.f90"
     !include "tstLoggerPrinters.f90"

end module tstLogger

