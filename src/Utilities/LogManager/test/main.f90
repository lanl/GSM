! =============================================================================
!
!> \file
!> \brief  Contains the main LogManager testing program
!> \author CMJ (XCP-3; LANL)
!
!> Program that calls the LogManager testing procedures
!
! =============================================================================
!
!> \fn     tstLogManager
!> \brief  Makes calls to all the LogManager tests for verification of usage
!> \author CMJ (XCP-3; LANL)
!
!> Implements the various test procedures for the LogManagers module
!
! =============================================================================
program tstLogManager
    use tstLogManager
    implicit none

    ! Call test procedures
    call tstSizing()
    call tstGetters()
    call tstMessaging()

end program tstLogManager
