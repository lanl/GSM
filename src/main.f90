
  program gsm1

! ====================================================================
!
! This is the main program for the Generalized Spallation Model (GSM).
!
! Here, we simply call a module which acts as a "driver" for GSM.
!
!
! ====================================================================

    use standaloneGSM, only: gsmDriver
    implicit none

! ====================================================================

    call gsmDriver()

! ====================================================================
  end program gsm1
