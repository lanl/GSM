! =============================================================================
!>
!> \file
!> \brief  Contains the Testing parameters
!> \author CMJ
!>
! =============================================================================

!> \brief Provides indicators to state a test passed, failed, or was not executed
enum, bind(C) :: TestStatus
    enumerator :: FAILED
    enumerator :: PASSED
    enumerator :: UNTESTED
end enum

