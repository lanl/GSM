!==============================================================================
!> Harness constants library
!>
!> Library to contain the various contants utilized by consumers. This should
!> contain several shortcuts to simplify usage, such as "0", "1", "1/2", pi,
!> etc.
!>
!> All constants should be publicly available. Arrayed constants are accessed
!> via functions that validate parameters.
!
!==============================================================================
module hrnConstants

  use, intrinsic :: iso_fortran_env, only: int32, real64, error_unit
  implicit none
  private

#include "numbers.f90"

#include "physics.f90"

contains

#include "physics_interface.f90"


end module hrnConstants

