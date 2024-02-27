
! ====================================================================
!
! This module contains the definitions of child classes for the
! Modified Dubna Cascade Model
!
! ====================================================================

module mDCMChildTypes

  use, intrinsic:: iso_fortran_env, only: int32, real64
  use numbers
  implicit none
  private


  ! Include various derived types utilized:
#include "types.f90"

contains

! ====================================================================

  ! Constructors:

  ! Querying the objects:
#include "access.f90"


! ====================================================================

end module mDCMChildTypes
