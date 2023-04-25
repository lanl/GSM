
! ==============================================================================
!
! Module file for output data
!
!
! Written by CMJ, XCP-3, 11/2019
!
! ==============================================================================
module OutputDataMod

  use, intrinsic:: iso_fortran_env, only: int32, int64, real64

  implicit none
  private

  ! General constructors
  public :: newOutputData
  interface newOutputData
     module procedure constructorMain
  end interface newOutputData

  ! Load output data data type
  include "outputData.f90"

contains

  ! Load constructors
  include "constructors.f90"

  ! Load getters
  ! include "getters.f90"

  ! Load setters
  ! include "setters.f90"

  ! Object introspection
  ! include "introspection.f90"

end module OutputDataMod

