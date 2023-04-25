
! ==============================================================================
!
! Module file for non-particle event specific data
!
!
! Written by CMJ, XCP-3, 11/2019
!
! ==============================================================================
module EventDataMod

  use, intrinsic:: iso_fortran_env, only: int32, int64, real64
  use EvaporationClass, only: EvapCompound => residual, FissionFragment
  use PreeqResiduals, only: PreeqResidual

  implicit none
  private

  ! General constructors
  public :: newEventData
  interface newEventData
     module procedure constructorMain
  end interface newEventData

  ! Load EventData data type
  include "eventData.f90"

contains

  ! Load constructors
  include "constructors.f90"

  ! Load getters
  ! include "getters.f90"

  ! Load setters
  ! include "setters.f90"

  ! Object introspection
  ! include "introspection.f90"

end module EventDataMod

