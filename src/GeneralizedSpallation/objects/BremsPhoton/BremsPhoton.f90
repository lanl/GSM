
! ==============================================================================
!
! Specifies data for a brems. photon particle
!
!
! Written by CMJ, XCP-3, 11/2019
!
! ==============================================================================
module BremsPhotonMod

  use, intrinsic:: iso_fortran_env, only: real64

  implicit none
  private

  ! General constructors
  public :: newBremsPhoton
  interface newBremsPhoton
     module procedure constructorMain
  end interface newBremsPhoton

  ! Load brems. photon data type
  include "bremsPhoton.f90"

contains

  ! Load constructors
  include "constructors.f90"

  ! Load getters
  include "getters.f90"

  ! Load setters
  include "setters.f90"

  ! Object introspection
  include "introspection.f90"

end module BremsPhotonMod

