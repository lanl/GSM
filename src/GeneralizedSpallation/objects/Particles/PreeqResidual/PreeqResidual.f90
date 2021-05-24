
! ==============================================================================
!
!> \file
!> \brief Module file for the PreeqResidual object
!>
!> This module defines the particle type used to obtain statistics on the
!> preequilibrium residual before and after de-excitation.
!
!
! Written by CMJ, XCP-3, 11/2019
!
! ==============================================================================
module PreeqResiduals

  use, intrinsic:: iso_fortran_env, only: int32, int64, real64

  implicit none
  private

  type, public:: PreeqResidual
     private

     !> \brief Number of nucleons in the nucleus (A number)
     real(real64), public:: numBaryons = 0.0_real64

     !> \brief Number of protons in the nucleus (Z number)
     real(real64), public:: numProtons = 0.0_real64

     !> \brief Thermalized energy of the nucleus, in [MeV]
     real(real64), public:: energy = 0.0_real64

     !> \brief Total linear momentum of the nucleus [GeV/c]
     real(real64), public:: linearMomentum = 0.0_real64

     !> \brief Angular momentum quantum number
     integer(int32), public:: angMom = 0_int32

   contains
     private
  end type PreeqResidual


contains

  ! Load constructors
  ! include "constructors.f90"

  ! Load getters
  ! include "getters.f90"

  ! Load setters
  ! include "setters.f90"

  ! Object introspection
  ! include "introspection.f90"

end module PreeqResiduals

