! =============================================================================
!
!> \file
!> \brief  Contains the Contracts module
!> \author CMJ (XCP-3; LANL)
!
!> The Contracts module contains the various functions and macros for the
!> design-by-contract specifications. This module relies heavily on the Fortran
!> preprocessor.
!
! =============================================================================
module Contracts

  implicit none
  private

  ! Contract functions
  public:: &
       & require, &
       & check, &
       & ensure, &
       & insist, &
       & not_implemented


#include "contract_macros.fpp"

    !> \brief Indicates if the line details should be shown
    logical, parameter, private :: &
        & showLine = (CONTRACTS_LEVEL > PRODUCTION_LEVEL)

    !> \brief Specifies a new line
    character(*), parameter, private:: &
        & newLine = achar(13) // achar(10)

    ! >>> DATA OBJECTS
contains

     ! >>> HELPERS
#include "append_location.f90"
#include "stop_execution.f90"

     ! >>> PROCEDURES
#include "require.f90"
#include "check.f90"
#include "ensure.f90"
#include "insist.f90"
#include "not_implemented.f90"

end module Contracts
