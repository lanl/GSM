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
        & insist,   &
        & validate, &
        & require

! Define various contract levels
#define PRODUCTION_LEVEL    0
#define PRECONDITION_LEVEL  1
#define INTRASCOPE_LEVEL    2
#define POSTCONDITION_LEVEL 3
#define DEBUG_LEVEL         4

! Define the Contracts level for the build
#ifndef Contracts_Level
#define Contracts_Level DEBUG_LEVEL
#endif

   logical, parameter, private:: &
       & showLine = (Contracts_Level > PRODUCTION_LEVEL)

! Define the various macros
#include "NewLine.f90"
#include "Insist.f90"
#include "Validate.f90"
!#include "Require.f90"
!#include "Check.f90"
!#include "Ensure.f90"
!#include "Remember.f90"
!#include "NotImplemented.f90"
!#include "NotConfigured.f90"
!#include "NotReachable.f90"

   ! >>> DATA OBJECTS
contains

   ! >>> CONSTRUCTORS

   ! >>> INTROSPECTION

   ! >>> SETTERS

   ! >>> PROCEDURES
#include "insist.f90"
#include "validate.f90"
#include "require.f90"

end module Contracts

