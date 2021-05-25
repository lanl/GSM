! =============================================================================
!>
!> \file
!> \brief  Contains the LiteStrings module
!> \author CMJ
!>
!> The LiteStrings module contains functionality to build messages, convert
!> basic numbers to a string format, and other similar functionality as
!> required for GSM usage.
!>
! =============================================================================
module LiteStrings

  use Types, only: gsmInt8, gsmInt16, gsmInt32, gsmInt64, &
       & gsmFloat, gsmDouble, gsmLongDouble

  implicit none
  private

  ! >>> PARAMETERIZED DEFAULTS
#include "parameters.F90"

  ! >>> INTERFACES
#include "conversionInterface.F90"

  ! >>> OBJECT CONSTRUCTION
!#include "constructorInterface.F90"

  ! >>> DATA OBJECTS
!#include "liteString.F90"

contains

#include "conversionFunctions.F90"


end module LiteStrings

