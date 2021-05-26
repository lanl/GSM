! =============================================================================
!>
!> \file
!> \brief  Contains the Basics module
!> \author CMJ
!>
!> The Basics module contains the absolute basic utilities that are used
!> in GSM. These utilities include things like abstracted data types, basic
!> string manipulation, and so forth. The module facilitates ease of use for
!> it's consumers. The module wraps several of these basic sub-modules
!> together to provide it's API / usage.
!>
! =============================================================================

#include "Types.F90"
#include "FExceptions/FExceptions.F90"
! #include "LiteStrings/LiteStrings.F90"


! =============================================================================
!>
!> \module
!> \brief  Provides basic abstraction and functionality for simple tasks
!> \author CMJ
!>
!> The Basics module contains the absolute basic utilities that are used
!> in GSM. These utilities include things like abstracted data types, basic
!> string manipulation, and so forth. The module facilitates ease of use for
!> it's consumers. The module wraps several of these basic sub-modules
!> together to provide it's API / usage.
!>
! =============================================================================

module Basics

    use Types
    use FExceptions
  
  implicit none
  public


end module Basics

