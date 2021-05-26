! =============================================================================
!>
!> \file
!> \brief  Contains the Exceptions module
!> \author CMJ
!>
!> The Exceptions module contains a set of objects to help facilitate clear
!> error handling that is abstracted from the consumers to provide more control
!> over the application when an unexpected event occurs.
!>
! =============================================================================
module FExceptions

  use Types, only: gsmInt8
  implicit none
  private

  !> \brief Quick reference for consumers to "throw" a FException
#define FThrow(fExcept, errText) call fExcept%throw( __FILE__, __LINE__, errText)

! Import module parameters and data types
#include "parameters.F90"
#include "fException.F90"
  

  !> \brief Default method to "throw" exceptions
  private :: throwExceptionBase
  
  !> \brief The procedure that is called to "throw" exceptions
  procedure(throwExceptionBase), private, pointer :: throwException => throwExceptionBase
  
  !> \brief Sets the default method that is called to "throw" exceptions
  public :: setExceptionThrower
  
contains

! Import the FException procedures
#include "baseProcedures.F90"

! Import derived class procedures
#include "derivedProcedures.F90"

end module FExceptions
