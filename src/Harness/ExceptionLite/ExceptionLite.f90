!==============================================================================
!> Harness ExceptionLite library
!>
!> Library to centralize error creation, management, and handling. This should
!> contain one base error type, several specializations, and special methods to
!> create, retrieve, and clear lite exceptions.
!>
!> Consumers can use these via any of the following methods:
!>    - throwRunTimeException - consumer "throws" (creates) a runtime exception
!>    - throwLogicException - consumer "throws" (creates) a logic exception
!>    - throwException - consumer "throws" (creates) an generic exception
!>    - retrieveException - allows consumer to "catch" (retrieve) an exception
!>    - handleException - allows consumer to flag an exception as handled
!>    - throwUnhandled - halts execution for an unhandled exception
!>
!> Note that the methods used to create and clear lite exceptions should manage
!> the private exception vector.
!
!==============================================================================
module hrnExceptionsLite

   use, intrinsic :: iso_fortran_env, only: int16, error_unit
   implicit none
   private

   ! Methods to "throw" an exception (wrappers)
   public :: throwRunTimeException
   public :: throwLogicException
   public :: throwException
   
   ! Methods to retrieve and handle exceptions
   public :: throwUnhandled
   public :: retrieveException
   public :: handleException

   ! Used to set the exception ID. These should be unique for every exception
   ! created.
   public :: getNewID

   ! Internal method to throw an exception by enum type
   private :: throwExceptionInternal

   ! >>> MODULE MEMBERS
   include "parameters.f90"

   ! >>> DATA OBJECTS
   include "objects.f90"

contains

   ! >>> PRIVATE MODULE / CLASS PROCEDURES
   include "exception_procedures.f90"

   ! >>> PUBLIC PROCEDURES
   include "public_procedures.f90"

end module hrnExceptionsLite

