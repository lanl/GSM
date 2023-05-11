!==============================================================================
!> Harness ExceptionLite library
!>
!> Library to centralize error creation, management, and handling. This should
!> contain one base error type, several specializations, and special methods to
!> create, retrieve, and clear lite exceptions.
!>
!> Consumers can use these via any of the following methods:
!>    - throwRunTimeException - consumer "throws" (creates) an exception
!>    - throwLogicException - consumer "throws" (creates) an exception
!>    - throwException - consumer "throws" (creates) an exception
!>    - retrieveException - allows consumer to "catch" (retrieve) an exception
!>    - handleException - allows consumer to flag an exception as handled
!>    - throwUnhandled - halts execution for an unhandled exception
!>
!> Note that the methods used to create and clear lite exceptions should manage
!> the private exception vector.
!
!==============================================================================
module hrnExceptionsLite

   use, intrinsic :: iso_fortran_env, only: int64, error_unit
   implicit none
   private

   public :: throwUnhandled
   public :: throwRunTimeException
   public :: throwLogicException
   public :: throwException
   public :: retrieveException
   public :: handleException

   private :: retrieveID

   private :: throwExceptionInternal

   ! >>> MODULE MEMBERS
   include "parameters.f90"

   ! >>> DATA OBJECTS
   include "exceptionLite.f90"

contains

   ! >>> PRIVATE MODULE / CLASS PROCEDURES
   include "procedures.f90"

   ! >>> PUBLIC PROCEDURES
   include "public_procedures.f90"

end module hrnExceptionsLite

