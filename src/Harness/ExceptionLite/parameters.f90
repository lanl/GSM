!==============================================================================
!> ExceptionLite module parameters
!>
!> Specifies all private and parameterized ExceptionLite module data
!>
!==============================================================================

!> Define whether or not exceptions will include details about the file
!> and line of the exception.
#ifndef ExLocationDetails
#define ExLocationDetails .true.
#endif

   !> Indicates what the next exception ID will be.
   integer(int64), private :: next_id = 1

   !> Defines a general exception type
   !>
   !> General exceptions should not include those of any other defined type.
   public :: GENERAL_EXCEPTION

   !> Defines a logic exception
   !>
   !> Define exceptions related to program logic, such as requirements on
   !> arguments, intermediate results, output results, lengths, domains, etc.
   public :: LOGIC_EXCEPTION

   !> Defines a runtime exception
   !>
   !> Runtime exceptions are those that occur at runtime. One such use case
   !> includes mathematicl under- or over-flow errors, system errors, file
   !> errors, and anything else that may occur at runtime not related to program
   !> logic and internal requirements.
   public :: RUNTIME_EXCEPTION

   !> Define several enums for exception handling (internally used)
   enum, BIND(C)
      enumerator :: GENERAL_EXCEPTION
      enumerator :: LOGIC_EXCEPTION
      enumerator :: RUNTIME_EXCEPTION
   end enum

