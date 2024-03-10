
! ====================================================================
!
! Module file that contains the parameterized data utilized by the
! modified  Dubna Cascade Model (DCM)
!
! Written by CMJ, XCP-3, 2/2024
!
! ====================================================================

module modifiedDCMDataClass

  ! Import general utilities used globally
  use, intrinsic:: iso_fortran_env
  use numbers
  use Contracts

  implicit none
  private

  ! For printing:
  integer(int32), private, parameter :: defaultMDCMDataVerbose = 4_int32
  integer(int32), public :: mDCMDataVerbose = defaultMDCMDataVerbose

  ! For construction errors:
  integer(int32), public, parameter :: objectNotConstructed = -1_int32
  integer(int32), public, parameter :: properlyConstructed  =  0_int32
  integer(int32), public, parameter :: invalidTargetFlag    =  1_int32
  integer(int32), public, parameter :: targetInitFailed     = 10_int32


  ! Default options:


  ! For divide by zero checks
  real(real64),  private, parameter :: div0Lim = 1.0d-15


  ! Interface to establish a new data class
  public :: newModifiedDCMData
  interface newModifiedDCMData
     module procedure :: mDCMDataMainConstructor
  end interface


  abstract interface
     ! Interface for I/O handling
     subroutine IOHANDLER(verbosity, type, text)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in) :: verbosity
       integer(int32),   intent(in) :: type
       character(len=*), intent(in) :: text
     end subroutine IOHANDLER
  end interface

  ! Include target data type
  include "mDCMDataClassTypes.f90"

  type, public :: ModifiedDCMData
     private
     ! Flag indicating if the data object was constructed
     logical, private :: objectConstructed = .FALSE.
     integer(int32), private :: constructed = objectNotConstructed

     ! I/O Handling (all messages filter through this type)
     type(mDCMDataIO),      private :: io

     ! Contains all options for the object
     type(mDCMDataOptions), private :: options

   contains
     private
     ! For querying object:
     procedure, public :: properlyConstructedObject
     procedure, public :: properlyConstructed => properlyConstructedObject
     procedure, public :: constructionState
     procedure, public :: queryOptions

     ! Private Member(s)
     procedure, private :: validateOptions

  end type ModifiedDCMData

  ! Non-class specific routines
  private :: mDCMDataMainConstructor   ! Gives client mDCM data class
  private :: printmDCMData         ! Handles all printing from mDCM data class

contains

! ====================================================================

  ! Include constructor and its associated procedures:
#include "mDCMDataMainConstructor.f90"
#include "validateOptions.f90"

  ! Include public procedures for clients to access data
#include "mDCMDataAccess.f90"      ! Gives a client access to the members of'target'
#include "constructionState.f90"   ! Tells client if construction was successful

  ! Include module-specific routines (used by data class, but class use not necessary

  ! Default I/O Handler for the data class
#include "printMDCMData.f90"

! ====================================================================
end module ModifiedDCMDataClass
