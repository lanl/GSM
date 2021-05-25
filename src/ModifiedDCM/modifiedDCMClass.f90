
! ====================================================================
!
! This module contains the modified Dubna Cascade Model objects
!
! ====================================================================

module modifiedDCMClass

  use, intrinsic:: iso_fortran_env, only: int32, real64
  implicit none
  private

  ! For verbosity:
  integer(int32), private, parameter :: defaultMDCMVerbose = 4_int32
  integer(int32), public :: mDCMVerbose = defaultMDCMVerbose

  ! Default options:
  real(real64), private, parameter :: defaultRM    = 0.2_real64
  real(real64), private, parameter :: defaultDelta = 1.3_real64

  ! For divide by zero limits:
  real(real64), private, parameter :: div0Lim = 1.0d-15


  ! Procedure interfaces:
  abstract interface
     subroutine IOHANDLER(verbosity, type, text)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in) :: verbosity
       integer(int32),   intent(in) :: type
       character(len=*), intent(in) :: text
     end subroutine IOHANDLER
     function RANDOM() result(rndmNum) BIND(C)
       use, intrinsic:: iso_C_binding, only: c_double
       implicit none
       real(c_double) :: rndmNum
     end function RANDOM
  end interface



  ! Constructor interfaces:
  ! (class constructor)
  public :: newModifiedDCM
  interface newModifiedDCM
     module procedure :: new_ModifiedDCM
  end interface
  ! (results constructor)
  public :: newMDCMResults
  interface newMDCMResults
     module procedure :: new_Results
  end interface

  ! Include various derived types utilized:
  include "mDCMClassTypes.f90"

  ! FOR CLIENTS TO USE MDCM [NOT PARALLELIZED]
  type(mDCMResults), public :: results
  

  type, public :: ModifiedDCM
     private
     ! Validates object as constructed:
     logical, private :: constructed = .FALSE.

     ! Controls messaging:
     type(mDCMIO), private :: io

     ! mDCM behavior options:
     type(mDCMOptions), private :: options

     ! Data objects:
     ! ---  N/A  --- (at this time)

     ! RNG:
     procedure(RANDOM), private, nopass, pointer :: rang => NULL()

     ! NOTE: Results and calculation data types are created in the main routine,
     !       in an attempt to minimize the state of the class

   contains
     private
     ! For simulations:
!     procedure, public  :: cascaw
!     procedure, public  :: simulate        => cascaw
!     procedure, public  :: execute         => cascaw
!     procedure, public  :: start           => cascaw
!     procedure, public  :: interact        => cascaw
!     procedure, public  :: collide         => cascaw
!     procedure, public  :: simulateINC     => cascaw
!     procedure, public  :: simulateCascade => cascaw

     ! For querying the object:
     procedure, public  :: properlyConstructed
     procedure, public  :: queryOptions
     procedure, public  :: queryRNG
     procedure, public  :: setOptions
     procedure, private :: validateOptions

     ! Internal to the simulation:
     ! ---  N/A  --- (for the time)

  end type ModifiedDCM



  private :: new_ModifiedDCM
  private :: new_Results
  public  :: setupMDCMReaction
  public  :: inigam0
  private :: printMDCM

contains

! ====================================================================

  ! Constructors:
  include "new_ModifiedDCM.f90"
  include "new_Results.f90"

  ! Querying the object:
  include "mDCMAccess.f90"

  ! For reaction setup:
  include "setupMDCMReaction.f90"
  include "inigam0.f90"

  ! For message handling:
  include "printMDCM.f90"

  ! For the simulation:
  ! --- N/A --- (at this time)

! ====================================================================

end module modifiedDCMClass
