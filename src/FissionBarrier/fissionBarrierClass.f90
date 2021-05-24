
! ======================================================================
!
! Module file for the obtaining of fission barriers.
!
! To use this physics class, users must "use" the module and the "Molnix"
! data structure/constructor below.
!
! Construction of a FissionBarrier object is discussed in the function 
! "new_FissionBarrier".
!
!
! Model written and edited by many others prior to class creation.
! Class written by CMJ, XCP-3, 8/2018 (Initial class creation)
!
! ======================================================================

module fissionBarrierClass

  use, intrinsic:: iso_fortran_env, only: int32, real64
  use molnixClass, only: Molnix

  implicit none
  private

  ! For printing:
  integer(int32), private, parameter :: defaultFBVerbose = 4_int32
  integer(int32), public :: fissBarrVerbose = defaultFBVerbose

  ! Class defaults
  real(real64),   private, parameter :: defaultR0m   = 1.2_real64


  ! For divide by zero limits:
  real(real64),   private, parameter :: div0Lim = 1.0d-15


  ! Polynomial orders
  integer(int32), private, parameter :: npa   = 5_int32
  integer(int32), private, parameter :: npz   = 8_int32
  integer(int32), private, parameter :: npl   = 5_int32
  real(real64),   private, parameter :: amm   = 0.704_real64


  ! Include data
  integer(int32), private, parameter :: fbDataDim = 51_int32
  include "camer70.f90"


  ! Class constructor
  public :: newFissionBarrier
  interface newFissionBarrier
     module procedure :: new_FissionBarrier
  end interface newFissionBarrier


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


  ! Provides numerics to the class:
  type, public :: fissionBarrierOptions
     private
     real(real64),   public :: r0m = defaultR0m
  end type fissionBarrierOptions


  ! Provides message handling to the class:
  type, private :: fissionBarrierIO
     private
     character(LEN=512),   private                  :: message =  ""
     procedure(IOHANDLER), private, pointer, nopass :: print   => printFissionBarrier
  end type fissionBarrierIO


  ! Class definition
  type, public :: FissionBarrier
     private

     ! To verify construction:
     logical,                     private :: constructed = .FALSE.

     ! Handles all messages:
     type(fissionBarrierIO),      private :: io

     ! Provides numerics:
     type(fissionBarrierOptions), private :: options

     ! For data:
     type(Molnix),                private,  pointer :: fbMolnix =>  NULL()

   contains
     private
     ! Client usable procedures:
     procedure, public  :: bf
     procedure, public  :: barfit
     procedure, public  :: bsfit

     ! Data procedures/interfaces:
     procedure, private :: bf2
!     procedure, private :: elmaxc
!     procedure, private :: lpoly2
     procedure, private :: subev

     ! Class member access function:
     procedure, public  :: properlyConstructed
     procedure, public  :: queryOptions
     procedure, public  :: queryMolnix

  end type FissionBarrier

  private :: elmaxc
  private :: lpoly2

contains
  ! ======================================================================

  ! Class constructor
  include "new_FissionBarrier.f90"

  ! Main top-level routine
  include "bf.f90"

  ! Majority of routines
  include "barfit.f90"

  ! Quadratic interpolation of values
  include "subev.f90"

  ! For client usage (verifies construction, returns 'options' data type, etc.)
  include "dataAccess.f90"

  ! For printing of messages:
  include "printFissionBarrier.f90"

! ======================================================================
end module fissionBarrierClass
