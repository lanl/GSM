

! ==============================================================================
!
! Module for the simulation for the coalescence of lightly bound nucleons 
!    around an excited, residual nucleus (typically from the IntraNuclear Cascade
!    process)
!
! To use this physics class, users must 'use' the module, the data structures
!    below (coales_pmemo, coales_imemo), the class itself (Coalescence), and the
!    class constructor (new_Coalescence). Details on creating an instance of
!    the coalescence class are shown in the "new_Coalescence" function.
!
! To modify how much information is printed to terminal, change value of
! 'coalesVerbose' (internal only). Current printout include warnings, errors, and
! comments (no particle information)
!
!
! Written by CMJ, XCP-3, July 2018 (Initial coalescence class creation)
!
! ==============================================================================

module coalescenceClass

  ! Parameters used by class
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use coalescenceParams, only : zro, two

  implicit none
  private

  ! For verbosity of messages:
  integer(int32), private, parameter :: defaultCoalesVerbose = 4_int32
  integer(int32), public :: coalesVerbose = defaultCoalesVerbose

  ! Default behaviors (i.e. options) used by all Coalescence objects
  integer(int32), private, parameter :: defaultExpandedCoalescence = 2_int32  ! Use (>0) or not use (<=0) expanded Coalescence model
  ! Default coalescence radii and transition energy range
  real(real64),   private, parameter :: lowerRadiiEBound = 0.3_real64   ! Lower bound at which to use the other default coalescence radii
  real(real64),   private, parameter :: upperRadiiEBound = 1.0_real64   ! Upper bound at which to use the other default coalescence radii
  real(real64),   private, parameter, dimension(2) :: &   ! Coalescence radii have units [GeV/c]
       & defaultCRDeut  = [ 0.090_real64, 0.150_real64 ], &
       & defaultCRTrit  = [ 0.108_real64, 0.175_real64 ], &
       & defaultCRAlpha = [ 0.130_real64, 0.205_real64 ], &
!       & defaultCRAlpha = [ 0.115_real64, 0.175_real64 ], & ! These values result in NO coalescence of alpha particles (and thus LF)
       & defaultCRLFrag = [ 0.175_real64, 0.250_real64 ]
  ! For "error protection" (mathematical errors)
  real(real64),   private, parameter :: div0Limit = 1.0d-20
  real(real64),   private, parameter :: sqrRootCorrection = 0.01d0

  ! Coalescence Class contructor
  public :: newCoalescence
  interface newCoalescence
     module procedure :: coalesConstructorMain
  end interface


  ! Coalescence Results contructor
  public :: newCoalescenceResults
  interface newCoalescenceResults
     module procedure :: coalesResultsConstructorMain
  end interface


  ! Coalescence Data Class contructor
  public :: newCoalescenceData
  interface newCoalescenceData
     module procedure :: coalesDataConstructorMain
  end interface


  ! I/O Handler interface
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


  ! Class data types
  include "coalesClassTypes.f90"


  ! Coalescence class
  type, public :: Coalescence
     private
     ! Ensures the object is constructed before allowing the simulation to proceed
     logical, private :: constructed = .FALSE.


     ! Handles all I/O in the class (messages and printing)
     type(coalescenceIO),      private :: io


     ! Contains all options used by the coalescence class
     type(coalescenceOptions), private :: options


     ! Data to be used with the class
     type(coalescenceData),    private :: data


     ! NOTE: Results for the coalescence class are NOT class-specific
     !       to minimize the class's state

   contains
     private
     ! For starting simulation:
     procedure, public  :: simulateCoalescence
     procedure, public  :: simulate  => simulateCoalescence
     procedure, public  :: execute   => simulateCoalescence
     procedure, public  :: start     => simulateCoalescence
     procedure, public  :: coalesce  => simulateCoalescence
     procedure, public  :: compound  => simulateCoalescence

     ! Access to class objects:
     procedure, public  :: properlyConstructedObject
     procedure, public  :: properlyConstructed => properlyConstructedObject
     procedure, public  :: queryOptions
     procedure, public  :: queryData

     ! For simulation:
     procedure, private :: coalesl              ! Performs simulation (main routine)
     procedure, private :: codir                ! Obtains coalesced particle direction
     procedure, private :: kinema               ! Kinematics Routine (performs Lorentz transformation)
!     procedure, private :: rijm                 ! Obtains distance of closest approach for 2 fragments
  end type Coalescence


  ! Procedures external to the coalescence class
  private :: coalesConstructorMain
  private :: coalesResultsConstructorMain
  private :: coalesDataConstructorMain

contains

! ======================================================================

  ! Class constructors
  include "coalesConstructorMain.f90"
  include "coalesResultsConstructorMain.f90"
  include "coalesDataConstructorMain.f90"

! ======================================================================

  ! Main calculation routine
  include "simulateCoalescence.f90"
  include "coalesl.f90"

  ! Class access:
  include "coalesAccess.f90"
  include "coalesDataAccess.f90"

  ! To obtain direction of coalescence particle
  include "codir.f90"

  ! Does some kinematics
  include "kinema.f90"

  ! Controls all of the classes I/O
  include "printCoales.f90"

  ! Distance of closest approach for 2 particles? [Not used]
!  include "rijm.f90"

! ======================================================================

end module coalescenceClass
