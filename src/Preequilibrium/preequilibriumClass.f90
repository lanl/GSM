
! ======================================================================
!
! Module file for the simulation of preequilibrium decay
!
! To use this physics class, users must "use" the module, the data structures
! "residualNucleus" and "preequilibriumFragment", and the "Preequilibrium" class.
!
! Construction of a Preequilbrium class is described in the
! "new_Preequilibrium" function (other file).
!
! All I/O is handled  by the preeqiulibrium class by default, however users
! have the option to point the preequilibrium class to a difference I/O
! handler, having control over the I/O produced. I/O is currently limited
! to error, warnings, and comments.
!
!
! Model written and edited by many others prior to class implementation.
! Written by CMJ, XCP-3, August 2018
!
! ======================================================================

module preequilibriumClass

  use, intrinsic:: iso_fortran_env, only: int32, real64
  use fermiBreakUpClass,       only:  FermiBreakup
  use molnixClass,             only:  Molnix
  use fissionBarrierClass,     only:  FissionBarrier
  use gammaJClass,             only:  GammaJ
  use preequilibriumDataClass, only:  PreequilibriumData
  use preequilibriumData,      only:  numAllowedFragments

  implicit none
  private

  ! For message printing:
  integer(int32), private, parameter :: defaultPreeqVerbose = 4_int32
  integer(int32), public :: preeqVerbose = defaultPreeqVerbose


  ! Module parameters
  real(real64),   private, parameter :: defaultEmissionWidth =  0.40_real64   ! Default exponential term for preequilibrium emission probability
  real(real64),   private, parameter :: defaultR0Multiplier  =  1.50_real64   ! Default radius parameter
  integer(int32), private, parameter :: defaultExcludePreeq  =  0_int32       ! Default for inclusion/exclusion of preequilibrium physics
  integer(int32), private, parameter :: defaultNumPreeqType  = 66_int32       ! Default number of particles to consider for preequilibrium decay
  integer(int32), private, parameter :: defaultLevelDenParam = 12_int32       ! Default value for which parameters to use for level density calculation
  integer(int32), public,  parameter :: minNumPreeqType = 6_int32               ! Min number of preeq. particles to be considered
  integer(int32), public,  parameter :: maxNumPreeqType = numAllowedFragments   ! Max number of preeq. particles to be considered
  integer(int32), private, parameter :: numCalcElements = numAllowedFragments + 1   ! Number of elements in the "preequilibriumCalculation" data type arrays

  ! For progeny origins...
  integer(int32), public,  parameter :: preequilibriumProgeny  = 0_int32
  integer(int32), public,  parameter :: fermiBreakUpProgenyType = 1_int32

  ! For state of residual:
  integer(int32), public,  parameter :: residualState     = 0_int32
  integer(int32), public,  parameter :: compoundState     = 1_int32
  integer(int32), public,  parameter :: stableState       = 2_int32
  integer(int32), public,  parameter :: fermiBreakupState = 3_int32

  ! For divide by zero checks:
  real(real64), private, parameter :: div0Lim = 1.0d-15


  ! Class Interface
  public :: newPreequilibrium
  interface newPreequilibrium
     module procedure :: new_Preequilibrium
  end interface


  ! Results object Interface
  public :: newPreequilibriumResults
  interface newPreequilibriumResults
     module procedure :: new_PreequilibriumResults
  end interface


  ! Load data types used by the preequilibrium class
  include "preequilibriumTypes.f90"


  abstract interface
     ! Random Number Generator Interface
     function RANDOM() result(rang) BIND(C)
       use, intrinsic:: iso_C_binding, only: c_double
       implicit none
       real(c_double) :: rang
     end function RANDOM
     ! Interface for I/O handling
     subroutine IOHANDLER(verbosity, type, text)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in) :: verbosity
       integer(int32),   intent(in) :: type
       character(len=*), intent(in) :: text
     end subroutine IOHANDLER
     ! For photon emission durring decay:
     subroutine PHOTOEMISSION(fragment, residual)
       import preequilibriumFragment, residualNucleus
       type(preequilibriumFragment), intent(in   ) :: fragment
       type(residualNucleus),        intent(in   ) :: residual
     end subroutine PHOTOEMISSION
  end interface




  ! Preequilibrium Class
  type, public :: Preequilibrium
     private

     ! Flags if the object was constructed or not:
     logical, private :: constructed = .FALSE.


     ! For all message handling:
     type(preequilibriumIO), private :: io


     ! Class options/numerics/behaviors:
     type(preequilibriumOptions), private :: options


     ! For the required data:
     ! (non-reaction specific)
     type(Molnix),             private, pointer :: molEnergy   => NULL()
     type(FissionBarrier),     private, pointer :: fissBarr    => NULL()
     ! (reaction specific)
     type(PreequilibriumData), private, pointer :: preeqData   => NULL()


     ! Random Number Generator:
     procedure(RANDOM),        private, pointer, nopass :: rng => NULL()


     ! For all external models (Fermi Breakup, photon emission):
     type(FermiBreakup),       private :: fbuObj
     procedure(PHOTOEMISSION), private, pointer, nopass :: photonEmission    => NULL()
     logical,                  private                  :: usePhotonEmission =  .FALSE.


   contains
     private
     ! For starting a simulation:
     procedure, public  :: decay       => simulatePreequilibrium
     procedure, public  :: deexcite    => simulatePreequilibrium
     procedure, public  :: equilibrate => simulatePreequilibrium 
     procedure, public  :: execute     => simulatePreequilibrium
     procedure, public  :: simulate    => simulatePreequilibrium
     procedure, public  :: simulatePreequilibrium
     procedure, public  :: start       => simulatePreequilibrium

     ! For querying object:
     procedure, public  :: properlyConstructed
     procedure, public  :: queryOptions
     procedure, public  :: queryMolnix
     procedure, public  :: queryFissionBarrier
     procedure, public  :: queryData
     procedure, public  :: queryRNG
     procedure, public  :: queryFermiBreakUp
     procedure, public  :: queryPhotonEmissionUse
     procedure, public  :: queryPhotonEmission

     ! For the simulation itself:
     procedure, private :: deltam
     procedure, private :: equilibration
     procedure, private :: fam
     procedure, private :: fermiPreeqInterface
     procedure, private :: kinema
     procedure, private :: peqemt
     procedure, private :: preqaux
     procedure, private :: rotor
!     procedure, private :: setALJ
     procedure, private :: trans8
     procedure, private :: validateOptions   ! Validates simulation options
  end type Preequilibrium


  public  :: preequilibriumInit
  private :: printPreeq
  private :: new_Preequilibrium
  private :: new_PreequilibriumResults
  private :: setALJ

contains
! ======================================================================

  ! Initialize preequilibrium model
  include "preequilibriumInit.f90"          ! Initializes preequilbrium data

  ! Class constructor
  include "new_Preequilibrium.f90"          ! Constructs a preequilibrium object
  include "new_PreequilibriumResults.f90"   ! Constructs a results object
  include "validateOptions.f90"             ! Validates options specified by the client

  ! Class access:
  include "preeqAccess.f90"

  ! For simulation:
  include "deltam.f90"                   ! Mass excesses
  include "equilibration.f90"            ! Performs equilbration of excited residual nucleus
  include "fam.f90"                      ! Obtains level density parameter
  include "fermiPreeqInterface.f90"      ! Fermi Break-up interface for preequilibrium physics
  include "kinema.f90"                   ! Kinematics routine
  include "peqemt.f90"                   ! Simulates single emission during preequilibrium decay
  include "preqaux.f90"                  ! Sets up auxiliary quantities for preequilibrium emission
  include "rotor.f90"                    ! Rotation of coordinate system(?)
  include "setALJ.f90"                   ! Sets up exciton phase factors for preequilbrium emission
  include "simulatePreequilibrium.f90"   ! Sets up residual and calls "equilibrate"
  include "trans8.f90"                   ! Exciton transition rates

  ! For printing to terminal (NOT class specific)
  include "printPreeq.f90"

! ======================================================================
end module preequilibriumClass
