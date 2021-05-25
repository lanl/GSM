
! ======================================================================
!
! Module for GSM Class
!
! Written by CMJ, XCP-3, August 2018 - Decemeber 2018
!
! ======================================================================

module GSMClass

  use, intrinsic:: iso_fortran_env, only: int32, int64, real64

  ! General information for classes
  use molnixClass,             only: Molnix, molnixOptions
  use fissionBarrierClass,     only: FissionBarrier, fissionBarrierOptions

  ! Pull in data objects for each submodel
  use StandardDCMDataClass,     only: StandardDCMData, sDCMDataOptions
!  use ModifiedDCMDataClass,     only: ModifiedDCMData
  use PreequilibriumDataClass,  only: PreequilibriumData
  use EvaporationDataClass,     only: EvaporationData, evaporationDataOptions

  ! Pull in model object (physics class) and their options type
  use standardDCMClass,     only: standardDCM, sDCMOptions, sDCMProjectile
!  use modifiedDCMClass,     only: modifiedDCM, mDCMOptions
  use coalescenceClass,     only: Coalescence, coalescenceData, &
       & coalescenceOptions, coalescenceDataOptions
  use preequilibriumClass,  only: Preequilibrium, preequilibriumOptions
  use evaporationClass,     only: Evaporation, evaporationOptions
  use fermiBreakUpClass,    only: FermiBreakup, fermiBreakUpOptions

  ! Import objects for GSM
  use BremsPhotonMod,  only: BremsPhoton
  use EventDataMod,    only: EventData
  use OutputDataMod,   only: OutputData

  implicit none
  private

  ! Iterator for parameterized arrays
  integer(int32), private :: iter = 0_int32

  ! Include all parameters for object defaults:
  include "defaultValues.f90"

  ! Non-model parameters:
  integer(int32), public :: gsmVerbose = defaultGSMVerbose


  ! GSM class constructor
  public :: newGSM
  interface newGSM
     module procedure :: new_GSM
  end interface

  ! GSM Results class constructor
  public :: newGSMResults
  interface newGSMResults
     module procedure :: gsmResultsMainConstructor
  end interface

  ! Include the objects utilized for a simulation
  include "objects/definitionsMain.f90"

  ! Procedure pointer interface (for I/O and RNG)
  abstract interface
     ! Interface for I/O handling
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
  include "photonEmissionInterfaces.f90"
  logical, private, protected :: useGammaCascade = .FALSE.
  procedure(preequilibriumPHOTONEMISSION), private, protected, pointer :: preeqGammaCascade => preeqPhotonEmission
  procedure(evaporationPHOTONEMISSION),    private, protected, pointer :: evapGammaCascade  => evapPhotonEmission
  procedure(PHOTOEMISSION),                private, protected, pointer :: gammaCascade      => NULL()


  ! GSM Class
  type, public :: GSM
     private

     ! Flag for proper construction
     logical, private :: constructed = .FALSE.
     type(GSMIO), private :: io

     ! GSM Class Simulation Options:
     type(GSMOptions), private  :: options

     ! RNG for GSM and its sub-models
     logical, private :: useDefaultRNG = .TRUE.
     procedure(RANDOM), private, nopass, pointer :: rang => NULL()

     ! Model data and model objects
     type(generalPhysicsData),   private :: genData
     type(generalPhysicsModels), private :: genModels

     ! For possible photon emission (necessary for consistent use of photon emission in evaporation/preequilibrium):
     logical, private :: usePhotonEmission = .FALSE.

     ! NOTE: Results and calculation data types are created in the main routine,
     !       minimizing the state of the GSM class

   contains
     private

     ! Establishing/changing state of GSM (without constructor):
     ! ------------------------
     procedure, public  :: updateOptions
     procedure, private :: validateOptions
     procedure, private :: validateGSMState
     procedure, public  :: updateRNG
     procedure, public  :: updateMessageHandler

     ! For querying the state of GSM:
     ! ------------------------
     procedure, public  :: properlyConstructed
     procedure, public  :: queryOptions
     procedure, public  :: queryRNG
     procedure, public  :: queryPhotoEmissionUse

     ! For generating an output file:
     ! ------------------------
     procedure, public  :: generateOutput
     procedure, public  :: output    => generateOutput
     procedure, private :: gsmMain   ! Procedure interfaced to for output generation

     ! For simulating an event:
     ! ------------------------
     procedure, public  :: simulateEvent
     procedure, public  :: performSimulation => simulateEvent
     procedure, public  :: collide  => simulateEvent
     procedure, public  :: interact => simulateEvent

     ! General procedures for simulating (output generation or event generation)
     ! ------------------------
     generic,   public  :: simulate => generateOutput, simulateEvent
     generic,   public  :: execute  => generateOutput, simulateEvent
     generic,   public  :: start    => generateOutput, simulateEvent

     ! Physics interfaces:
     ! ------------------------
     procedure, public  :: standardDCMInterface
     procedure, private :: setMDCMReaction   ! Establishes mDCM defaults and the nuclei for an event
     procedure, public  :: modifiedDCMInterface
     procedure, public  :: coalescenceInterface
     procedure, public  :: fermiBreakUpInterface
     procedure, public  :: preequilibriumInterface
     procedure, public  :: evaporationInterface
     procedure, private :: simulateDecay   ! Performs preeq./evap. physics

     ! Misc. (clients may test the results of GSM internally this way to see its behaviors)
     ! ------------------------
     procedure, public  :: formNuclei
     procedure, public  :: buildNuclei => formNuclei
     procedure, private :: swapNuclei
     procedure, private :: updateLevelDensities
     procedure, public  :: inquireINCModel
     procedure, private :: sampleEnergy   ! Used for sampling transition energy (when smooth desired)

     ! =======================================================
     ! PRIVATE PROCEDURES INTERNAL TO GSM:
     ! =======================================================
     ! For the GSMReaction object (convenient way to carry the reaction info around):
     ! ------------------------
     procedure, private :: setupReaction   ! NOTE: This is the GSMReaction constructor! (needs GSM internals)
     procedure, private :: constructGeneralData
     procedure, private :: constructGeneralPhysics
     procedure, private :: constructSpecificData
     procedure, private :: constructSpecificPhysics
     procedure, private :: constructGeneralModels

     ! For simulation of spallation physics:
     ! ------------------------
     procedure, private :: renorm   ! Renormalizes p/E for residual target (sDCM only)
     procedure, private :: checkMomentum
     procedure, private :: simulateEvents
     procedure, private :: eventLoop
     procedure, private :: restor1
     procedure, private :: ststcs

     ! For event accounting and accumulation (tallying):
     ! ------------------------
     procedure, private :: vlobd
     procedure, private :: resdist
     procedure, private :: opandist
     procedure, private :: disnmul

     ! For output information:
     ! ------------------------
     procedure, private :: prinp     ! Prints license, etc.
     procedure, private :: initial   ! Resets output file variables
     procedure, private :: typeout   ! Prints output (main procedure for this)
     ! Regarding inverse cross section calculation:
     procedure, private :: sinvla
     procedure, private :: xabs
     procedure, private :: cskalb
     procedure, private :: radius
     ! Regarding printing of distributions:
     procedure, private :: prtmult
     procedure, private :: prtdadz
     procedure, private :: prrdis
     procedure, private :: propan
     procedure, private :: pdisnm
     procedure, private :: prtdist
     ! Regarding pisa information:
     procedure, private :: pisaInit
     procedure, private :: PisaSpectra
     procedure, private :: pisaPrint

     ! Misc. Procedures:
     ! ------------------------
     procedure, private :: printPartProp   ! For printing particle properties [for debugging]

  end type GSM


  ! For including photon emission physics
  ! NOTE: This, in its current implementation, must be GENERAL FOR ALL GSM OBJECTS
  public  :: setPhotonEmission
  private :: preeqPhotonEmission
  private :: evapPhotonEmission

  ! Non-GSM Class procedures
  private :: new_GSM          ! Returns a 'GSM' object
  private :: gsmResultsMainConstructor   ! Returns a 'results' object
  private :: printGSM         ! Filters and prints all messages as appropriate

  ! Misc.
  public  :: determineRestMass

contains

! ======================================================================

  ! Include procedures for all sub-objects of GSM
  include "objects/functionsMain.f90"


  ! For establishing/changing/querying the state of GSM:
  ! ------------------------
  include "new_GSM.f90"                    ! GSM Constructor
  include "dataAccess.f90"                 ! Contains most query/update proceduress
  include "validateOptions.f90"            ! Corrects all invalid options set
  include "validateGSMState.f90"           ! Validates the state of GSM and the results object
  include "validateOutputObj.f90"          ! Validates the GSMOutput object
  include "constructGeneralData.f90"       ! Constructs general  physics data  objects
  include "constructGeneralPhysics.f90"    ! Constructs general  physics model objects
  include "constructGeneralModels.f90"     ! Calls to both general data/physics constructors

  ! Regarding the output file:
  ! ------------------------
  include "generateOutput.f90"   ! Interfaces to "gsmMain"
  include "gsmMain.f90"          ! Main procedure to generate output files
  ! (regarding physics):
  include "simulateEvents.f90"   ! Simulates all desired events
  ! (regarding variable modification):
  include "initial.f90"          ! Resets all variables used by the output file
  include "ststcs.f90"           ! Accumulates statistics at various phases of the reaction
  ! (regarding output objects):
  include "outputObjects/main.f90"
  ! (regarding event tallying):
  include "tally/main.f90"
  ! (regarding writing of the output):
  include "output/main.f90"

  ! For simulating an event:
  ! ------------------------
  include "gsmResultsMainConstructor.f90"   ! GSM Results Constructor
  include "simulateEvent.f90"   ! Main procedure to simulate a single event (client usage ONLY)
  include "eventLoop.f90"       ! Main event simulating procedure
  include "validResidual.f90"   ! Determines if a residual is valid or not
  include "determineRestMass.f90" ! Determines rest mass of a nucleus (given A/Z)
  include "formNuclei.f90"      ! Forms/establishes values for the projectile and target objects
  include "swapNuclei.f90"      ! Swap a projectile and target nucleus
  include "updateLevelDensities.f90"       ! Updates Af and Cz multipliers
  include "inquireINCModel.f90" ! Inquires the INC model to be used for the reaction (sDCM or mDCM)
  include "sampleEnergy.f90"    ! Samples a value around a provided central value and std. dev.
  include "renorm.f90"          ! Renormalizes post sDCM physics for a valid event
  include "restor1.f90"         ! Adds the stable residual to the progeny array
  include "checkMomentum.f90"   ! Checks momentum conservation for all progeny/residuals

  ! Physics interfaces:
  ! ------------------------
  include "standardDCMInterface.f90"      ! INC physics (nucleon or smaller projectile)
  include "modifiedDCMInterface.f90"      ! INC physics (any projectile, large kinEnergy)
  include "setMDCMReaction.f90"           !    Sets mDCM defaults/nuclei information
  include "coalescenceInterface.f90"      ! Coalescence physics (i.e. formation of fast compounds)
  include "fermiBreakUpInterface.f90"     ! Fermi BreakUp physics (i.e. disintegration)
  include "preequilibriumInterface.f90"   ! Preequilibrium physics (i.e. equilibration)
  include "evaporationInterface.f90"      ! Evaporation/fission physics (i.e. de-excitation)
  include "simulateDecay.f90"             ! Performs preeq./evap. physics
  include "photonEmission.f90"            ! Interfaces to a client photon emission procedure

  ! For the GSMReaction object:
  ! ------------------------
  include "resetEvent.f90"                 ! Resets the data from an event in the results object
  include "setupReaction.f90"              ! Establish all values needed prior to a reaction
  include "constructSpecificData.f90"      ! Constructs specific physics data  objects
  include "constructSpecificPhysics.f90"   ! Constructs specific physics model objects

  ! Misc. class procedures:
  ! ------------------------
  include "residualMomenta.f90"   ! Calculates total momentum of a residual (linear and angular)

  ! All other procedures:
  ! ------------------------
  include "printPartProp.f90"   ! For debugging (prints all particle information)
  include "printGSM.f90"        ! Handles all I/O (messages)

! ======================================================================

end module GSMClass
