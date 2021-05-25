! ===================================================================================
!
! This module contains the class for the evaporation/fission model, GEM2.
!
! To use this physics class, users must 'use' the module. Explicitly, users MUST use
! the interface and data type 'Evaporation' and 'evaporationFragment'. Physics
! simulation is done by invoking "call myEvaporationObject%simulate()".
! Users have the option to also use the data types 'evaporationFragment', 'residual',
! and 'fissionFragment'.
!
! Results from the evaporation class simulation are obtained by 'getNumProgeny',
! 'getResidualExists', 'getResidual' (if it exists), 'getFissionOccurred', 'getFissionFragments'
! (if exists), 'getPreFissionNucleus' (if exists), and 'getFissionProbability'.
!
! Construction of an Evaporation class is described later in the 'new_Evaporation' function.
!
! Verbosity of class is determined below by the 'evapVerbose' parameter. The default value
! is 2 (developers may change later if desired). With this verbosity level, all errors,
! warnings, and comments of sufficient worth (regarding class construction) are printed.
! By default, NO divide by zero errors or square root errors are printed.
!
!
! Model written and edited by many others prior to class creation.
! Written by CMJ, XCP-3, July 2018 (evaporation/fission class creation).
!
! ===================================================================================

module evaporationClass

  ! Parameters used by the evaporation/fission model (GEM2)
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use evaporationParams, only: zro, one, two
  use evaporationFissionData, only: maxSize, aMax, iiz, inn
  use molnixClass, only: Molnix
  use fermiBreakUpClass, only: FermiBreakup
  use evaporationDataClass, only: EvaporationData

  implicit none
  private


  ! For message verbosity:
  integer(int32), private, parameter :: defaultEvapVerbose = 4_int32
  integer(int32), public :: evapVerbose = defaultEvapVerbose

  ! Flags for compound, fission fragments, and progeny physics underwent/origin:
  integer(int32), private, parameter :: fissionFragFlag   = -1_int32   ! Starting flag for non-existant fission fragments
  integer(int32), public,  parameter :: evaporationFlag   =  0_int32   ! Flags particle underwent or originated from Evaporation
  integer(int32), public,  parameter :: fermiBreakUpFlag  =  1_int32   ! Flags particle underwent or originated from Fermi Breakup
  integer(int32), public,  parameter :: fissionFlag       =  2_int32   ! Flags particle underwent or originated from Fission


  ! Default calculation options
  integer(int32), private, parameter    :: defaultNumEvapType      = 66_int32    ! Number of particles considered for evaporation
  integer(int32), private, parameter    :: defaultFissParameter    =  0_int32    ! =0 (Reevaluation parameters in the RAL fission model), else (Original parameters)
  integer(int32), private, parameter    :: defaultRedCompTime      =  0_int32    ! =0 (Normal calculation), Else (CPU time saving calculation)
  real(real64),   private, parameter    :: defaultInverseParameter =  0.0_real64   ! Parameter for inv. react. cross section; =0 (Dostrovsky+Matsuse), =10 (Simple parameter set w/ r0=1.5), 0 < rcal < 10 (Simple parameter set w/ r0=rcal)
  real(real64),   private, parameter    :: defaultAfMultipler      =  1.0_real64   ! Multipler for A_f
  real(real64),   private, parameter    :: defaultCzMultipler      =  1.0_real64   ! Multipler for C(Z)
  ! For class procedures:
  real(real64),   private, parameter    :: defaultRRVcoul          =  1.5_real64   ! Default inverse cross section parameter used in "vcoul" procedure

  ! Valid range of evaporated particles
  integer(int32), private, parameter    :: maxNumEvapType = 66_int32
  integer(int32), private, parameter    :: minNumEvapType =  1_int32


  ! Valid range of InverseParameter option:
  real(real64),   private, parameter    :: maxInverseParameter = 10.0_real64


  ! For divide by zero errors:
  real(real64), private, parameter :: div0Lim = 1.0d-15



  ! Routines accessible by outside calls
  public :: setLevelDensity


  ! Obtain the data types used by the evaporation class
  include "evaporationTypes.f90"


  ! Procedure pointer interface for message handling, random number generator, and photon emission
  abstract interface
     ! Interface for I/O handling
     subroutine IOHANDLER(verbosity, type, text)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in) :: verbosity
       integer(int32),   intent(in) :: type
       character(len=*), intent(in) :: text
     end subroutine IOHANDLER
     function RANDOM() result(rndm) BIND(C)
       use, intrinsic:: iso_C_binding, only: c_double
       real(c_double) :: rndm
     end function RANDOM
     subroutine PHOTOEMISSION(fragment, parentA, parentZ, parentKinE)
       use, intrinsic:: iso_fortran_env, only: real64
       import evaporationFragment
       implicit none
       type(evaporationFragment), intent(in   ) :: fragment
       real(real64),              intent(in   ) :: parentA
       real(real64),              intent(in   ) :: parentZ
       real(real64),              intent(in   ) :: parentKinE
     end subroutine PHOTOEMISSION
  end interface


  ! Evaporation results object constructor
  public :: newEvaporationResults
  interface newEvaporationResults
     module procedure :: new_EvaporationResults
  end interface newEvaporationResults


  ! Evaporation class constructor
  public :: newEvaporation
  interface newEvaporation
     module procedure :: new_Evaporation
  end interface newEvaporation


  ! Class object
  type, public :: Evaporation
     private

     ! Flag to ensure class was constructed prior to simulation
     logical,                  private :: constructed = .FALSE.


     ! Handles all messages
     type(evaporationIO),      private :: io


     ! Simulation options/behaviors
     type(evaporationOptions), private :: options
     real(real64),             private :: rrVcoul = defaultRRVcoul

     ! Random number generator:
     procedure(RANDOM),        private, pointer, nopass :: rang => NULL()


     ! External models (fermi break up, photon emission, etc.)
     type(FermiBreakup), private :: fbuObj
     procedure(PHOTOEMISSION), private, pointer, nopass :: photonEmission => NULL()
     logical, private :: usePhotoEmission = .FALSE.


     ! Simulation data:
     type(Molnix),             private, pointer :: evapMolnix => NULL()
     type(EvaporationData),    private, pointer :: data => NULL()

   contains
     private
     ! To start a simulation:
     procedure, public  :: simulate  => gemdec
     procedure, public  :: execute   => gemdec
     procedure, public  :: start     => gemdec
     procedure, public  :: evaporate => gemdec
     procedure, public  :: decay     => gemdec
     procedure, public  :: deexcite  => gemdec
     procedure, public  :: GEMDecay  => gemdec

     ! For client access to class members:
     procedure, private :: validateOptions      ! Ensures the options passed by the client are all valid
     procedure, public  :: properlyConstructed
     procedure, public  :: queryOptions
     procedure, public  :: queryRNG
     procedure, public  :: queryMolnix
     procedure, public  :: queryData
     procedure, public  :: queryFermiBreakup
     procedure, public  :: queryPhotoEmissionUse
     procedure, public  :: queryPhotoEmission

     ! For the simulation itself:
     procedure, private :: gemdec
     procedure, private :: calr
     procedure, private :: cam
!     procedure, private :: ckcal
!     procedure, private :: dost
     procedure, private :: drein1
     procedure, private :: drein2
     procedure, private :: efms
     procedure, private :: energy
     procedure, private :: ey2
     procedure, private :: ey3
     procedure, private :: eye
     procedure, private :: evapFermiInterface
     procedure, private :: fisgem
     procedure, private :: fprob
     procedure, private :: gamma
     procedure, private :: gauss2
     procedure, private :: pe
     procedure, private :: progenyMass
!     procedure, private :: radgem
     procedure, private :: rb
     procedure, private :: rho
     procedure, private :: selecte
     procedure, private :: stdcay
     procedure, private :: vcoul
  end type Evaporation



  ! Constructor(s):
  private :: new_Evaporation
  private :: new_EvaporationResults


  ! Misc.:
  private :: f1
  private :: f2
  private :: printEvaporation

contains

! ===================================================================================

  ! Constructor
  include "new_Evaporation.f90"
  include "new_EvaporationResults.f90"
  include "validateOptions.f90"

  ! Accessing members of the evaporation object:
  include "evapAccess.f90"

  ! Sets level densities (available to user and class)
  include "setLevelDensities.f90"

  ! Routines for the simulation
  include "drein.f90"
  include "energy.f90"
  include "eye.f90"
  include "fermiEvapInterface.f90"
  include "fisgem.f90"
  include "fita.f90"
  include "fprob.f90"
  include "gamma.f90"
  include "gauss2.f90"
  include "gem2.f90"
  include "gemdec.f90"
  include "geta.f90"
  include "selecte.f90"
  include "stdcay.f90"
  include "vcoul.f90"

  ! Handles the printing of all messages
  include "printEvaporation.f90"

! ===================================================================================

end module evaporationClass
