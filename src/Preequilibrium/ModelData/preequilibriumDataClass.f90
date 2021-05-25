
! ======================================================================
!
! Module file for the preequilbrium DATA class
!
! This class is used HEAVILY by the preeuqilibrium class. It contains
! within it:
! (1) Ground state and Fission barrier energies
! (2) The class for GammaJ [fitting coefficients]
! (3) The class for Fission Barrier calculations
!    (3a) The Molnix class (obtains mass excess, pairing gap energies,
!         and ground state energy corrections (microscopic)
!
! To use this class, users must "use" the module and the class/constructor
! "PreequilibriumData".
!
! Construction of the class is described in further detail below, in the
! constructor "new_PreequlibriumData"
!
!
!
! Written by CMJ, XCP-3, 8/2018
!
! ======================================================================

module preequilibriumDataClass

  use, intrinsic:: iso_fortran_env, only: int32, real64
  use molnixClass,           only:   Molnix
  use fissionBarrierClass,   only:   FissionBarrier
  use gammaJClass,           only:   GammaJ

  implicit none
  private

  ! For printing:
  integer(int32), private, parameter :: defaultPreeqDataVerbose = 4_int32
  integer(int32), public :: preeqDataVerbose = defaultPreeqDataVerbose

  ! For array sizes:
  integer(int32), parameter, private :: dataEgsEbDim1 =  30_int32   ! Regarding  Z index
  integer(int32), parameter, private :: dataEgsEbDim2 =  70_int32   ! Regarding  A index
  integer(int32), parameter, private :: dataEgsEbDim3 = 100_int32   ! Regarding ln index [angular momentum quantum number]


  ! Procedure pointer interface
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



  ! Class constructor/interface
  public :: newPreequilibriumData
  interface newPreequilibriumData
     module procedure :: new_PreequilibriumData
  end interface newPreequilibriumData


  ! For message handling:
  type, private :: preequilibriumDataIO
     private
     character(LEN=512),   private                  :: message =  ""
     procedure(IOHANDLER), private, pointer, nopass :: print   => printPreeq
  end type preequilibriumDataIO


  ! Contains ground state and fission barrier energies for a compound nucleus:
  type, private :: compoundEnergies
     private
     real(real64), private, dimension(dataEgsEbDim1, dataEgsEbDim2, dataEgsEbDim3) :: eb    = 0.0_real64   ! Fission barrier height for all possible nuclei
     real(real64), private, dimension(dataEgsEbDim1, dataEgsEbDim2, dataEgsEbDim3) :: egs   = 0.0_real64   ! Ground state energy    for all possible nuclei
  end type compoundEnergies


  ! Contains the A/Z/kinEnergy pairing for the compound
  type, private :: preequilibriumCompound
     private
     real(real64), private :: numBaryons = 0.0_real64   ! Total neutrons present
     real(real64), private :: numProtons = 0.0_real64   ! Total protons  present
     real(real64), private :: kinEnergy  = 0.0_real64   ! Initial energy [GeV]
  end type preequilibriumCompound



  ! Preequilibrium data class
  type, public :: PreequilibriumData
     private

     ! To ensure data class was constructed before use:
     logical, private :: constructed = .FALSE.


     ! For message handling:
     type(preequilibriumDataIO), private :: io


     ! Non-Specific Data:
     type(Molnix),         private, pointer  :: molEnergy => NULL()   ! For energies (related to masses)
     type(FissionBarrier), private, pointer  :: fissBarr  => NULL()   ! For the fission barrier


     ! Reaction Specific Data:
     type(preequilibriumCompound), private :: compound
     type(compoundEnergies),       private :: compEnergy   ! Ground state and fission barrier energies for the compound
     type(GammaJ),                 public  :: gammaJObj    ! For F_j coefficients


   contains
     private

     ! To obtain information on the compound nucleus that was used to establish the data:
     procedure, public  :: numBaryons
     procedure, public  :: numProtons
     procedure, public  :: kinEnergy


     ! To obtain ground state and fissin barrier energies of the compound:
     procedure, public :: auxl


     ! To obtain the values of the ground state and fission barrier arrays:
     procedure, private :: checkIndex   ! Verifies that the requested valued of 'eb' or 'egs' exists
     procedure, public  :: eb
     procedure, public  :: egs


     ! Various ways clients can access the information here:
     procedure, public  :: properlyConstructed
     procedure, public  :: getMolnix
     procedure, public  :: getFissionBarrier

  end type PreequilibriumData
     

contains


! ======================================================================

  ! Class consructor(s):
  include "new_PreequilibriumDataClass.f90"


  ! Auxl routine (book keeping)
  include "auxl.f90"


  ! So clients can access the data object's information:
  include "dataAccess.f90"


  ! For Handling I/O
  include "printPreeq.f90"


! ======================================================================
end module preequilibriumDataClass
