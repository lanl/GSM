
! ==============================================================================
!
! Data type for the calculation options of GSM (when generating an output file)
! and for the variables utilized to generate the output file
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

  ! Establishes a histogram bin
  !    (note negative value uses GSM defaults)
  type, private :: histogramBin
     private
     real(real64), public :: lowerBound = -1.0_real64
     real(real64), public :: upperBound = -1.0_real64
  end type histogramBin

  ! Include all default values:
  include "defaultOutputOptions.f90"

  ! Data type containing all of the options used to create the output file
  type, public :: GSMOutput
     private

     ! File names
     character(2*nameLength), public :: inputFile     = defaultInpFile
     character(  nameLength), public :: outputFile    = defaultOutFile
     character(  nameLength), public :: auxiliaryFile = defaultAuxFile

     ! Date of initial write to the output file
     character(24), private :: date = ""

     ! Input comments
     integer(int32), public :: numComments = 0_int32
     character(nameLength), public, dimension(maxCommentLines) :: comments = ""

     ! Number of inelastic events to simulate
     integer(int64), public :: numInelasticEvents = defaultInelasticEvents
     integer(int64), private:: maxEventAttempts = 10_int64 * defaultInelasticEvents

     ! Calculation options
     integer(int32), public :: &
          & angularSpectra      = defaultAngSpec, & ! Flags to calculate angular spectra [0, 1, 2] MANG
          & energySpectra       = defaultEnergySpec, & ! Flags to calculate energy spectra [0, 1, 2] MSPEC
          & doubleDiffSpectra   = defaultDoubDiffSpec, & ! Flags to calculate double diff. spectra [0, 1, 2] MDUBL
          & multiplicities      = defaultMultiplicity, & ! Flags to calculate particle multiplicities [0, 1, 2] MPYLD
          & channelCrossSection = defaultChannelXSCalc, & ! Flags to calculate reaction channel cross sections [0, 1] MCHY
          & nuclideCrossSection = defaultNuclideXSCalc    ! Flags to calculate nuclide cross section [0, 1, 2, 3] MISY
     logical, public :: &
          & printPISA           = defaultPISAPrint, & ! Flags to print PISA cross section tables
          & printHist           = defaultHistPrint    ! Flags to print HIST residual tables

     ! Angular bins for spectra calculations
     type(histogramBin), public, dimension(numAnglesConsidered) :: &
          & angleBins = defaultAngleBins
     ! Angular bin width for spectra calculations
     real(real64), public :: deltaTheta = defaultDTheta

     ! Energy bins for spectra calculations
     type(histogramBin), public, dimension(numEnergyBins) :: &
          & energyBins = defaultEnergyBins
     ! Energy bin sub-steps for spectra calculations
     real(real64), public, dimension(numEnergyBins) :: &
          & energyBinSubStep = defaultEnergyBinSubStep

     ! Central bin values for PISA calculations
     real(real64), public, dimension(numAnglesConsidered) :: &
          & pisaAngles = defaultPisaAngles
     ! Angular bin width for PISA calculations
     real(real64), public :: pisaDTheta = defaultDTheta

     ! For ejectiles to include in spectra calculations
     logical, public, dimension(numSpecEjectiles) :: &
          & spectraEjectiles(9) = defaultEjectiles

     ! To be deprecated...
     integer(int32), public :: minEjectileRange = 1_int32
     integer(int32), public :: maxEjectileRange = 9_int32

     ! ========================================================================
     ! Regarding what's calculated - set_results
     integer(int32), public :: mspec = 0_int32   ! [0, 1, 2]
     integer(int32), public :: mpyld = 0_int32   ! [0, 1, 2]
     integer(int32), public :: mchy  = 0_int32   ! [0, 1]
     integer(int32), public :: misy  = 0_int32   ! [0, 1, 2, 3]
     integer(int32), public :: mdubl = 0_int32   ! [0, 1, 2]
     integer(int32), public :: mang  = 0_int32   ! [0, 1, 2]

     ! Momentum bins for calculation - set_KE_bins
     integer(int32), public :: numTBins = 4_int32
     real(real64),   public, dimension(:), pointer :: tMin   => NULL()
     real(real64),   public, dimension(:), pointer :: tMax   => NULL()
     real(real64),   public, dimension(:), pointer :: deltaT => NULL()


     ! Regarding angle bins - set_angles
     integer(int32), public :: numAngleBins = 0_int32
     real(real64),   public, dimension(:), pointer :: thetaMin => NULL()
     real(real64),   public, dimension(:), pointer :: thetaMax => NULL()
     real(real64),   public, dimension(:), pointer :: dTheta   => NULL()
     ! (for double differential cross sections)
     integer(int32), public :: includePisa = 0_int32   ! Include pisa spectra (/=1 excludes it)
     integer(int32), public :: numPisaAngleBins = 0_int32
     real(real64),   public, dimension(:), pointer :: aveTheta => NULL()
     real(real64),   public                        :: dThetaP  = 5.0_real64

     ! For inclusion of histograms of exciton data in output
     integer(int32), public :: includeHist = 0_int32


     ! Regarding what particles are tracked/allowed during the simulation
     !    NOTE: particle ID ranges from 1 to 9, being [n, p, d, t, He-3, He-4, pi-, pi0, pi+]
     integer(int32), public :: minParticleID = 1_int32   ! Smallest particle ID to consider 
     integer(int32), public :: maxParticleID = 9_int32   ! Largest  particle ID to consider (pi+)

     ! For fission behaviors
     ! Flag to take into account fission processes during Evaporation
     integer(int32), public :: accountFission = 1_int32
     ! Use fission events ONLY for multiplicities
     logical,        public :: fisOnly        = .false.
     
   contains
     private
     procedure, private :: validateOutputObj
  end type GSMOutput
