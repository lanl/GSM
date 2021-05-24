
! ======================================================================
!
! Contains all of the default values regarding output file generation
!
! ======================================================================

  ! For naming
  integer(int32), parameter, private :: nameLength = 128
  character(*), parameter, private :: &
       & defaultInpFile = "gsm.inp", & ! Default input file name
       & defaultOutFile = "gsm.out", & ! Default output file name (contains simulation results)
       & defaultAuxFile = "gsm.aux"    ! Default aux. file name (contains simulation warnings and errors)

  ! Maximum allowed input file comments
  integer(int32), parameter, private :: maxCommentLines = 10_int32

  ! Default number of inelastic events to simulate
  integer(int64), parameter, private :: defaultInelasticEvents = 100000_int64

  ! Output options
  integer(int32), parameter, private :: &
       & defaultAngSpec       = 0_int32, & ! Default for angular integrated spectra
       & defaultEnergySpec    = 0_int32, & ! Default for energy integrated spectra
       & defaultDoubDiffSpec  = 0_int32, & ! Default for double differential spectra
       & defaultMultiplicity  = 0_int32, & ! Default for calculating particle multiplicities
       & defaultChannelXSCalc = 0_int32, & ! Default for calculating reaction channel cross sections
       & defaultNuclideXSCalc = 0_int32    ! Default for calculation nuclide cross section
  logical, parameter, private :: &
       & defaultPISAPrint     = .FALSE., & ! Default for printing PISA cross section tables
       & defaultHistPrint     = .FALSE.    ! Default for printing HIST residual tables

  ! Default angular bin width for spectra calculation
  !    (energy int., angle int., double diff.)
  real(real64), parameter, private :: &
       & defaultDTheta = 5.0_real64, &
       & halfDTheta = defaultDTheta / 2.0_real64

  ! The number of angles considered for particle emission
  integer(int32), parameter, public :: numAnglesConsidered = 10_int32

  ! The default angles considered for PISA histograms
  real(real64), parameter, private, dimension(numAnglesConsidered) :: &
       & defaultPisaAngles = [ &
       &  15.0_real64,  30.0_real64,  45.0_real64,  60.0_real64, &
       &  75.0_real64,  90.0_real64, 105.0_real64, 120.0_real64, &
       & 135.0_real64, 150.0_real64 &
       & ]

  ! The default number of angles considered for ejectile spectra
  type(histogramBin), parameter, private, dimension(numAnglesConsidered) :: &
       & defaultAngleBins = [ (histogramBin( &
       &       defaultPisaAngles(iter) - halfDTheta, defaultPisaAngles(iter) + halfDTheta), &
       &    iter=1,numAnglesConsidered) ]

  ! Number of grouped energy bins considered
  integer(int32), parameter, private :: numEnergyBins = 4_int32

  ! Default grouped energy bins
  type(histogramBin), parameter, private, dimension(numEnergyBins) :: &
       & defaultEnergyBins = &
       & [ histogramBin(   0.0_real64,  250.0_real64), &
       &   histogramBin( 250.0_real64,  500.0_real64), &
       &   histogramBin( 500.0_real64, 1000.0_real64), &
       &   histogramBin(1000.0_real64, 2500.0_real64)  &
       & ]
  ! Default sub-step size for energy bins
  real(real64), parameter, private, dimension(numEnergyBins) :: &
       & defaultEnergyBinSubStep = &
       & [  2.0_real64, &
       &    5.0_real64, &
       &   10.0_real64, &
       &   25.0_real64  &
       & ]
  
  ! Number of ejectiles considered for spectra calculations
  integer(int32), parameter, private :: numSpecEjectiles = 9_int32

  ! Default flags for ejectile includsion regarding spectra calculations
  logical, parameter, private, dimension(numSpecEjectiles) :: &
       & defaultEjectiles = [(.TRUE., iter=1, 9)]
