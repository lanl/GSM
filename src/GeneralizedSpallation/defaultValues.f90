
! ======================================================================
!
! This file contains all default values utilized by GSM and its sub-objects
!
! Written by CMJ, XCP-3, August 2018 - Decemeber 2018
!
! ======================================================================

  ! For printing:
  integer(int32), private, parameter :: defaultGSMVerbose = 4_int32


  ! For projectile flagging:
  integer(int32),  public, parameter :: nucleusProjFlag  = 0_int32   ! Flags projectile as a nucleus
  integer(int32),  public, parameter :: photonProjFlag   = 1_int32   ! Flags projectile as a mono energetic photon
  integer(int32),  public, parameter :: bremsProjFlag    = 2_int32   ! Flags projectile as a brems. photon
  integer(int32),  public, parameter :: pionProjFlag     = 3_int32   ! Flags projectile as a pion
  ! For stating system calculation is done in:
  integer(int32), private, parameter :: labSystem       = 0_int32   ! Flags use of the lab. system
  integer(int32), private, parameter :: antilabSystem   = 1_int32   ! Flags use of the anti-lab. system
  ! For flagging INC usage:
  integer(int32), public,  parameter :: sDCMFlagged     = 0_int32   ! Flags to use the sDCM INC model
  integer(int32), public,  parameter :: mDCMFlagged     = 1_int32   ! Flags to use the mDCM INC model
  integer(int32), public,  parameter :: defaultINCFlag  = sDCMFlagged


  ! -----------------------------------------------------------------------------
  ! For scaling the combined baryon number to estimate maximum progeny bank size:
  real(real64),   private, parameter :: defaultBankScaling = 1.30_real64
  real(real64),   public :: bankScaling = defaultBankScaling
  integer(int32), private, parameter :: defaultMinBankSize = 50_int32
  integer(int32), public :: minBankSize = defaultMinBankSize

  ! -----------------------------------------------------------------------------
  ! Regarding the end-simulation state (in the results object):
  integer(int32),  public, parameter :: successfulSingleEvent =   0_int32
  integer(int32),  public, parameter :: failedSingleEvent     =   5_int32   
  integer(int32),  public, parameter :: noDataConstruction    =  50_int32
  integer(int32),  public, parameter :: noResultsConstruction =  60_int32
  integer(int32),  public, parameter :: noObjectConstruction  =  70_int32


  ! -----------------------------------------------------------------------------
  ! Used for warning users when outside applicable energy regimes for GSM
  real(real64),   private, parameter :: minRecommendedEnergy = 100.0_real64   ! Minimum energy GSM is recommended for simulation [MeV]
  real(real64),   private, parameter :: maxRecommendedEnergy =   1.0_real64   ! Maximum energy GSM is recommended for simulation [TeV/A]


  ! -----------------------------------------------------------------------------
  ! Default values for GSM Options
  ! (regarding physics handling):
  logical,        private, parameter :: defaultConserveTotalMomentum   = .FALSE.   ! Conserve momentum on total or on average after INC
  logical,        private, parameter :: defaultTallyResidualProjectile = .TRUE.   ! Tally residual projectile nuclei
  logical,        private, parameter :: defaultSDCMAfCzMultipliers     = .TRUE.   ! Use the sDCM multipliers (else outdated mDCM multipliers)
  ! (misc.):
  integer(int32), private, parameter :: defaultPrintIncrement = -1   ! <=0 GSM determines; >0 use provided value
  ! (Regarding the transition between sDCM and mDCM physics):
  logical,        private, parameter :: smoothINCTransitionDefault = .TRUE.     ! Flags use of a smooth transition INC energy (not hard values)
  real(real64),   private, parameter :: transitionWidthDefault     =  100.0_real64   ! Default width to sample INC usage for
  real(real64),   public,  parameter :: nucleonTransEDefault       = 4490.0_real64   ! Default INC transition energy for incident nucleons [MeV]
  real(real64),   public,  parameter :: pionTransEDefault          = 2490.0_real64   ! Default INC transition energy for incident pions [MeV]
  real(real64),   public,  parameter :: monoPhotonTransEDefault    = 1200.0_real64   ! Default INC transition energy for incident monoenergetic photons [MeV]
  real(real64),   public,  parameter :: bremPhotonTransEDefault    = 5000.0_real64   ! Default INC transition energy for incident bremstrahhlung photons [MeV]

  ! For internal error protection:
  real(real64),   private, parameter :: div0Lim = 1.0d-15


  ! -----------------------------------------------------------------------------
  ! Misc. Values:
  character(len=*), public, parameter :: stopName = "stop"
