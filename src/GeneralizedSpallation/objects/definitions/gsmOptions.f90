
! ==============================================================================
!
! This file contains all options specifications for each of the submodels and for GSM
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

  ! Contains all simulation options used by GSM (for itself and its submodels)
  type, public :: GSMOptions
     private

     ! Base Physics models:
     type(molnixOptions),          public :: molnixOpts
     type(fissionBarrierOptions),  public :: fissBarOpts

     ! Sub-model options:
     type(sDCMDataOptions),        public  :: sDCMDataOpts
     type(sDCMOptions),            public  :: sDCMOpts
!     type(mDCMDataOptions),        public  :: mDCMDataOpts
!     type(mDCMOptions),            public  :: mDCMOpts
     type(coalescenceDataOptions), public  :: coalesDataOpts
     type(coalescenceOptions),     public  :: coalesOpts
     type(fermiBreakUpOptions),    public  :: fbuOpts
     type(preequilibriumOptions),  public  :: preeqOpts
     type(evaporationDataOptions), public  :: evapDataOpts
     type(evaporationOptions),     public  :: evapOpts

     ! For some physics handling:
     logical, public :: conserveTotalMomentum   = defaultConserveTotalMomentum
     logical, public :: tallyResidualProjectile = defaultTallyResidualProjectile
     logical, public :: useSDCMAfCzMultipliers  = defaultSDCMAfCzMultipliers

     ! For general client interactions:
     integer(int64), public  :: printIncrement = defaultPrintIncrement


     ! Transition energy information
     logical,        public  :: smoothTransition      = smoothINCTransitionDefault
     real(real64),   public  :: transitionWidth       = transitionWidthDefault
     real(real64),   public  :: nucleonTransitionE    = nucleonTransEDefault
     real(real64),   public  :: pionTransitionE       = pionTransEDefault
     real(real64),   public  :: monoPhotonTransitionE = monoPhotonTransEDefault
     real(real64),   public  :: bremPhotonTransitionE = bremPhotonTransEDefault

  end type GSMOptions


  ! Contains all information needed for the I/O handler
  ! NOTE: CHANGE ALL TO PRIVATE AFTER COMPLETING THE GSM CLASS!
  type, private :: GSMIO
     private
     procedure(IOHANDLER), public,  nopass, pointer :: print   => printGSM
     character(LEN=512),   public                   :: message =  ""
  end type GSMIO
