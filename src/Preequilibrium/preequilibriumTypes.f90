

! ============================================================================================
!
! File contains the definition of various derived/data types used by the preequilibrium model
!
! ============================================================================================


  ! For all message handling:
  type, private :: preequilibriumIO
     private
     character(LEN=512),   private                  :: message =  ""
     procedure(IOHANDLER), private, pointer, nopass :: print   => printPreeq
  end type preequilibriumIO


  ! Contains options available to the preequilibrium simulation (ALL user specifiable)
  type, public  :: preequilibriumOptions
     private
     real(real64),   public  :: r0Mult         = defaultR0Multiplier    != Radius multipler [fm]
     integer(int32), public  :: numPreeqType   = defaultNumPreeqType    != Number of preequilibrium particles to consider   [Unitless]
     integer(int32), public  :: levelDenParam  = defaultLevelDenParam   != Flag for which parameters to use for level density calculation
     real(real64),   public  :: emissionWidth  = defaultEmissionWidth   != Standard deviation preeq. emission probability (follows gaussian curve)
     integer(int32), public  :: excludePreeq   = defaultExcludePreeq    != Flags to NOT simulate preequilibrium emission ( /= 0 is 'on')
  end type preequilibriumOptions



  ! For the excited residual nucleus
  type, public :: residualNucleus
     private
     ! Particle A, Z
     real(real64),   public :: numBaryons  = 0.0_real64   != Mass number (Z+N)        [Unitless]
     real(real64),   public :: numProtons  = 0.0_real64   != Atomic Number (Charge)   [Unitless]
     ! Particle Energies
     real(real64),   public :: kinEnergy   = 0.0_real64   != Total Kinetic Energy     [MeV]
     real(real64),   public :: thermEnergy = 0.0_real64   != Thermal Energy           [MeV]
     real(real64),   public :: recEnergy   = 0.0_real64   != Recoil Energy            [MeV]
     real(real64),   public :: rotEnergy   = 0.0_real64   != Rotational Energy        [MeV]
     real(real64),   public :: restMass    = 0.0_real64   != Rest Mass                [GeV/c**2]
     ! (fission barrier)
     real(real64),   public :: fissBarr    = 0.0_real64   != Fission Barrier Height   [MeV]
     ! Linear Momentum
     real(real64),   public :: linearMomX  = 0.0_real64   != X-momentum               [GeV/c]
     real(real64),   public :: linearMomY  = 0.0_real64   != Y-momentum               [GeV/c]
     real(real64),   public :: linearMomZ  = 0.0_real64   != Z-momentum               [GeV/c]
     ! Angular Momentum
     integer(int32), public :: angMomFlag  = 0_int32      != Angular Momentum         [Quantum Number; Unitless]
     real(real64),   public, dimension(3) :: angularMom = 0.0_real64   != Angular momentum [GeV * (time)]
     real(real64),   public, dimension(3) :: normSpeed  = 0.0_real64   != Normalized momentum OR speed [Unitless or ?]
     ! Physics flag
     integer(int32), public :: state = residualState      != Flag for state of compound
                                  !  =0 for residual (undergo preequilibrium decay)
                                  !  =1 for compound (undergo evaporation)
                                  !  =2 for stable nucleus
                                  ! <=0 non-existent (underwent Fermi Break-up)
  end type residualNucleus


  type, public :: preeqExcitonData
     private
     integer(int32), public :: numTotal    = 0_int32
     real(real64),   public :: numNucleons = 0_int32
     real(real64),   public :: numProtons  = 0_int32
     real(real64),   public :: numHoles    = 0_int32
  end type preeqExcitonData



  ! For storage of fragments
  type, public :: preequilibriumFragment
     private
     real(real64),   public :: numBaryons  = 0.0_real64   != Mass number (Z+N)     [Unitless]
     real(real64),   public :: numProtons  = 0.0_real64   != Atomic Number (Z)     [Unitless]
     real(real64),   public :: kinEnergy   = 0.0_real64   != Kinetic Energy        [MeV]
     real(real64),   public :: restMass    = 0.0_real64   != Rest Mass of fragment [MeV/c**2]
     real(real64),   public :: theta       = 0.0_real64   != Theta of fragment     [Degrees]
     real(real64),   public :: phi         = 0.0_real64   != Phi   of fragment     [Degrees]
     integer(int32), public :: origin      = preequilibriumProgeny   != Fragment Origin [Unitless];
                                  ! =0 from preequilibrium decay;
                                  ! =1 from Fermi Break-up
  end type preequilibriumFragment



  ! For containg all of the results:
  type, public :: preequilibriumResults
     private

     ! Flags to the object that progeny array is NOT associated with any location in memory
     logical, private :: constructed = .FALSE.


     ! Progeny information:
     integer(int32), public  :: numProgeny = 0_int32
     integer(int32), public  :: maxProgeny = 0_int32
     type(preequilibriumFragment), public, dimension(:), pointer :: progenyBnk => NULL()


     ! Nucleus information:
     type(residualNucleus), public :: initResidual
     type(residualNucleus), public :: residual


     ! Indicator of how the simulation ended:
     integer(int32), public :: simState = 0_int32

  end type preequilibriumResults



  ! Contains all arrays used for the preequilibrium simulation
  type, private :: preequilibriumCalculation
     private
     ! NOTE: "Emission channels" are associated with a preequilbrium fragment type (i.e. channel j corresponds to emission of particle j)
     ! (Temporary values for the residual)
     integer(int32), private :: in      =  0_int32      ! Number of neutrons in   the residual [Unitless]
     integer(int32), private :: iz      =  0_int32      ! Number of protons in    the residual [Unitless]
     real(real64),   private :: athrd   =  0.0_real64   ! A**(1/3) for            the residual [Unitless]

     ! (Information regarding selected emission channel?)
     real(real64),   private :: ac      =  0.0_real64   ! Auxiliary quantity for preequilibrium emission width  [Unitless]
     real(real64),   private :: exn     =  0.0_real64   ! Exciton nucleons + holes                              [Unitless]
     real(real64),   private :: dl      =  0.0_real64   ! ?Mass Excess for                    current emission? [Unitless]

     ! (Information regarding emission channels or characteristics of the fragment corresponding to the emission channel)
     real(real64),   private, dimension(numCalcElements)     :: afj     =  0.0_real64   ! A (mass number) for                        fragment j [Unitless]
     real(real64),   private, dimension(numCalcElements)     :: zfj     =  0.0_real64   ! Z (atomic number) for                      fragment j [Unitless]
     real(real64),   private, dimension(numCalcElements)     :: vj      =  0.0_real64   ! V (coulomb barrier) for                    fragment j [MeV]
     real(real64),   private, dimension(numCalcElements)     :: rj      =  0.0_real64   ! Free thermal energy at the Coulomb barrier for fragment j [MeV]
     real(real64),   private, dimension(numCalcElements)     :: uej     =  0.0_real64   ! Free energy for                    emission channel j [MeV]
     real(real64),   private, dimension(numCalcElements)     :: pevapj  =  0.0_real64   ! Pairing energy for                 emission channel j [MeV]
     real(real64),   private, dimension(numCalcElements)     :: afjthr  =  0.0_real64   ! A**(1/3) for the particle in       emission channel j [Unitless]
     real(real64),   private, dimension(numAllowedFragments) :: alj     =  0.0_real64   ! Exciton phase space factors for    emission channel j [?]
     real(real64),   private, dimension(numAllowedFragments) :: ami     =  0.0_real64   ! Level density parameter for        emission channel j [Unitless]
     real(real64),   private, dimension(numCalcElements)     :: bj      =  0.0_real64   ! Binding energy of the particle for emission channel j [MeV]
     real(real64),   private, dimension(numAllowedFragments) :: gb      =  0.0_real64   ! Gamma-Beta multiplier for          emission channel j [Unitless?]
     real(real64),   private, dimension(numAllowedFragments) :: redpre  =  0.0_real64   ! Reduced mass for                   emission channel j [MeV/c**2]
  end type preequilibriumCalculation
