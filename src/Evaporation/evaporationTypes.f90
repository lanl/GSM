
! ===================================================================================
!
! This file contains the types used by the evaporation class.
!
! ===================================================================================

  ! Data structure for evaporation options
  type, public :: evaporationOptions
     private
     integer(int32), public  :: numEvapType      = defaultNumEvapType
     real(real64),   public  :: inverseParameter = defaultInverseParameter
     integer(int32), public  :: fissParameter    = defaultFissParameter
     integer(int32), public  :: redCompTime      = defaultRedCompTime
!     type(fermiBreakupOptions), public :: fermiOptions
  end type evaporationOptions


  ! For handling messages
  type, private :: evaporationIO
     private
     character(LEN=512),   private                  :: message =  ""
     procedure(IOHANDLER), private, pointer, nopass :: print   => printEvaporation
  end type evaporationIO



  ! Data structure for emitted fragment list
  type, public :: evaporationFragment
     private
     real(real64),   public  :: numBaryons  = 0.0_real64   != Baryon number         of fragment [Unitless]
     real(real64),   public  :: numProtons  = 0.0_real64   != Charge                of fragment [Unitless]
     real(real64),   public  :: kinEnergy   = 0.0_real64   != Kinetic Energy        of fragment [MeV]
     real(real64),   public  :: restMass    = 0.0_real64   != Rest mass             of fragment [GeV/c**2]
     real(real64),   public  :: linearMomX  = 0.0_real64   != Normalized X-momentum of fragment [Unitless]
     real(real64),   public  :: linearMomY  = 0.0_real64   != Normalized Y-momentum of fragment [Unitless]
     real(real64),   public  :: linearMomZ  = 0.0_real64   != Normalized Z-momentum of fragment [Unitless]
     integer(int32), public  :: physicsFlag = evaporationFlag   != Flags fragment origen
                                ! =0 for Evaporation;
                                ! =1 for Fermi Breakup;
                                ! =2 for Evaporation from Fission Fragment
  end type evaporationFragment


  ! Data structure for residual particle
  type, public :: residual
     private
     real(real64),   public  :: numBaryons   = 0.0_real64   != Baryon number   of residual [Unitless]
     real(real64),   public  :: numProtons   = 0.0_real64   != Charge          of residual [Unitless]
     real(real64),   public  :: kinEnergy    = 0.0_real64   != Kinetic energy  of residual [MeV]
     real(real64),   public  :: recEnergy    = 0.0_real64   != Recoil energy   of residual [MeV]
     real(real64),   public  :: linearMomTot = 0.0_real64   != Linear momentum of residual [MeV/c]
     integer(int32), public  :: physicsFlag  = evaporationFlag   != Flags the processes underwent by the residual
                                ! =0 for Evaporation;
                                ! =1 for Evaporation, Fermi Breakup, or both;
                                ! =2 for Evaporation, Fission, or both
  end type residual


  ! Data structure for a fission fragment
  type, public :: fissionFragment
     private
     real(real64),   public  :: numBaryons = 0.0_real64   != Baryon number       of fission fragment [Unitless]
     real(real64),   public  :: numProtons = 0.0_real64   != Charge              of fission fragment [Unitless]
     real(real64),   public  :: kinEnergy  = 0.0_real64   != Kinetic energy      of fission fragment [MeV]
     real(real64),   public  :: rotEnergy  = 0.0_real64   != Rotational energy   of fission fragment [MeV]
     real(real64),   public, dimension(3) :: linearMomFrac  = 0.0_real64    != Normalized momentum of fission fragment [Unitless]
     integer(int32), public  :: physicsFlag  = fissionFragFlag   != Flags the processes underwent by the fission fragment
                                ! =-1 for not being created (i.e. it doesn't exist)
                                ! = 0 for Evaporation;
                                ! = 1 for Evaporation, Fermi Breakup, or both;
                                ! = 3 for Evaporation, and SHOULD undergo fission [the model cannot accomodate this currently]
                                ! NOTE: =3 could be changed to mean "Evaporation, Fission, or both" if the algorithm is updated to include this feature
  end type fissionFragment


  ! Data structure for a the compound (for starting the simulation)
  type, public :: evapCompound
     private
     real(real64),   public  :: numBaryons   = 0.0_real64   != Baryon number       of compound [Unitless]
     real(real64),   public  :: numProtons   = 0.0_real64   != Charge              of compound [Unitless]
     real(real64),   public  :: kinEnergy    = 0.0_real64   != Kinetic energy      of compound [MeV]
     real(real64),   public  :: recEnergy    = 0.0_real64   != Recoil              of compound [MeV]
     real(real64),   public  :: linearMomTot = 0.0_real64   ! Total linear momentum [MeV/c]
     real(real64),   public, dimension(3) :: linearMomFrac  = 0.0_real64    != Normalized momentum of compound [Unitless]
     ! For multipliers of fission parameters:
     real(real64),   public  :: afMultiplier = defaultAfMultipler ! Multiplier for A_f fission parameter
     real(real64),   public  :: czMultiplier = defaultCzMultipler ! Multiplier for C(Z) fission parameter
  end type evapCompound


  type, public :: evaporationResults
     private
     ! Regarding progeny:
     integer(int32), public  :: numProgeny        = 0_int32   ! Number of progeny existing in the bank
     integer(int32), private :: progenyBnkSize    = 0_int32   ! Maximum bank size
     integer(int32), private :: progenyBnkSizeM1  = 0_int32   ! Maximum bank size, minus 1
     type(evaporationFragment), public, dimension(:), pointer :: progenyBnk => NULL()   ! Array of evaporation fragments

     ! Regarding the residual:
     type(residual), public  :: compound

     ! Regarding fission:
     real(real64),   public  :: fissionProbability  = 0.0_real64   ! Fission probability
     integer(int32), public  :: numFissionFragments = 0_int32      ! Number of fission fragments existing
     integer(int32), private :: fissionBnkSize      = 0_int32      ! Size of the fission fragment array
     type(fissionFragment), public, dimension(:), pointer :: fissionBnk => NULL()   ! Array of fission fragments

     ! Regarding simulation state
     integer(int32), public  :: simState = 0_int32

     ! Misc.
     logical,        private :: constructed = .FALSE.   ! Flags to evaporation model if results were established
  end type evaporationResults


  ! For calculation values that aren't used by the client:
  type, private :: evaporationCalculation
     private
     ! Internal calculation array(s)
     ! NOTE: These variables used to be embedded into common blocks
     real(real64), private                     :: smalla0 = 0.0_real64
     real(real64), private, dimension(maxSize) :: r       = 0.0_real64   != Decay Width
     real(real64), private, dimension(maxSize) :: rr      = 1.0_real64   != Decay Width enhancement facor by excited state particle emission
     real(real64), private                     :: sigma   = 0.0_real64   != Toal decay width
     real(real64), private, dimension(maxSize) :: gj      = 0.0_real64   != Spin Multiplicity Factor times mu (see gj in Eqn 39 of CEM03 Manual)
     real(real64), private, dimension(maxSize) :: q       = 0.0_real64   != Q-value
     real(real64), private, dimension(maxSize) :: v       = 0.0_real64   != Coulomb Barrier
     real(real64), private, dimension(maxSize) :: delta   = 0.0_real64   != Pairing Energy
     real(real64), private, dimension(maxSize) :: smalla  = 0.0_real64   != Level density parameter
     ! Used only in "gamma" routine:
     real(real64), private, dimension(maxSize) :: couk    = 1.0_real64
     real(real64), private, dimension(maxSize) :: couc    = 0.0_real64
  end type evaporationCalculation
