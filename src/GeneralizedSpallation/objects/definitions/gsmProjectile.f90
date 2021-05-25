
! ==============================================================================
!
!> \file
!> \brief   Contains the implementation of \c gsmProjectile
!> \author  CMJ, XCP-3 (LANL)
!
! ==============================================================================
!  DEFAULT VALUES
! ==============================================================================

! --------------------------------------------- PROJECTILE NAME ----------------

  !> Max allowed length of a particle name
  integer(int32), private, parameter :: maxPartNameLen   = 6

  !> The total number of projectile names defined
  integer(int32), public,  parameter :: numDefaultProjNames = 7

  !> Defines the available projectile names for a simulation
  character(len=4), public, parameter, dimension(numDefaultProjNames) :: &
       & projNames = [ "prot", "neut", "pipl", "pimi", "pize", "gamm", "gamb" ]

  !> Default projectile name
  character(len=maxPartNameLen), private, parameter :: defaultProjName = projNames(1)


! --------------------------------------------- PROJECTILE COMPOSITION ---------

  !> Default projectile baryon number (A number)
  integer(int32), private, parameter :: defaultProjBaryons  =  1_int32

  !> Default projectile proton number (Z number; charge) 
  integer(int32), private, parameter :: defaultProjProtons  =  1_int32

  !> Default projectile quantum decay number (not presently used)
  integer(int32), private, parameter :: defaultDecayNumber  =  0_int32

  !> Default projectile particle flag (nucleus, pion, or photon)
  integer(int32), private, parameter :: defaultProjPartFlag = nucleusProjFlag


! --------------------------------------------- PROJECTILE ENERGY ----------------

  !> Default projectile rest mass [GeV/c**2] (negative flags GSM to calculate) 
  real(real64),   private, parameter :: defaultRestMass     = -1.0_real64

  !> Default projectile kinetic energy [GeV] 
  real(real64),   private, parameter :: defaultKinEnergy    =  1.0_real64

  !> Default projectile maximum kinetic energy [GeV]
  real(real64),   private, parameter :: defaultKinEnergyMax = 1.0E+09

  !> Default projectile incident energy increment [MeV]
  real(real64),   private, parameter :: defaultDKinEnergy   = -50.0_real64


! --------------------------------------------- PROJECTILE MULTIPLIERS ---------

  !> Default projectile \f$ a_{f} \f$ multiplier
  real(real64),   private, parameter :: defaultAfMultiplier = -1.0_real64

  !> Default projectile \f$ C\left(Z\right) \f$ multiplier
  real(real64),   private, parameter :: defaultCzMultiplier = -1.0_real64

  
! --------------------------------------------- PROJECTILE MISC. DATA ----------

  !> Default projectile simulation system (lab or antilab.)
  integer(int32), private, parameter :: defaultSystem       = labSystem


! ==============================================================================
!
!  DATA TYPE DESCRIPTION:
!
!> \class  gsmProjectile
!> \brief  Defines the projectile object utilized by clients for GSM simulations
!
!> Defines the projectile object utilized by clients for simulations of GSM. The
!> \c gsmProjectile object defines characteristics of the incident particle during
!> a high energy nuclear event.
!
!> Supported incident particles (projectiles) include neutrons, protons, singly
!> or neutrally-charged pions, photons (mono. energetic and bremsstrahlung), and
!> light- and heavy-ions.
!
!> The data type is passed throughout the simulation object, being stored only by
!> local variables. The object is only pointed to by these variables, thus
!> minimizing any data replication where possible. Note the provided \c
!> gsmProjectile object may be altered by GSM during the simulation depending on
!> the specification of the object by the software client.
!
! ==============================================================================
type, public :: gsmProjectile
     private

     !> The name of the particle
     character(LEN=maxPartNameLen), public :: particleName = defaultProjName

     !> The number of nucleons (baryons) contained within the particle
     integer(int32), public :: numBaryons   = defaultProjBaryons

     !> The number of protons (or charge) of the particle
     integer(int32), public :: numProtons   = defaultProjProtons

     !> The particle's rest energy [GeV/c**2]
     real(real64),   public :: restMass     = defaultRestMass

     !> The particle's incident (kinetic) energy [GeV]
     real(real64),   public :: kinEnergy    = defaultKinEnergy

     !> The particle's maximum allowed incident (kinetic) energy [GeV]
     real(real64),   public :: kinEnergyMax = defaultKinEnergyMax

     !> The incremental step size of the particle's incident energy [MeV]
     !> Note: Will not be utilized if the value is (<= 0).
     real(real64),   public :: dKinEnergy   = defaultDKinEnergy

     !> Quantum decay number of the particle
     integer(int32), public :: decayNumber  = defaultDecayNumber

     !> Flags what type of particle the projectile is (nucleus, pion, or photon
     !> [mono. and brems. available])
     integer(int32), public :: particleFlag = defaultProjPartFlag

     !> Data for a bremstrahlung photon (if so)
     type(BremsPhoton), public, pointer :: brems => NULL()

     !>  Particle's \f$ a_{f} \f$ fission multiplier
     real(real64),   public :: afMultiplier = defaultAfMultiplier

     !> The particle's \f$ C\left(Z\right) \f$ fission multiplier
     real(real64),   public :: czMultiplier = defaultCzMultiplier

     !> Flags what system the simulation occurs in (lab. or antilab. available)
     integer(int32), private :: system       = defaultSystem

   contains
     private

     !> @{
     !> Describes the particle to a provided I/O stream
     procedure, private :: describeProjectile
     procedure, public  :: describe => describeProjectile
     !> @}

  end type GSMProjectile

