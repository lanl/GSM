
! ==============================================================================
!
! This file contains all of the data types used by the standard Dubna Cascade
! Model's Data Class
!
!
! Written by CMJ, XCP-3, 12/2018
!
! ==============================================================================

  ! Object used to set the behavior of the sDCM simulation
  type, public :: sDCMOptions
     private
     ! Includes or excludes reflection/refraction of particles at nuclear zone boundaries (default = 0 [off])
     integer(int32), public :: boundaryEffects = boundaryEffectsDefault

     ! Use new angular distribution method (in costa) where applicable
     integer(int32), public :: newAngularDistributions = newAngularDistributionsDefault
  end type sDCMOptions


  ! Object used to contain all I/O options for GSm
  type, private :: sDCMIO
     private
     procedure(IOHANDLER), private, nopass, pointer :: print   => printSDCM
     character(LEN=512),   private                  :: message =  ""
  end type sDCMIO



  ! Progeny tracking (information tracked by the S-DCM during simulation)
  type, public :: sDCMProgeny
     private
     real(real64),   public  :: xCoord      = 0.0_real64   ! X location of particle [fm]
     real(real64),   public  :: yCoord      = 0.0_real64   ! Y location of particle [fm]
     real(real64),   public  :: zCoord      = 0.0_real64   ! Z location of particle [fm]
     real(real64),   public  :: sinTheta    = 0.0_real64   ! sin(theta) of direction of motion
     real(real64),   public  :: cosTheta    = 0.0_real64   ! cos(theta) of direction of motion
     real(real64),   public  :: sinPhi      = 0.0_real64   ! sin(phi) of direction of motion
     real(real64),   public  :: cosPhi      = 0.0_real64   ! cos(phi) of direction of motion
     real(real64),   public  :: kinEnergy   = 0.0_real64   ! Kinetic energy of the particle [GeV]
     real(real64),   public  :: restMass    = 0.0_real64   ! Rest energy of the particle [GeV/c**2]
     integer(int32), public  :: numProtons  = 0_int32      ! Number of protons  in the nucleus
     integer(int32), public  :: numBaryons  = 0_int32      ! Number of nucleons in the nucleus
     integer(int32), public  :: strangeness = 0_int32      ! Quantum decay number of particle
     integer(int32), public  :: photonFlag  = 0_int32      ! Flags if particle is a photon (>0) or not (=0)
     integer(int32), public  :: nuclearZone = 0_int32      ! Zone number where particle is located
     integer(int32), public  :: index       = 0_int32      ! Index of particle's creation
  end type sDCMProgeny



  ! For the resulting compound nucleus
  type, public :: sDCMCompound
     private
     real(real64), public :: numBaryons  = 0.0_real64   ! Number of nucleons in the nucleus
     real(real64), public :: numProtons  = 0.0_real64   ! Number of protons  in the nucleus
     real(real64), public :: kinEnergy   = 0.0_real64   ! Kinetic energy of the particle [GeV]
     real(real64), public, dimension(3) :: linearMom  = 0.0_real64   ! Linear  momentum of the particle [GeV/c]
     real(real64), public, dimension(3) :: angularMom = 0.0_real64   ! Angular momentum of the particle [GeV/c]
  end type sDCMCompound



  ! Create a "projectile" object
  type, public :: sDCMProjectile
     private

     ! Projectile nucleus composition
     integer(int32), public :: numBaryons   = 0_int32     ! Number of baryon (i.e. nucleons) in the nucleus
     integer(int32), public :: numProtons   = 0_int32     ! Number of protons in the nucleus

     ! Projectile energy information
     real(real64),   public :: kinEnergy    = 0.0_real64  ! Kinetic energy of projectile (in target's rest frame) [GeV]
     real(real64),   public :: restMass     = 0.0_real64  ! Rest energy of projectile [GeV]

     ! Projectile quantum number and flag
     integer(int32), public :: decayNumber  = 0_int32     ! Decay quantum number of projectile
     integer(int32), public :: gammaFlag    = 0_int32     ! Flags if projectile is photon (>0) or not (=0)
  end type sDCMProjectile



  ! Contains all data regarding exciton information
  type, private :: sDCMExcitons
     private
     integer(int32), public :: numExcProt    = 0_int32   ! Number of exciton charges
     integer(int32), public :: numExcNeut    = 0_int32   ! Number of exciton neutrons
     integer(int32), public :: numExcHoles   = 0_int32   ! Number of holes
  end type sDCMExcitons



  ! Create a "results" object (passed in to class for simulation)
  type, public :: StandardDCMResults
     private
     ! For progeny...
     integer(int32),     public  :: numProgeny   = 0_int32   ! Tracks number of progeny created during sDCM simulation
     integer(int32),     private :: maxProgeny   = 0_int32   ! Stores max allowed number of progeny
     integer(int32),     private :: maxProgenyM3 = 0_int32   ! Stores max allowed number of progeny, minus 3
     type(sDCMProgeny),  public,  pointer, dimension(:) :: progenyBnk  => NULL()

     ! These arrays contain the bank of cascade particles that are still interacting
     real(real64),       private, pointer, dimension(:,:) :: pmemo => NULL()
     integer(int32),     private, pointer, dimension(:,:) :: imemo => NULL()

     ! For compound nucleus
     type(sDCMCompound), public  :: residual

     ! For the created excitons
     type(sDCMExcitons), public  :: excitons

     ! Misc.
     integer(int32),     public  :: numElastic   = 0_int32   ! Number of elastic interactions that occurred
     logical,            private :: initialized  = .FALSE.   ! Flags that object has NOT been initialized (for particle bank)

     ! Flags regarding the simulation
     integer(int32),     public  :: simState     = 0_int32   ! Flags state of the sDCM simulation
     integer(int32),     public, dimension(2) :: &
          & numRestarts = 0_int32         ! Flags how many times the INC was restarted due to either
                                          ! (1) an issue in GEOM8 or
                                          ! (2) if 'typint' couldn't determine the number of secondaries
  end type StandardDCMResults


  ! Create a "calculation variables" object
  type, private :: sDCMPhotonCrossSections
     private
     ! For photons (fine angular photon cross sections [si] and CDF at each angle [ri])
     ! (From old /isecgpn/ common block)
     real(real64),   private, dimension(2, 182) :: si = 0.0_real64
     real(real64),   private, dimension(2, 181) :: ri = 0.0_real64
  end type sDCMPhotonCrossSections


  type, private :: sDCMPauliInfo
     private

     ! Flags for N N -> pi N N reactions
     integer(int32), private :: indi = 1_int32

     ! Generation information for a particle produced during the SDCM simulation
     ! (from old /gener/ common block)
     integer(int32), private :: ing = 0_int32
     integer(int32), private :: meso = 0_int32
     integer(int32), private, dimension(300) :: ngen = 0_int32
     integer(int32), private, dimension(300) :: igs  = 0_int32
  end type sDCMPauliInfo
