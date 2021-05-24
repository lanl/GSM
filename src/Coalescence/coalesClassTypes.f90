
! ==============================================================================
!
! This file contains the various data types used by the Coalescence class. These
! are:
!
! (public)
!    1. Options type
!    2. Data type
!    3. Progeny type
!    4. Results type (includes progeny type)
! (private)
!    5. I/O
!
! 
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

  type, public :: coalescenceOptions
     private
     integer(int32), public  :: expandedCoalescence = defaultExpandedCoalescence
  end type coalescenceOptions


  type, public :: coalescenceDataOptions
     private
     ! NOTE: Values < 0 will use the default specified radius!
     !       Clients need only alter the desired value!
     real(real64), public  :: coalesRadiiDeut  = -1.0_real64   ! Coalescence radius for deuterons       (D)
     real(real64), public  :: coalesRadiiTrit  = -1.0_real64   ! Coalescence radius for tritium/Helion  (T/He-3)
     real(real64), public  :: coalesRadiiAlpha = -1.0_real64   ! Coalescence radius for helium          (He-4)
     real(real64), public  :: coalesRadiiLFrag = -1.0_real64   ! Coalescence radius for light fragments (He-6, Li-6/7, Be-7)
  end type coalescenceDataOptions



  type, private :: coalescenceIO
     private
     character(LEN=512),   private                  :: message =  ""
     procedure(IOHANDLER), private, nopass, pointer :: print   => printCoales
  end type coalescenceIO


  type, private :: effectiveCoalescenceRadii
     private
     ! The effective coalescence radii:
     real(real64), private :: coalesRadiiDeut  = defaultCRDeut(1)
     real(real64), private :: coalesRadiiTrit  = defaultCRTrit(1)
     real(real64), private :: coalesRadiiAlpha = defaultCRAlpha(1)
     real(real64), private :: coalesRadiiLFrag = defaultCRLFrag(1)

     ! The effective coalescence radii squared:
     real(real64), private :: coalesRadiiDeutSqrd  = defaultCRDeut(1)  **2
     real(real64), private :: coalesRadiiTritSqrd  = defaultCRTrit(1)  **2
     real(real64), private :: coalesRadiiAlphaSqrd = defaultCRAlpha(1) **2
     real(real64), private :: coalesRadiiLFragSqrd = defaultCRLFrag(1) **2
  end type effectiveCoalescenceRadii


  type, public :: CoalescenceData
     private
     ! Verify proper initialization
     logical, private :: constructed = .FALSE.
     type(coalescenceDataOptions), private :: options
     type(effectiveCoalescenceRadii) :: radii
   contains
     private
     ! Querying object:
     procedure, public  :: properlyConstructedData
     procedure, public  :: properlyConstructed => properlyConstructedData
     procedure, public  :: queryDataOptions
     procedure, public  :: queryOptions => queryDataOptions
     ! Querying data:
     procedure, public  :: coalesRadiiDeut
     procedure, public  :: coalesRadiiTrit
     procedure, public  :: coalesRadiiAlpha
     procedure, public  :: coalesRadiiLFrag
     procedure, public  :: coalesRadiiDeutSqrd
     procedure, public  :: coalesRadiiTritSqrd
     procedure, public  :: coalesRadiiAlphaSqrd
     procedure, public  :: coalesRadiiLFragSqrd
  end type CoalescenceData


  type, public :: coalescenceParticle
     private
     ! Composition of particle
     integer(int32), public  :: numBaryons  = 0_int32      ! Number of baryons in the particle
     integer(int32), public  :: charge      = 0_int32      ! Charge of the particle (number of protons in case of baryon)
     integer(int32), public  :: strangeness = 0_int32      ! Quantum decay number
     ! Energy information of particle
     real(real64),   public  :: kinEnergy   = 0.0_real64   ! Kinetic energy of particle [GeV]
     real(real64),   public  :: restMass    = 0.0_real64   ! Rest mass      of particle [GeV/c**2]
     ! Position (in X/Y/Z coordinates) relative to C.M. system
!     real(real64),   public  :: xCoord      = 0.0_real64   ! X-position [fm]
!     real(real64),   public  :: yCoord      = 0.0_real64   ! Y-position [fm]
!     real(real64),   public  :: zCoord      = 0.0_real64   ! Z-position [fm]
     ! Momentum (in X/Y/Z coordinates)
     real(real64),   public  :: linearMomX  = 0.0_real64   ! Momentum in X [GeV/c]
     real(real64),   public  :: linearMomY  = 0.0_real64   ! Momentum in Y [GeV/c]
     real(real64),   public  :: linearMomZ  = 0.0_real64   ! Momentum in Z [GeV/c]
     ! (Flags for the particles)
     real(real64),   public  :: coalesceFlag = 0.0_real64  ! Flags a fragment that was formed from coalescence
     integer(int32), private :: coalesceNum  = 0_int32     ! Flags pairs of particles that are to coalesce (i.e. when /=0, particle coalesced)
  end type coalescenceParticle


  type, public :: coalescenceResults
     private
     ! Verify proper construction before use
     logical,        private :: constructed = .FALSE.   ! Flags if progeny array was created
     integer(int32), public  :: simState    = 0_int32   ! Flags how the simulation ended

     ! For the particle bank
     type(coalescenceParticle), public, pointer, dimension(:) :: partBnk => NULL()   ! Particle bank
     integer(int32),            public :: numParticles = 0_int32   ! Number of particles existing in the particle bank
     integer(int32),            public :: partBnkSize  = 0_int32   ! Maximum amount of partiles existing
     ! General accounting information
     integer(int32),            public :: numCoalesced       = 0_int32   ! Number of particles that coalesced
     integer(int32),            public :: numFormedFragments = 0_int32   ! Number of fragments formed from coalescing particles
  end type coalescenceResults
