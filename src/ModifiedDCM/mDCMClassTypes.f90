
! ====================================================================
!
! This file contains the derived types used by the mDCM
!
! ====================================================================

  type, private :: mDCMIO
     private
     character(LEN=512), private :: message = ""
     procedure(IOHANDLER), private, nopass, pointer :: print => printMDCM
  end type mDCMIO


  type, public :: mDCMOptions
     private
     real(real64), public :: rm    = defaultRm
     real(real64), public :: delta = defaultDelta
  end type mDCMOptions


  ! Replaces the "memoryLAQ" common block (once in array)
  type, public :: mDCMProgeny
     private
     ! Particle ID:
     integer(int32), public  :: numProtons  = 0_int32      ! Number of protons  in the nucleus
     integer(int32), public  :: numBaryons  = 0_int32      ! Number of nucleons in the nucleus
     integer(int32), public  :: strangeness = 0_int32      ! Quantum decay number of particle
     integer(int32), public  :: photonFlag  = 0_int32      ! Flags if particle is a photon (>0) or not (=0)
     ! Energy:
     real(real64),   public  :: kinEnergy   = 0.0_real64   ! Kinetic energy of the particle [GeV]
     real(real64),   public  :: restMass    = 0.0_real64   ! Rest energy of the particle [GeV/c**2]
     ! Position:
     real(real64),   public, dimension(3) :: position  = 0.0_real64   ! 3-D position, in cartesian coordiantes [fm]
     ! Linear Momentum:
     real(real64),   public, dimension(3) :: linearMom = 0.0_real64   ! 3-D momentum, in cartesian coordinates
     ! Emission Angles:
     real(real64),   public  :: theta         = 0.0_real64
     real(real64),   public  :: phi           = 0.0_real64
     ! Time Values:
     real(real64),   public  :: formationTime = 0.0_real64
     integer(int32), public  :: timeParameter = 0_int32   ! Some time parameter (specifics unknown)
  end type mDCMProgeny


  ! NOTE: This will replace the "resultLAQ" common block
  type, public :: mDCMResidual
     private
     real(real64), public :: numBaryons = 0.0_real64
     real(real64), public :: numProtons = 0.0_real64
     real(real64), public :: kinEnergy  = 0.0_real64   ! Kinetic energy of fragment [GeV]
     real(real64), public, dimension(3) :: linearMom  = 0.0_real64
     real(real64), public, dimension(3) :: angularMom = 0.0_real64
  contains
     private
     ! angular momentum modifer
  end type mDCMResidual


  ! Review "holpt" common block - this may need incorporated here or a similar
  ! object
  type, public :: mDCMExciton
     private
     integer(int32), public :: numTotal    = 0_int32
     integer(int32), public :: numProtons  = 0_int32
     integer(int32), public :: numHoles    = 0_int32
  end type mDCMExciton


  type, public :: mDCMResults
     private
     ! Progeny information:
     integer(int32), public  :: numProgeny = 0_int32
     integer(int32), public  :: maxProgeny = 0_int32   ! CHANGE TO PRIVATE LATER!
     integer(int32), public  :: maxProgenyM1 = 0_int32   ! CHANGE TO PRIVATE LATER!
     type(mDCMProgeny), public, pointer, dimension(:) :: progenyBnk => NULL()

     ! Residual information:
     type(mDCMResidual), public  :: targRes
     type(mDCMResidual), public  :: projRes

     ! Exciton information:
     type(mDCMExciton), public  :: targExc
     type(mDCMExciton), public  :: projExc

     ! Misc. simulation data:
     integer(int32), public  :: numElastic = 0_int32
     integer(int32), public  :: simState = 0_int32

     ! Misc. information:
     logical,        private :: constructed = .FALSE.
  end type mDCMResults
