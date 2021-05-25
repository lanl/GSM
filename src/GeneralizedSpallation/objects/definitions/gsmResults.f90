
! ==============================================================================
!
! This data type contains all of the results objects and variables that a single
! simulation of GSM would need.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

  ! Information tracked for a single particle
  type, public :: GSMProgeny
     private
     ! Particle composition
     real(real64), public :: numBaryons = 0.0_real64   ! Number of baryons in fragment (A number)
     real(real64), public :: numProtons = 0.0_real64   ! Charge of fragment (Z number)
     real(real64), public :: kinEnergy  = 0.0_real64   ! Kinetic energy of fragment [GeV]
     real(real64), public :: restMass   = 0.0_real64   ! Rest mass of particle [GeV/c**2]

     ! Particle direction
     real(real64), public :: phi        = 0.0_real64   ! Phi   of direction of movement (in p-space)
     real(real64), public :: theta      = 0.0_real64   ! Theta of direction of movement (in p-space)
     real(real64), public :: sinTheta   = 0.0_real64   ! sin(theta)
     real(real64), public :: cosTheta   = 0.0_real64   ! cos(theta)

     ! Other particle information
     real(real64), public :: typeID     = 0.0_real64   ! Type #; 1-9 for n-pi+, else 1000*Z + N
     real(real64), public :: prodMech   = 0.0_real64   ! Origin of fragment (production mechanism)
                             ! <    0 for hole
                             ! <  100 for INC particle
                             ! =  100 for pre-equilibrium particles
                             ! = 1000 for evaporation from spallation residue (and residue itself)
                             ! = 1500 for Fermi Break-Up
                             ! = 2000 for evaporation from fission product (and fragments themselves)
  end type GSMProgeny


  ! Regarding created exciton information:
  type, private :: excitonData
     private
     integer(int32), public :: numTotal    = 0_int32   ! Total number of created excitons (sum of three below)
     integer(int32), public :: numProtons  = 0_int32   ! Number of exciton charges
     integer(int32), public :: numNeutrons = 0_int32   ! Number of exciton neutrons
     integer(int32), public :: numHoles    = 0_int32   ! Number of holes
  end type excitonData


  ! Information tracked in all of GSM as results
  type, public :: GSMResults
     private

     ! Used for easy access to the initial nuclei and current compound nuclei:
     ! (intial reactants)
     type(GSMProjectile), public,  pointer :: initialProj
     type(GSMTarget),     public,  pointer :: initialTarg
     ! (residual nuclei during/after simulation)
     type(GSMResidual),   public  :: projRes
     type(GSMResidual),   public  :: targRes
     ! (exciton information from the residual)
     type(excitonData),   private :: projExc
     type(excitonData),   private :: targExc


     ! Progeny information:
     type(GSMProgeny),    public, dimension(:), pointer :: progenyBnk => NULL()
     integer(int32),      public  :: numProgeny   =   0_int32
     integer(int32),      private :: maxProgeny   =   0_int32
     integer(int32),      private :: maxProgenyM1 =   0_int32


     ! Misc. reaction information:
     integer(int32),      public  :: numElasticEvents = 0_int32
     type(EventData),     public  :: info


     ! For counting of model usage and other misc. information:
     type(GSMModelUsage), public  :: modelUsage

     ! Regarding simulation state (warnings, errors) and crash-protection:
     integer(int32),      public  :: simState         = successfulSingleEvent
     logical,             private :: constructed      = .FALSE.   ! Flags if object was initialized properly
     logical,             private :: tallySim         = .FALSE.

   contains
     procedure, private :: resetEvent
  end type GSMResults
