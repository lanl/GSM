
! ==============================================================================
!
!> \file
!> \brief   Contains the implementation of \c gsmTarget
!> \author  CMJ, XCP-3 (LANL)
!
! ==============================================================================
!  DEFAULT VALUES
! ==============================================================================

! --------------------------------------------- TARGET NAME --------------------
  
  !> Default target name
  character(len=*), public, parameter :: defaultTargName     = "Pb-208"

  !> Default number of baryons (nucleons) in the target nucleus
  real(real64),     public, parameter :: defaultTargBaryons  = 208.0_real64

  !> Default number of protons (i.e. charge) in the target nucleus
  real(real64),     public, parameter :: defaultTargProtons  =  82.0_real64

! ==============================================================================
!
!  DATA TYPE DESCRIPTION:
!
!> \class  gsmTarget
!> \brief  Defines the target object utilized by clients for GSM simulations
!
!> Defines the target object utilized by clients for simulations of GSM. The
!> \c gsmTarget object defines characteristics of the static particle (in the
!> lab. system) during a high energy nuclear event.
!
!> Supported target particles include any nucleus with a real, positive nucleus
!> (being defined as having a positive number of protons and neutrons), namely
!> consisting of light- and heavy-ions. Note the \c StandardDCM object does NOT
!> well handle light nuclei (\f$ A_{t} \leq 4 \f$)
!
!> The data type is passed throughout the simulation object, being stored only by
!> local variables. The object is only pointed to by these variables, thus
!> minimizing any data replication where possible. Note the provided \c
!> gsmTarget object may be altered by GSM during the simulation depending on
!> the specification of the object by the software client.
!
! ==============================================================================
  type, public :: gsmTarget
     private

     !> Name of the target particle
     character(LEN=maxPartNameLen), public :: particleName = defaultTargName

     !> Number of baryons (nucleons, A number) in the target
     real(real64),     public :: numBaryons   = defaultTargBaryons

     !> Number of protons (charge, Z number) in the target
     real(real64),     public :: numProtons   = defaultTargProtons

     !> Rest mass of the target nucleus
     real(real64),     public :: restMass     = defaultRestMass

     !> The target particle's \f$ a_{f} \f$ fission scaling multiplier
     real(real64),     public :: afMultiplier = defaultAfMultiplier

     !> The target particle's \f$ C \left(Z\right) \f$ fission multiplier
     real(real64),     public :: czMultiplier = defaultCzMultiplier

   contains
     private

     !> @{
     !> Describes the particle to a provided I/O stream
     procedure, private :: describeTarget
     procedure, public  :: describe => describeTarget
     !> @}

  end type gsmTarget

