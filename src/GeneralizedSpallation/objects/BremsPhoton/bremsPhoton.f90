
! ==============================================================================
!
! Specifies the data type 'BremsPhoton'
!
!
! Written by CMJ, XCP-3, 11/2019
!
! ==============================================================================

  ! Tracks event-specific data for output files
  type, public :: BremsPhoton
     private
     
     real(real64), public:: &
          & tgmin   = 0.0_real64, &  ! Min KE
          & tgmax   = 0.0_real64, &  ! Max KE
          & tEquiv  = 0.0_real64, &  ! Equivalent T (for a set of events)
          & absSX   = 0.0_real64     ! Absolute scattering cross section (for a set of events)

   contains
     private

     ! Getters
     procedure, public  :: getTMin
     procedure, public  :: getTMax
     procedure, public  :: getTEqv
     procedure, public  :: getSXAbs
     generic,   public  :: tMin  => getTMin
     generic,   public  :: tMax  => getTMax
     generic,   public  :: TEqv  => getTEqv
     generic,   public  :: sxabs => getSXAbs

     ! Setters
     procedure, private :: setTMin
     procedure, private :: setTMax
     procedure, private :: resetTEqv
     procedure, private :: resetSXAbs

     ! Introspection
     procedure, public  :: incrementTEQV
     procedure, public  :: incrementSXABS
     generic,   public  :: addTEQV  => incrementTEQV
     generic,   public  :: addSXABS => incrementSXABS
     procedure, public  :: sampleEnergy

  end type BremsPhoton

