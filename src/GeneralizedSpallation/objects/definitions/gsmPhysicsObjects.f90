 
! ==============================================================================
!
! This data type contains all of the model objects for construction of GSM
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================


  ! This object contains all models that are general (require no rxn-specific data)
  type, public :: generalPhysicsModels
     private
     integer(int32),       public :: numConstructionErrors = 3_int32

     type(StandardDCM),    private :: sDCM
!     type(ModifiedDCM),    private :: mDCM
     type(Evaporation),    private :: evap
     type(FermiBreakup),   private :: fbu
  end type generalPhysicsModels


  ! This object contains all reaction-specific models (Preeq., Coales.)
  type, private :: rxnSpecificModels
     private
     integer(int32), private :: numConstructionErrors = 0_int32

     type(Coalescence),    private :: coales
     type(Preequilibrium), private :: preeq
  end type rxnSpecificModels
