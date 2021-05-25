
! ==============================================================================
!
! This file contains all data objects for each of the submodels
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

  ! Contains all general physics data:
  type, private :: generalPhysicsData
     private
     integer(int32), private :: numConstructionErrors = 3_int32

     type(Molnix),              private :: molnixEnergies
     type(FissionBarrier),      private :: fissBarrier
     type(EvaporationData),     private :: evap
  end type generalPhysicsData


  ! Contains all reaction specific data:
  type, private :: rxnSpecificData
     private
     integer(int32), private :: numConstructionErrors = 3_int32

     type(StandardDCMData),     private :: sDCMTargetData
     type(CoalescenceData),     private :: coales
     type(PreequilibriumData),  private :: preeq
  end type rxnSpecificData
