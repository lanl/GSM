
! ==============================================================================
!
! Contains information about the reaction data class.
!> \todo Create inidivual modules for each of the GSM objects (attr too)
!> \todo Parse \c gsmResults into several objects, pass around only needed info
!> \todo Make \c Reactions module that points to proj/targ and 
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

  ! NOTE: This object is meant to simply POINT to the various reaction-specific variables
  !       to simplify the passing of argumetns, and (hopefully) without adding much memory
  !       requirements (i.e. POINT to all [large] variables/objects!)
  ! NOTE: GSM will query this object DIRECTLY (not through 'getter' methods).
  type, private :: GSMReaction
     private
     ! Inputs to GSM when simulating:
     type(GSMResults),    private, pointer :: results  => NULL()

     ! Model physics:
     type(rxnSpecificData),   private :: data
     type(rxnSpecificModels), private :: models

     ! Interim values:
     integer(int32),       private :: incFlag = defaultINCFlag
     integer(int64),       private, pointer :: eventNum => NULL()
     type(sDCMProjectile), private, pointer :: sDCMProj => NULL()

     ! Tallying
     integer(int64),      private :: eventAttempt = 0_int64

     ! Output
     type(GSMOutput),     private, pointer :: output => NULL()
     type(OutputData),    private, pointer :: outData  => NULL()

     ! Misc:
     logical, private :: constructed = .FALSE.

   contains
     private
     procedure, public :: properlyConstructed => reactionProperlyConstructed
  end type GSMReaction
