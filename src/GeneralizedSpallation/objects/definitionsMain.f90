
! =============================================================================
!
!> \file
!> \brief   Contains the main include files for all small derived type definitions
!> \author  CMJ, XCP-3 (LANL)
!
! ==============================================================================

   ! Counting of model usage
   include "objects/definitions/gsmModelUsage.f90"

   ! Contains data types for particle definitions
   include "objects/definitions/gsmProjectile.f90"
   include "objects/definitions/gsmTarget.f90"
   include "objects/definitions/gsmResidual.f90"

   ! Contains the GSM object attributes
   include "objects/definitions/gsmOptions.f90"

   ! Contains the GSM output attributes
   include "objects/definitions/gsmOutput.f90"

   ! Contains the composite GSM data objects
   include "objects/definitions/gsmDataObjects.f90"

   ! Contains the composite GSM physics objects
   include "objects/definitions/gsmPhysicsObjects.f90"

   ! Contains the main results object
   include "objects/definitions/gsmResults.f90"

   ! Contains the reaction object (composite for results, output, etc.)
   include "objects/definitions/gsmReaction.f90"
