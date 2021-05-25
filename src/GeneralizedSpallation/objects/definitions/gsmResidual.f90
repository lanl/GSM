
! ==============================================================================
!
! This file contains the nuclei data types for GSM regarding initial projectile,
! target, and the resulting residual.
!
!
! Written by CMJ, XCP-3, 01/2019
!
! ==============================================================================

  type, public :: gsmResidual
     private
     ! Regarding nucleonic composition...
     real(real64),   public :: numBaryons = 0.0_real64
     real(real64),   public :: numProtons = 0.0_real64

     ! Regarding kinetic energy AFTER collission in each event
     real(real64),   public :: kinEnergy = 0.0_real64   ! [MeV]

     ! Regarding momentum (both linear and angular)...
     real(real64),   public, dimension(3) :: linearMom  = 0.0_real64   ! [GeV/c]
     real(real64),   public, dimension(3) :: angularMom = 0.0_real64   ! [GeV/c]
   contains
     private
     procedure, public  :: totalLinearMomentum
     procedure, public  :: totalAngularMomentum
  end type gsmResidual
