
! ====================================================================
!
! This file contains the derived types used by the mDCM
!
! ====================================================================


  ! NOTE: This will replace the "resultLAQ" common block
  type, public :: mDCMResidual
     private
     real(real64), private :: numBaryons = 0.0_real64
     real(real64), private :: numProtons = 0.0_real64
     real(real64), private :: kinEnergy  = 0.0_real64   ! Kinetic energy of fragment [GeV]
     real(real64), private, dimension(3) :: linearMom  = 0.0_real64   ! Units?
     real(real64), private, dimension(3) :: angularMom = 0.0_real64   ! Units?
  contains
     private
     ! Methods to retrieve values
     procedure, public :: getNumBaryons => residualBaryons
     procedure, public :: getNumProtons => residualProtons
     procedure, public :: getKinEnergy => residualKinEnergy
     procedure, public :: getLinearMom => residualLinearMom
     procedure, public :: getAngularMom => residualAngularMom

     ! Methods to adjust values
     procedure, public  :: adjustResidual
     procedure, public  :: adjustMomentum
  end type mDCMResidual

