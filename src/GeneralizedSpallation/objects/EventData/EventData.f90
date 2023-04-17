
! ==============================================================================
!
! Specifies the data type 'EventData'
!
!
! Written by CMJ, XCP-3, 11/2019
!
! ==============================================================================

  ! Tracks event-specific data for output files
  type, public :: EventData
     private

     ! From the /dele/ block
     real(real64), public :: &
          & wf     =  0.0_real64, &   ! Fission probability
          & fusion = -1.0_real64      ! Flags fission (1), not (0)

     !> \brief Residual prior to/post Preequilibrium
     !> NOTE: Data from the /resid0/ block
     type(PreeqResidual), public, dimension(2):: residual

     !> \brief Compound nuclues from Evaporation process
     !> NOTE: Data from the /residf/ block
     type(EvapCompound), public:: compound

     !> \brief Flags if fission occurred (Flagged also in \c fissFrag below
     integer(int32), public:: ifiss = 0_int32

     !> \brief Fission fragments from Evaporation physics
     !> NOTE: Data from the /ifissc/ block
     type(FissionFragment), public, dimension(2):: fissFrag

   contains
     private

     ! Getters

     ! Setters

     ! Introspection

  end type EventData

