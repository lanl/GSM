
! ==============================================================================
!
! Specifies the data type 'OutputData'
!
!
! Written by CMJ, XCP-3, 11/2019
!
! ==============================================================================

  ! Tracks event-specific data for output files
  type, public :: OutputData
     private

     ! From /dele/ block
     real(real64), public :: &
          & sfu    = 0.0_real64   ! Counter for average fission probability

     ! From the /fiscem/ block
     integer(int64), public :: nfis = 0_int64

     ! From /vul/ block
     real(real64), public :: sigom = 0.0_real64
     integer(int64), public :: &
          & ncas  = 0_int64, &   ! Num. successful inelastic events
          & intel = 0_int64, &   ! Num. total elastic events
          & limcc = 0_int64      ! Total inelastic events to try
     
     ! For coalescence information
     integer(int32), dimension(8), public :: ncoal = 0_int32
     real(real64),   dimension(4), public :: pcoal = 0_int32

     ! Counts the number of times Fermi Break-up was used
     integer(int64), public :: ifermi = 0_int64

   contains
     private

     ! Getters

     ! Setters

     ! Introspection

  end type OutputData

