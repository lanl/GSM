
! ==============================================================================
!
! This file contains all of the data types used by the modified Dubna Cascade
! Model's Data Class
!
!
! Written by CMJ, XCP-3, 02/2024
!
! ==============================================================================

  ! Create an IO type
  type, private :: mDCMDataIO
     private
     procedure(IOHandler), private, nopass, pointer :: print   => printMDCMData
     character(LEN=512),   private                  :: message =  ""
  end type mDCMDataIO


  ! Options for the data class
  type, public :: mDCMDataOptions
     private


  end type mDCMDataOptions

