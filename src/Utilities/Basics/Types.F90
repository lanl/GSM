! =============================================================================
!>
!> \file
!> \brief  Contains the Types module
!> \author CMJ
!>
!> The Types module contains an enumeration of types available in Fortran.
!> The module exists to alias the available types, specifically those utilized
!> by GSM and its sub-libraries.
!>
! =============================================================================
module Types

  ! Import standard Fortran types used (int, float)
  use, intrinsic :: iso_fortran_env, only: &
       fInt8             => int8, &
       fInt16            => int16, &
       fInt32            => int32, &
       fInt64            => int64, &
       fFloat            => real32, &
       fDouble           => real64, &
       fLongDouble       => real128

  ! Import standard C types used (int, float)
  use, intrinsic :: iso_c_binding, only: &
       cInt8             => c_int8_t, &
       cInt16            => c_int16_t, &
       cInt32            => c_int32_t, &
       cInt64            => c_int64_t, &
       ! cInt128           => c_int128_t, &
       cFloat            => c_float, &
       cDouble           => c_double, &
       cLongDouble      => c_long_double
  ! , cFloat128

  ! Import custom types that may be used
  !> \todo Move these to a FileStream module
  use, intrinsic :: iso_fortran_env, only: &
       STAT_END           => iostat_end, &  ! Stream end
       STAT_EOR           => iostat_eor, &  ! Stream record end
       STAT_IS_INTERNAL   => iostat_inquire_internal_unit, & ! Flag for an "internal" file unit
       error_unit, &
       output_unit, &
       input_unit

  use, intrinsic :: iso_c_binding, only: &
       c_bool, &
       c_char, &
       ! Character-related items:
       c_alert, &
       c_new_line, &
       c_horizontal_tab, &
       c_vertical_tab
  
  implicit none
  public

  ! GSM Types
  integer(fInt8),       parameter :: gsmInt8         = fInt8
  integer(fInt16),      parameter :: gsmInt16        = fInt16
  integer(fInt32),      parameter :: gsmInt32        = fInt32
  integer(fInt64),      parameter :: gsmInt64        = fInt64
  integer(fFloat),      parameter :: gsmFloat        = fFloat
  integer(fDouble),     parameter :: gsmDouble       = fDouble
  integer(fLongDouble), parameter :: gsmLongDouble   = fLongDouble

end module Types

