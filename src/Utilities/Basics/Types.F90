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
       f_int8             => int8, &
       f_int16            => int16, &
       f_int32            => int32, &
       f_int64            => int64, &
       f_float            => real32, &
       f_double           => real64, &
       f_long_double      => real128

  ! Import standard C types used (int, float)
  use, intrinsic :: iso_c_binding, only: &
       c_int8             => c_int8_t, &
       c_int16            => c_int16_t, &
       c_int32            => c_int32_t, &
       c_int64            => c_int64_t, &
       ! c_int128           => c_int128_t, &
       c_float, &
       c_double, &
       c_long_double
       ! , c_float128

  ! Import custom types that may be used
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
  integer(f_int8),       parameter :: gsm_int8 = f_int8
  integer(f_int16),      parameter :: gsm_int16 = f_int16
  integer(f_int32),      parameter :: gsm_int32 = f_int32
  integer(f_int64),      parameter :: gsm_int64 = f_int64
  real(f_float),         parameter :: gsm_float = f_float
  real(f_double),        parameter :: gsm_double = f_double
  real(f_long_double),   parameter :: gsm_long_double = f_long_double

end module Types

