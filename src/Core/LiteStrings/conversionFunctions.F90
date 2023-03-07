! =============================================================================
!
!> \file
!> \brief  Contains the numerical conversion functions
!> \author CMJ
!
! =============================================================================


! =============================================================================
!> \brief Returns the format an integer should be displayed in
function int8ToString(val) result(str)
  implicit none
  integer(gsmInt8), intent(in   ) :: val
  character(:), allocatable :: str

  write(str, intFormat) val

end function int8ToString


! =============================================================================
!> \brief Returns the format an integer should be displayed in
function int16ToString(val) result(str)
  implicit none
  integer(gsmInt16), intent(in   ) :: val
  character(:), allocatable :: str

  write(str, intFormat) val

end function int16ToString


! =============================================================================
!> \brief Returns the format an integer should be displayed in
function int32ToString(val) result(str)
  implicit none
  integer(gsmInt32), intent(in   ) :: val
  character(:), allocatable :: str

  write(str, intFormat) val

end function int32ToString


! =============================================================================
!> \brief Returns the format an integer should be displayed in
function int64ToString(val) result(str)
  implicit none
  integer(gsmInt64), intent(in   ) :: val
  character(:), allocatable :: str

  write(str, intFormat) val

end function int64ToString


! =============================================================================
!> \brief Returns the format a real number should be displayed in
function floatToString(val) result(str)
  implicit none
  real(gsmFloat), intent(in   ) :: val
  character(:), allocatable :: str

  if (useScientificForm) then
     write(str, scientificFormat) val
  else
     write(str, realFormat) val
  end if
end function floatToString


! =============================================================================
!> \brief Returns the format a real number should be displayed in
function doubleToString(val) result(str)
  implicit none
  real(gsmDouble), intent(in   ) :: val
  character(:), allocatable :: str

  if (useScientificForm) then
     write(str, scientificFormat) val
  else
     write(str, realFormat) val
  end if
end function doubleToString


! =============================================================================
!> \brief Returns the format a real number should be displayed in
function longDoubleToString(val) result(str)
  implicit none
  real(gsmLongDouble), intent(in   ) :: val
  character(:), allocatable :: str

  if (useScientificForm) then
     write(str, scientificFormat) val
  else
     write(str, realFormat) val
  end if
end function longDoubleToString


! =============================================================================
!> \brief Returns the format a logical should be displayed in
function logicalToString(val) result(str)
  implicit none
  logical, intent(in   ) :: val
  character(:), allocatable :: str

  write(str, logicalFormat) val
end function logicalToString


! =============================================================================
!> \brief Returns the format a string should be displayed in
function stringToString(val) result(str)
  implicit none
  character(LEN = *), intent(in   ) :: val
  character(:), allocatable :: str

  write(str, stringFormat) val
end function stringToString
