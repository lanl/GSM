! =============================================================================
!
!> \file
!> \brief  Contains the Insist macro
!> \author CMJ (XCP-3; LANL)
!
!> \def Insist
! =============================================================================

#ifndef Insist
#define Insist(COND, MSG) \
   if (.not. (COND)) insist(COND, MSG, __FILE__, __LINE__)
#endif
