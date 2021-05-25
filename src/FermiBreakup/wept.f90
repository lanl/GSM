
  function wept (fbuObj, ia, iz)
 
! ======================================================================
!
!    Called by: RAZVAL, WECHAN
!
!    Last change: 13-Aug-2003 BY NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Changed by KKG 25 Apr 2006
!    Modified 03 May 2006 by R E Prael to set mass of a "proton star",
!       thus allowing breakup under all conditions.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fermiBreakupParams, only : zro, one, thousandth, &
         & nucleon_mass, electron_mass, proton_mass, unkn_mass
    use fermiBreakUpData, only : dm, min_fermi_AData

    implicit none
    class(FermiBreakup), intent(inout) :: fbuObj
    integer(int32),      intent(in   ) :: ia
    integer(int32),      intent(in   ) :: iz

    real(real64) :: wept

    integer(int32) :: i, j
    real(real64)   :: dmp

! ======================================================================
 
    j = iz
    i = ia - j
    if (i == zro .and. j > one) then      ! ia=iz
!    Modified 03 May 2006 by R E Prael (IA*proton mass)
       wept = dble(ia)*(unkn_mass + thousandth*dm(1,2))
    elseif (i < 0 .or. i > 11 .or. j <= 0 .or. j > 10) then
! Changed by KKG 25 Apr 2006
       wept = proton_mass*dble(ia)
    else
       if ( i+1 > min_fermi_AData .or. j+1 > (min_fermi_AData-1) ) then
          ! dm(X,Y) will be exceeded; warn client and approximate
          write(fbuObj%io%message, 2100) i+1, j+1
          call fbuObj%io%print(3, 3, fbuObj%io%message)
          write(fbuObj%io%message, 2110)
          call fbuObj%io%print(3, 3, fbuObj%io%message)
          if ( i+1 > min_fermi_AData     ) i = min_fermi_AData-1
          if ( j+1 > (min_fermi_AData-1) ) j = min_fermi_AData-2
       end if
       dmp = dm(i+1,j+1)
       if (dmp <= 90.d0) then
          wept = nucleon_mass*dble(ia) + thousandth*dmp - electron_mass*dble(j)
       else
          wept = proton_mass*dble(ia)
       endif
    endif

    return
! ======================================================================
2100 format("The Fermi break-up modle's 'dm' data array was ", &
          & "exceeded (N=", i3, ", Z=", i3, ").")
2110 format("   Using the closest valid value.")
! ======================================================================
  end function wept
