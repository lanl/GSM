
  function vcoul (evapObj, z1, a1, iz2, ia2, ck, j)

! ======================================================================
!
!  VCOUL
!    Calculate Coulomb potential
!
! <variables>
!     a1  :   mass of nucleus #1                    (IN)
!     z1  :   charge  of nucleus #1                 (IN)
!    ia2  :   mass of nucleus #2                    (IN)
!    iz2  :   charge  of nucleus #2                 (IN)
!     j   :   index of nucleus #2                   (IN)
!    ck   :   transmission probability              (IN)
!  vcoul  :   Coulomb potential  [MeV]              (OUT)
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, September, 2003.
!    Modified by A. J. Sierk, LANL T-16, January, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, hbarc, ato3rd

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64),       intent(in   ) :: z1
    real(real64),       intent(in   ) :: a1
    integer(int32),     intent(in   ) :: iz2
    integer(int32),     intent(in   ) :: ia2
    real(real64),       intent(in   ) :: ck
    integer(int32),     intent(in   ) :: j
    real(real64)                      :: vcoul

    integer(int32) :: ia1
    real(real64)   :: a1thrd, a2thrd, r0, r1, r2, rzero, temp, ziz

! ======================================================================

    real(real64), parameter :: rc     = 1.70d0
    real(real64), parameter :: ee     = 137.0359895d0
    real(real64), parameter :: vcoul0 = hbarc / ee
    real(real64), parameter :: vcoul1 = vcoul0 / rc

! ======================================================================

!  No Coulomb potential for neutron emission:
    if (j == 1) then
       vcoul = zro
       return
    endif

    ia1 = nint(a1)
    a1thrd = ato3rd(ia1)
    a2thrd = ato3rd(ia2)
    ziz = z1*dble(iz2)
    if (evapObj%options%inverseParameter == -one .and. j > 6) then        
!   Expression in NP A475 (1987) 663 (for heavy ions)
       vcoul = vcoul0
       temp = one + 9.443d-3*ziz
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(evapObj%io%message,1000) '1269'
          call evapObj%io%print(4, 3, evapObj%io%message)
       end if
       rzero = 2.173d0*(one + 6.103d-3*ziz)/(temp)
       r1 = a1thrd
       r2 = a2thrd
       r0 = rzero*(r1 + r2)
    elseif (evapObj%options%inverseParameter < zro .and. j <= 6) then    
!  Formula implemented in FB model             
       vcoul = ck*vcoul0
       r1 = radgem(ia1)                      
       r2 = radgem(ia2)                     
       r0 =  r1 + r2                       
    elseif (evapObj%options%inverseParameter > zro) then
!  Simple parameter set:
       vcoul = vcoul0
       r1 = evapObj%rrVcoul*a1thrd
       r2 = evapObj%rrVcoul*a2thrd
       r0 = r1 + r2
    elseif (evapObj%options%inverseParameter == zro) then
       if (j <= 6) then
! Dostrovsky's parameter set
          vcoul = ck*vcoul1
          temp = rc
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(evapObj%io%message,1000) '1294, 1297'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          r2 = 1.2d0/temp
          if (j == 2) r2 = zro
          r1 = a1thrd
          if (ia1 <= 4 .and. nint(z1) <= 2) r1 = 1.2d0/temp
          r0 = r1 + r2
       else
! Matsuse's parameter set ...PRC 26 (1982) 2338
          vcoul = vcoul0
          temp = a1thrd
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(evapObj%io%message,1000) '1307'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          r1 = 1.12d0*a1thrd - 0.86d0/temp
          temp = a2thrd
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(evapObj%io%message,1000) '1313'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          r2 = 1.12d0*a2thrd - 0.86d0/temp
          r0 = r1 + r2 + 3.75d0
       endif
    endif
    temp = r0
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(evapObj%io%message,1000) '1322'
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    vcoul = vcoul*ziz/temp
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'vcoul.f90', line ", A)
! ======================================================================
  end function vcoul
