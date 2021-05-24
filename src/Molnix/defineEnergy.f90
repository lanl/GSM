
  function defineEnergy (molObj, iz, in, itype) result(newEnergy)

! ======================================================================
!
!   This function extracts either the ground-state microscopic 
!   correction, for itype = 1; the experimental (or theoretical, if
!   Mex is not measured) mass excess for itype = 2; or the energy
!   shift from the Moller, Nix & Kratz calculated pairing gaps, for
!   itype = 3  from the tables contained in the block data subroutine
!   MOLLNIX. The inputs are the proton and neutron numbers, and the 
!   value of itype.
!
!   Calls: MNMACRO
!
!    Written by A. J. Sierk  LANL  T-2  May, 1996.
!    Modifications (AJS); July, October, 1997.
!    Corrected bug (AJS); September, 1998
!   "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use molnixParams, only: zro, one, two

    implicit none
    class(Molnix),  intent(inout) :: molObj
    integer(int32), intent(in   ) :: iz
    integer(int32), intent(in   ) :: in
    integer(int32), intent(in   ) :: itype
    real(real64)                  :: newEnergy

    integer(int32) :: in1, ina, iz1, iz2, nma, nmi, nmx, iz3, iz4
    real(real64)   :: a, eprh1, eprh2, eprl1, eprl2, pair, temp, un, z
    logical        :: nok

! ======================================================================

    un = dble(in)
    z = dble(iz)
    a = z + un
    iz1 = iz - 7
    iz2 = iz - 34
    iz3 = iz - 47
    iz4 = iz - 64
    newEnergy = zro

    ! Determine if within array bounds
    nok = .true.
    if (iz <= 0) then
       return
    elseif (iz < 8) then
       nok = in >= nminaData(iz) .and. in <= nmaxaData(iz)
       ina = in - nminaData(iz) + 1
       nmi = nminaData(iz)
       nma = nmaxaData(iz)
       nmx = nmaxaData(iz) - nminaData(iz) + 1
    elseif (iz > 100) then
       nok = .false.
    else
       nok = in >= nminData(iz1) .and. in <= nmaxData(iz1)
       in1 = in - nminData(iz1) + 1
       nmi = nminData(iz1)
       nma = nmaxData(iz1)
       nmx = nmaxData(iz1) - nminData(iz1) + 1
    endif
    
    if (nok) then
       if (iz < 8) then
          if (itype == 1) then
             newEnergy = emica(ina,iz)
          elseif (itype == 2) then
             newEnergy = emxa(ina,iz)
          elseif (itype == 3) then
             newEnergy = epaira(ina,iz)
          endif
       elseif (iz < 35 .and. iz >= 8) then
          if (itype == 1) then
             newEnergy = emic1(in1,iz1)
          elseif (itype == 2) then
             newEnergy = emx1(in1,iz1)
          elseif (itype == 3) then
             newEnergy = epair1(in1,iz1)
          endif
       elseif (iz < 48 .and. iz >= 35) then
          if (itype == 1) then
             newEnergy = emic2(in1,iz2)
          elseif (itype == 2) then
             newEnergy = emx2(in1,iz2)
          elseif (itype == 3) then
             newEnergy = epair2(in1,iz2)
          endif
       elseif (iz < 65 .and. iz >= 48) then
          if (itype == 1) then
             newEnergy = emic3(in1,iz3)
          elseif (itype == 2) then
             newEnergy = emx3(in1,iz3)
          elseif (itype == 3) then
             newEnergy = epair3(in1,iz3)
          endif
       elseif (iz < 101 .and. iz >= 65) then
          if (itype == 1) then
             newEnergy = emic4(in1,iz4)
          elseif (itype == 2) then
             newEnergy = emx4(in1,iz4)
          elseif (itype == 3) then
             newEnergy = epair4(in1,iz4)
          endif
       endif
    else
!  (nok is false):if IN is outside the table for this value of IZ
!  or for IZ > 100:

       if (in < 0) return
       if (itype == 2) then
          newEnergy = molObj%mnmacro(iz, in)
       elseif (itype == 1) then
!  For nuclei outside the table, set the microscopic correction to that
!  for the first or last nuclide (for this value of iz) in the table.
          if (iz < 8) then
             if (in < nminaData(iz) .and. in > 0) newEnergy = emica(ina, iz)
             if (in > nmaxaData(iz)) newEnergy = emica(nmx, iz)
          elseif (iz >= 8 .and. iz < 35) then
             if (in < nminData(iz1) .and. in > 0) newEnergy = emic1(1, iz1)
             if (in > nmaxData(iz1)) newEnergy = emic1(nmx, iz1)
          elseif (iz >= 35 .and. iz < 48) then
             if (in < nminData(iz1) .and. in > 0) newEnergy = emic2(1, iz2)
             if (in > nmaxData(iz1)) newEnergy = emic2(nmx, iz2)
          elseif (iz >= 48 .and. iz < 65) then
             if (in < nminData(iz1) .and. in > 0) newEnergy = emic3(1, iz3)
             if (in > nmaxData(iz1)) newEnergy = emic3(nmx, iz3)
          elseif (iz >= 65 .and. iz < 101) then
             if (in < nminData(iz1) .and. in > 0) newEnergy = emic4(1, iz4)
             if (in > nmaxData(iz1)) newEnergy = emic4(nmx, iz4)
          elseif (iz >= 101) then
             newEnergy = zro
          endif
       elseif (itype == 3) then
          if (iz < 8) then
             eprl1 = epaira(1, iz)
             eprl2 = epaira(2, iz)
             eprh1 = epaira(nmx-1, iz)
             eprh2 = epaira(nmx, iz)
          elseif (iz < 35 .and. iz >= 8) then
             eprl1 = epair1(1, iz1)
             eprl2 = epair1(2, iz1)
             eprh1 = epair1(nmx-1, iz1)
             eprh2 = epair1(nmx, iz1)
          elseif (iz < 48 .and. iz >= 35) then
             eprl1 = epair2(1, iz2)
             eprl2 = epair2(2, iz2)
             eprh1 = epair2(nmx-1, iz2)
             eprh2 = epair2(nmx, iz2)
          elseif (iz < 65 .and. iz >= 48) then
             eprl1 = epair3(1, iz3)
             eprl2 = epair3(2, iz3)
             eprh1 = epair3(nmx-1, iz3)
             eprh2 = epair3(nmx, iz3)
          elseif (iz < 101 .and. iz >= 65) then
             eprl1 = epair4(1, iz4)
             eprl2 = epair4(2, iz4)
             eprh1 = epair4(nmx-1, iz4)
             eprh2 = epair4(nmx, iz4)
          elseif (iz >= 101) then
             pair = (one - z + two*dble(iz/2)) +  &
                  & (one - un + two*dble(in/2))
             temp = a
             if (temp < div0Lim .and. temp > -div0Lim) then
                temp = div0Lim
                write(molObj%io%message,1000) "183"
                call molObj%io%print(4, 3, molObj%io%message)
             end if
             newEnergy = molObj%options%cevap*pair/sqrt(abs(temp))
             return
          endif
          if (in < nmi .and. in > 0) then
             if (in == nmi-1 .or. in == nmi-3 .or. in == nmi-5) then
                newEnergy = eprl2
             elseif (in == nmi-2.or.in == nmi-4 .or. in == nmi-6) then
                newEnergy = eprl1
             else
                pair = (one - z + two*dble(iz/2)) +  &
                     & (one - un + two*dble(in/2))
                temp = a
                if (temp < div0Lim .and. temp > -div0Lim) then
                   temp = div0Lim
                   write(molObj%io%message,1000) "199"
                   call molObj%io%print(4, 3, molObj%io%message)
                end if
                newEnergy = molObj%options%cevap*pair/sqrt(abs(temp))
             endif
          elseif (in > nma) then
             if (in == nma+1 .or. in == nma+3 .or. in == nma+5) then
                newEnergy = eprh1
             elseif (in == nma+2.or.in == nma+4 .or. in == nma+6) then
                newEnergy = eprh2
             else
                pair = (one - z + two*dble(iz/2)) +   &
                     (one - un + two*dble(in/2))
                temp = a
                if (temp < div0Lim .and. temp > -div0Lim) then
                   temp = div0Lim
                   write(molObj%io%message,1000) "214"
                   call molObj%io%print(4, 3, molObj%io%message)
                end if
                newEnergy = molObj%options%cevap*pair/sqrt(abs(temp))
             endif
          endif
       endif
    endif
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'defineEnergy.f90', ", &
          & "line(s) ", A)
! ======================================================================
  end function defineEnergy
