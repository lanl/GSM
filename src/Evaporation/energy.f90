
  function energy (evapObj, iz, ia)

! ======================================================================
!
!****************This routine was originally in the LAHET code********
!  ENERGY
!   Calculates mass excess using OLD Cameron formula when outside
!   Wapstra table!!!!!  (Based on O-16 = 0.00 binding energies!)
!
!    Modified to use Moller-Nix Masses for N, Z >= 8; macroscopic
!    M-N masses outside the table, and Wapstra tables only for N,Z
!   < 8 and masses not measured.  AJS
!
!   For light nuclei outside the Wapstra systematics table, relies on
!   the old Cameron calculation.....VERY IMPROBABLE!
!  Z = 1, A > 12;  Z = 2, A > 14;  Z = 3, A > 17;  Z = 4, A > 22;
!  Z = 5, A > 25;  Z = 6, A > 28;  Z = 7, A > 31.
! 
! =====================================================================
! <variables>
!     ia   :   mass number of nucleus               (IN)
!     iz   :   charge of nucleus                    (IN)
!  energy  :   Mass excess [MeV]                    (OUT)
!
!   Last change: 13-AUG-2003 by NVMokhov
!   Modified by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by CMJ, XCP-3, July 2018 (Evap. class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, ato3rd
    use evaporationFissionData, only: wapsm, inn, iiz, sz, sn

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    integer(int32),     intent(in   ) :: iz
    integer(int32),     intent(in   ) :: ia
    real(real64)                      :: energy

    integer(int32) :: in
    real(real64)   :: a, a3, z
    logical        :: nsmal, zsmal

! ======================================================================

    energy = 1.d10
    in = ia - iz

    ! Limits in program
    ! IA = 4, IZ = 1, IN = 3
    if ( (iz < zro .OR. iiz < iz) .OR. (in < zro .OR. inn < in) ) then
       ! Outside function limits, warn user
       write(evapObj%io%message,1000) ia, iz, in
       call evapObj%io%print(2, 3, evapObj%io%message)
    endif
    if (iz < 0 .or. in < 0) return

    zsmal = iz < 8 .or. in < 8
    if (iz == 0) then
       nsmal = .true.
       if (in == 0) return
    elseif (iz < 8) then
       nsmal = in < evapObj%evapMolnix%nmina(iz) .or. in > evapObj%evapMolnix%nmaxa(iz)
    elseif (iz >= 8) then
       nsmal = in < evapObj%evapMolnix%nmin(iz-7) .or. in > evapObj%evapMolnix%nmax(iz-7)
    endif

    ! Obtain mass excess values
    if (.not.zsmal .or. (zsmal .and. .not.nsmal)) then
       energy = evapObj%evapMolnix%defineEnergy (iz, in, 2)
    elseif (zsmal .and. nsmal) then
       energy = wapsm(iz,in)
       if (energy == 0.d0) then
          if (iz == 6 .and. ia == 12) return
          if (iz > 0 .and. in > 0) then
!   Cameron's mass formula for nuclides outside the Wapstra table:
!   ONLY for light nuclei outside the Wapstra table in mass.tbl.
             a = dble(ia)
             z = dble(iz)
             a3 = ato3rd(ia)
             energy = evapObj%cam (a, z, a3) + sz(iz) + sn(in)
          endif
       endif
    endif

    return
! ======================================================================
1000 format("Secondary array exceeded in 'energy.f90', ", &
          & "line 55 (A=", i3, ", Z = ", i3, " [N = ", i3, "]).")
! ======================================================================
  end function energy


  function cam (evapObj, a, z, b)

! ======================================================================
!
!   Cameron mass-excess formula; extracted from inline function from
!   ENERGY. Formula originally in LAHET. b = A**(1/3).
!
!   Written by A. J. Sierk, LANL T-16, September, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection).
!   Edited by CMJ, XCP-3, July 2018 (Evap. class creation)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: one, two, ato3rd

    implicit none
    class(Evaporation), intent(inout) :: evapObj
    real(real64),       intent(in   ) :: a
    real(real64),       intent(in   ) :: z
    real(real64),       intent(in   ) :: b
    real(real64)                      :: cam

    integer(int32) :: iz
    real(real64)   :: abi, zft, tempa, tempb

! ======================================================================

    iz = nint(z)
    tempa = a
    if ( tempa < div0Lim .and. tempa > -div0Lim) then
       tempa = div0Lim
       write(evapObj%io%message,1000) "134, 139-141"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    abi = one - two*z/tempa
    zft = ato3rd(iz)**4
    tempb = b
    if ( tempb < div0Lim .and. tempb > -div0Lim) then
       tempb = div0Lim
       write(evapObj%io%message,1000) "139-141"
       call evapObj%io%print(4, 3, evapObj%io%message)
    end if
    cam = 8.071323d0*tempa - 0.782354d0*z - 17.0354d0*tempa*(one - &
         & 1.84619d0*abi**2) + 25.8357d0*tempb**2*(one - &
         & 1.71219d0*abi**2)*(one - 0.62025d0/tempb**2)**2 + &
         & 0.779d0*z*(z - one)*(one - 1.5849d0/tempb**2 + &
         & 1.2273d0/tempa + 1.5772d0/(tempa*tempb))/tempb - 0.4323d0 * &
         & zft*(one - 0.57811d0/tempb - &
         & 0.14518d0/tempb**2 + 0.49597d0/tempa)/tempb

    return

! ======================================================================
1000 format("Divide by zero error prevented in 'energy.f90,' line ", A)
! ======================================================================
  end function cam

