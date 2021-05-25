
  subroutine setLevelDensity(apar, zpar, epar, afMult, czMult, levelDensityType)

! ===================================================================================
!
! Interface to set level density parameters for the user if they wish to use class
! specific data (highly recommended)
!
! ===================================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64

    implicit none
    real(real64),   intent(in   ) :: apar
    real(real64),   intent(in   ) :: zpar
    real(real64),   intent(in   ) :: epar
    real(real64),   intent(  out) :: afMult
    real(real64),   intent(  out) :: czMult
    integer(int32), intent(in   ), optional :: levelDensityType

! ===================================================================================

    ! NOTE:
    ! afMult  = a_f(Nucleon)/a_f(RAL)
    ! czMult = C(Z)[Nucleon]/C(Z)[RAL]
    if ( present(levelDensityType) ) then
       if ( levelDensityType <= 0 ) then
          ! Using level density parameters for nucleon-induced reactions
          call setLevelDensityNucleon(apar, zpar, epar, afMult, czMult)
       else
          ! Use level density parameters for light-, heavy-ion induced reactions or for very high energies
          ! These are the LAQGSM values [old, needs updating]
          call setLevelDensityIon(apar, zpar, epar, afMult, czMult)
       end if
    else
       ! Using level density parameters for nucleon-induced reactions
       call setLevelDensityNucleon(apar, zpar, epar, afMult, czMult)
    end if


    return
! ===================================================================================
  end subroutine setLevelDensity


  subroutine setLevelDensityNucleon(apar, zpar, epar, afMult, czMult)

! ===================================================================================
!
! Sets level density parameters for nucleon induced reactions at smaller energies
!
! ===================================================================================

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none
    real(real64), intent(in   ) :: apar
    real(real64), intent(in   ) :: zpar
    real(real64), intent(in   ) :: epar
    real(real64), intent(  out) :: afMult
    real(real64), intent(  out) :: czMult

    real(real64), parameter :: zMin       = 67_real64
    real(real64), parameter :: zThreshold = 88_real64

! ===================================================================================

    if (zpar > zThreshold ) then 
       call fitafac (apar, zpar, epar, afMult, czMult)
    elseif (zpar <= zThreshold .and. zpar >= zMin ) then
       call fitafpa (apar, zpar, epar, afMult, czMult)
    else
       ! No fission will occur
       afMult = one
       czMult = one
    endif

    return
! ===================================================================================
  end subroutine setLevelDensityNucleon


  subroutine setLevelDensityIon(apar, zpar, epar, afMult, czMult)

! ===================================================================================
!
! Sets level density parameters for light-/heavy- ion induced reactions or for 
! highly energetic reactions
!
! ===================================================================================

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none
    real(real64), intent(in   ) :: apar
    real(real64), intent(in   ) :: zpar
    real(real64), intent(in   ) :: epar
    real(real64), intent(  out) :: afMult
    real(real64), intent(  out) :: czMult

    real(real64), parameter :: zMin       = 67_real64
    real(real64), parameter :: zThreshold = 88_real64

! ===================================================================================

    if (zpar > zThreshold ) then
       call fitafacq (apar, zpar, epar, afMult, czMult)
    elseif (zpar <= zThreshold .and. zpar >= zMin ) then
       call fitafpaq (apar, zpar, epar, afMult, czMult)
    else
       ! No fission will occur
       afMult = one
       czMult = one
    endif

    return
! ===================================================================================
  end subroutine setLevelDensityIon
