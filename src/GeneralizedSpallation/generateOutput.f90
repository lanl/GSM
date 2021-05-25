
  subroutine generateOutput( gsmObj, projNucleus, targNucleus, &
       & output, bankSize )

! ====================================================================
!
! This procedure acts as the top-level interface by which driver
! clients call to generate output files
!
! Sets up array sizes to be used and verifies construction state of GSM
!
!
! Written by CMJ, XCP-3 (03/2019)
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32

    implicit none
    class(GSM),           intent(inout) :: gsmObj
    type(GSMProjectile),  intent(inout) :: projNucleus
    type(GSMTarget),      intent(inout) :: targNucleus
    class(GSMOutput),     intent(inout) :: output
    integer(int32),       intent(in   ), optional :: bankSize

    integer(int32) :: maxBankSize = 0_int32

! ====================================================================

    ! Set the max allowed progeny:
    maxBankSize = nint( bankScaling * &
         & (projNucleus%numBaryons + targNucleus%numBaryons) )
    maxBankSize = max( minBankSize, maxBankSize )   ! Require at least 50 progeny (helpful for simulations with only light particles)
    if ( present(bankSize) ) maxBankSize = bankSize

    ! Perform simulation:
    call gsmObj%gsmMain( projNucleus, targNucleus, output, maxBankSize )

    return
! ====================================================================
  end subroutine generateOutput
