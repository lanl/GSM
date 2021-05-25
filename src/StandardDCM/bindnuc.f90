
  function bindnuc (sDCM, ipart, a, z)

! ========================================================================
!
!  Calculates binding energy of a neutron (ipart=1) or proton (ipart=2)
!  in a nucleus with A=a and Z=z (in GeV) using Moller's et al. nuclear
!  masses (1995) or experimental masses, when available
!  Incorrect call to DELTAM, since the ipart argument is not operative.
!  AJS  11/14/03.
!
!    Written by S. G. Mashnik, LANL, T-2, 1998 
!    Modified by SGM in April 2001 to use the new (1999) DELTAM
!    Corrected by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ========================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: one, two, thousandth
    use standardDCMData,   only: dlmn

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    integer(int32), intent(in) :: ipart
    real(real64),   intent(in) :: a
    real(real64),   intent(in) :: z
    real(real64)               :: bindnuc

    real(real64) :: afj, dl, zfj

! ======================================================================

    dl  = sDCM%molnixE%massExcess (a, z)
    afj = a - one
    zfj = z
    if (ipart == two) zfj = z - one
    bindnuc = thousandth * (  sDCM%molnixE%massExcess(afj, zfj) - &
         & (dl - dlmn(ipart) )  )
    return

! ======================================================================
  end function bindnuc
