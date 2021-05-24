
! ======================================================================
!
! This module contains data that is used for the Fermi Breakup (FBU) library.
! The data is establish via the "setFermiData" routine, which the FBU library utilizes
! through its "fermi_init" routine, which simply calls this routine.
!
! Created by CMJ, XCP-3, July 2018
!
! ======================================================================
module fermiBreakUpData

  ! Parameters used by fermi break up routine
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use fermiBreakupParams, only : zro, one, two, twelve, thrd, pi

  implicit none
  private

  ! Flags if class data was initialized:
  logical,        public, protected :: fermiDataInitialized = .FALSE.


  ! For array sizes:
  ! [For storage of produced fragments (data type)]
  integer(int32), public, parameter :: min_fermi_AData = nint( twelve )
  ! [For calculation (can be allocated when FBU used)]
  integer(int32), public, parameter :: mp      = 8000 ! Size for iaz array
  integer(int32), public, parameter :: mf      = 1999 ! Used to ensure within array limits
  integer(int32), public, parameter :: mf_wnq  = 2000 ! Size for wnq array
  ! [used in 'crack' and 'divaz' routines for mpa/mpz/mra/mrz arrays]
  integer(int32), public, parameter :: mpraz   =  200 ! Size for mpa/mpz/mra/mrz arrays


  ! Model data
  integer(int32), private, parameter :: defaultGafSize   = 30_int32
  integer(int32), public,  parameter :: gafSize          = max(defaultGafSize, min_fermi_AData)
  real(real64),   public,  protected, dimension(min_fermi_AData, min_fermi_AData-1), save    :: ms = zro
  real(real64),   public,  protected, dimension(min_fermi_AData, min_fermi_AData-1), save    :: dm = zro
  real(real64),   public,  protected, dimension( gafSize ),                          save    :: gaf = zro
  ! (perhaps vak is the system volume, in energy or mass, of a single particle?)
  ! Could this be: vak = (1/3)Sqrt(2/pi) [5.07 * Sqrt(2 E_0) ]**3 the actual equation (not seen in literature)??
  !    NOTE: The sDCM uses a similar factor, 5.06768... to convert angular momentum to units of (h-bar c).
  !          It seems there may be some parallels here from statistical physics
  real(real64),   public,  parameter :: vak              = &
       & thrd * sqrt( two/pi ) * ( (1.4d0*sqrt(0.94d0)*5.07d0)**3 )



  public :: initializeFermiBreakUpData

contains

  subroutine initializeFermiBreakUpData

! ======================================================================
!
! Sets up the ms/dm/gaf arrays
!
!    Last change: 13-Aug-2003 BY NVMokhov
!    Modified by A. J. Sierk, LANL T-16, September-November, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fermiBreakupParams, only : one, one_eighty, &
         & thousandth

    implicit none
    integer(int32) :: i, j, k, md, ms1
    real(real64)   :: gq, qk, qki

! ======================================================================

!  2*S + 1  with some empirical multipliers for common isotopes:
    dimension ms1(min_fermi_AData, min_fermi_AData-1)
    data ms1 / &
         & 0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, & ! Column  1
         & 0,  3,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, & ! Column  2
         & 0,  2,  1,  4,  1,  0,  1,  0,  0,  0,  0,  0, & ! Column  3
         & 0,  0,  0,  4,  6,  5,  4,  0,  4,  0,  0,  0, & ! Column  4
         & 0,  0,  0,  6,  1,  4,  7,  2,  1,  0,  1,  0, & ! Column  5
         & 0,  0,  0,  5,  4, 14, 16, 18,  4,  5,  0,  0, & ! Column  6
         & 0,  0,  0,  4,  1,  4, 11, 14, 13,  8,  1,  0, & ! Column  7
         & 0,  0,  0,  0,  0,  3,  2, 18, 10,  5,  2,  3, & ! Column  8
         & 0,  0,  0,  0,  0,  4,  1, 10, 17,  6,  1,  6, & ! Column  9
         & 0,  0,  0,  0,  0,  0,  0,  1,  6,  3,  2,  5, & ! Column 10
         & 0,  0,  0,  0,  0,  0,  0,  2,  1,  2,  1,  4 /   ! Column 11

!  Mass excesses in keV:
    dimension md(min_fermi_AData, min_fermi_AData-1)
    data ((md(i,j),i=1,12), j=1,9)/   &
         99000,  8071, 99999, 99999, 99999, 99999, 99999, 99999, 99999,   &
         99999, 99999, 99999,    &
         99000, 13136, 14950, 25900, 36800, 41900, 99999, 99999, 99999,   &
         99999, 99999, 99999,    &
         99000, 14931,  2425, 11390, 17594, 26110, 31598, 40820, 48810,   &
         65000, 75240, 89260,    &
         99000, 25300, 11680, 14086, 14908, 20946, 24954, 33050, 40800,   &
         52940, 61570, 72280,    &
         99000, 38000, 18375, 15770,  4942, 11348, 12607, 20174, 25080,   &
         33700, 39900, 51210,    &
         99000, 99999, 27870, 22921, 12416, 12051,  8668, 13369, 16562,   &
         23660, 28970, 37080,      &
         99000, 99999, 35090, 28914, 15699, 10651,   000,  3125,  3020,   &
         9873, 13694, 21040,      &
         99000, 99999, 99999, 39700, 25300, 17338,  5345,  2863,   101,   &
         5683,  7870, 13120,    &
         99000, 99999, 99999, 42700, 32050, 23111,  8007,  2855,  4737,   &
         809,   782,  3334/

    data ((md(i,j), i=1,12),j=10,11)/   &
         99000, 99999, 99999, 99999, 39700, 33600, 16800, 10680,  1952,   &
         873,  1487,    17,      &
         99000, 99999, 99999, 99999, 49400, 36400, 23990, 16490,  5307,   &
         1751,  7042,  5732/

! ======================================================================

    ! Flag data as constructed:
    fermiDataInitialized = .TRUE.

    ! -----------------------
    ! Set up 'ms', 'dm' arrays
    ! -----------------------
    ! NOTE:  proton number is j - 1
    ! NOTE:  neutron number is i - 1
    do i = 1,min_fermi_AData
       do j = 1,min_fermi_AData - 1
          ms(i,j) = ms1(i,j)
          dm(i,j) = thousandth*dble(md(i,j))
       end do
    end do
    ms(1,2) = 2
    dm(1,2) = 7.289d0
    dm(9,9) = -dm(9,9)
    dm(10,9) = -dm(10,9)
    dm(11,9) = -dm(11,9)
    dm(11,10) = -dm(11,10)
    dm(12,10) = -dm(12,10)
    dm(11,11) = -dm(11,11)
    dm(12,11) = -dm(12,11)


    ! -----------------------
    ! Set up 'gaf' array
    ! -----------------------
    gaf(1) = one
    gaf(2) = one
    do k = 3, gafSize
       qki = 1.5d0*dble(k-1) - one   ! Exponent for energy in M_{N} equation
       qk = one/qki
       gq = one + qk*(one + qk*(one - qk*139.d0/one_eighty)/24.d0)/12.d0
       gaf(k) = sqrt(0.1591549d0 * qk)/gq
    end do

    return
! ======================================================================
  end subroutine initializeFermiBreakUpData

end module fermiBreakUpData
