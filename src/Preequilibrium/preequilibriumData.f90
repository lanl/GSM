
! ======================================================================
!
! Module containing the preequilibrium data
!
!
! Written by CMJ, XCP-3, 8/2018
!
! ======================================================================

module preequilibriumData

  use, intrinsic:: iso_fortran_env, only: int32, real64
  use preequilibriumParams, only: ato3rd

  implicit none
  private


  ! Constants
  logical,        public, protected :: preeqDataDeclared    = .false.       != Flag to determine if data is established
  real(real64),   public, parameter :: t0logicSwitch        = 0.210_real64  != For particles w/ T0 < 210 [MeV], use old logic
  integer(int32), public, parameter :: numAllowedFragments  = 66_int32      ! Number of allowed/considered fragments during Preequilibrium simulation


  ! Data
  integer(int32), public, parameter :: alj0D1Limit = 60_int32
  integer(int32), public, parameter :: alj0D2Limit = 90_int32
  ! Note: the alj0 array uses ~170 kB memory
  !> \brief Exciton phase space factors (size of [alj0D1Limit, alj0D2Limit, 4])
  real(real64),   public, protected, dimension(:, :, :), allocatable :: alj0
  ! Note: the following use ~131 kB memory
  !> \brief Reduced mass (size of [300, 28])
  real(real64),   public, protected, dimension(:, :), allocatable :: emured
  !> \brief Reduced mass relative to A_j (size of [300, 28])
  real(real64),   public, protected, dimension(:, :), allocatable :: emuredf


  ! Include data regarding the 28 allowed fragments (aj, zj, dlmn, ajthr)
  include "cameronData.f90"


  public :: initPreequilibriumData


contains


! ======================================================================

  subroutine initPreequilibriumData()

! ======================================================================
!
! Sets and establishes data for the preequilibrium model
!
!
! Written by CMJ, XCP-3, 8/2018 (Preequil. class creation)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use preequilibriumParams, only: zro, two, six, twelve

    implicit none

    integer(int32) :: ip, l, lmax, k
    real(real64)   :: aaj, aal, aajd

! ======================================================================

    if (preeqDataDeclared) return

    ! Establish exciton phase space factors [used when establishing the "alj" array]
    if (.not.allocated(alj0)) allocate(alj0(alj0D1Limit, alj0D2Limit, 4))
    do l = 2, alj0D2Limit
       lmax = min(l, alj0D1Limit)
       do ip = 1,lmax
          alj0(ip, l, 1) = dble(ip*(l-1))
          alj0(ip, l, 2) = dble((ip-1)*(l-2))*alj0(ip, l, 1)/two
          alj0(ip, l, 3) = dble((ip-2)*(l-3))*alj0(ip, l, 2)/six
          alj0(ip, l, 4) = dble((ip-3)*(l-4))*alj0(ip, l, 3)/twelve
       end do
    end do


    ! Establish reduced mass and reduced mass ratios (neglecting mass defect due to speed)
    if (.not.allocated(emured))  allocate(emured(300, 28))
    if (.not.allocated(emuredf)) allocate(emuredf(300, 28))
    do k = 1,300
       aaj = dble(k)   !  Mass number of parent:
       do l = 1, numAllowedFragments
          aal = aj(l)          !  Mass number of evaporated particle:
          aajd = aaj - aal     !  Mass number of daughter:
          if (aajd > zro) then
             ! Reduced mass /= 0 for aajd > 0
             emuredf(  k, nint( aj(l) )  ) = aajd/aaj       ! Reduced mass fraction
             emured (  k, nint( aj(l) )  ) = aajd*aal/aaj   ! Reduced mass
          else
             ! Reaction is NOT allowed physically; unphysical daughter nucleus (A <= 0)
             ! NOTE: This only occurs for Residual nuclei w/ A <= aj(max) [for ALL A>aj(max), no issue]
          endif
       end do
    end do

    ! Data for preequilibrium class has been declared
    preeqDataDeclared = .true.

    return
! ======================================================================
  end subroutine initPreequilibriumData


end module preequilibriumData
