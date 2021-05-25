
! ====================================================================
!
! Module file that contains the parameterized data utilized by the
! standard Dubna Cascade Model (DCM)
!
! Written by CMJ, XCP-3, 12/2018
!
! ====================================================================
module standardDCMData

  use, intrinsic:: iso_fortran_env, only: int32, real64
  use photonEventGeneratorClass,    only: PhotonEventGenerator
  use standardDCMParams, only: zro, one, two

  implicit none
  private

  ! Flag stating if data was established
  logical, public, protected :: sDCMDataEstablished = .FALSE.

  ! For message filtering:
  integer(int32), private, parameter :: defaultSDCMVerbose = 4_int32
  integer(int32), public :: sDCMVerbose = defaultSDCMVerbose

  ! Regarding data initialization
  integer(int32), public, parameter :: sDCMDataSetup   =   0_int32   ! Signals no error in data setup
  integer(int32), public, parameter :: sDCMQintError   =   1_int32   ! Signals error from establishing "qints" data
  integer(int32), public, parameter :: sDCMGammaError  =  10_int32   ! Signals error from reading gamma file
  integer(int32), public, parameter :: sDCMDataError   =  -1_int32   ! Signals unknown error



  ! Procedure pointer interface for message handling
  abstract interface
     ! Interface for I/O handling
     subroutine IOHANDLER(verbosity, type, text)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in) :: verbosity
       integer(int32),   intent(in) :: type
       character(len=*), intent(in) :: text
     end subroutine IOHANDLER
  end interface



  ! Import 1st set of data (angular distribution data, momentum
  !    distributions of final state products in pion-production reactions)
  include "sDCMdata1.f90"


  ! Import 2nd set of data (tables of cross sections and energies)
  include "sDCMdata2.f90"


  ! Import 3rd set of data (misc. table and data)
  include "sDCMdata3.f90"


  ! For the photon event generator
  type(PhotonEventGenerator), public :: photonEG

  ! Data for the "qints" function
  ! Note: these arrays require ~20 kB of memory
  real(real64), public, protected, dimension(30, 28) :: aa = zro
  real(real64), public, protected, dimension(30, 28) :: bb = zro
  real(real64), public, protected, dimension(30, 28) :: cc = zro
  real(real64), public, parameter, dimension(    28) :: iiQints = &
       & [ 1, 1, 1, 1, 2, 2, 2, &
       &   2, 2, 2, 3, 3, 3, 3, &
       &   3, 3, 3, 3, 3, 4, 4, &
       &   5, 6, 6, 6, 6, 6, 4  ]


  ! For photon cross section channel data
  ! Data file information (NOTE: This is the same as for mDCM, but with only 4 channels in file)
  ! Note this file is used because it is more likely to open faster.
  character(LEN=*), public,  parameter :: defaultGammaFile = "gamman.tbl"
  ! Data file information (NOTE: Same file as for mDCM)
  ! character(LEN=*), public,  parameter :: defaultGammaFile = "channel1.tab"

  ! Angle bins (course and fine):
  real(real64),   public, parameter :: dtheta  = 10.0_real64   ! Course angle bin width
  real(real64),   public, protected, dimension( 19) ::  theta   = zro   ! Theta value [course]
  real(real64),   public, protected, dimension( 19) :: ctheta   = zro   ! cosine of theta value [course]
  real(real64),   public, protected, dimension( 19) ::  domo2   = zro   ! pi*(bin width by ctheta) [course]
  real(real64),   public, parameter :: dthetai =  1.0_real64   ! Fine   angle bin width
  real(real64),   public, protected, dimension(181) ::  thetai  = zro   ! Theta value [fine]
  real(real64),   public, protected, dimension(181) :: cthetai  = zro   ! cosine of theta value [fine]
  real(real64),   public, protected, dimension(181) ::  domo2i  = zro   ! pi*(bin width by ctheta) [fine]
  ! Data read from photon cross section file
  ! Note: The following arrays require ~18 kB memory when allocated.
  real(real64),   public, protected, dimension(:, :), allocatable:: &
       & gppipn, &   ! gamma + p -> n + pi+ channel data [size=(50,20)]
       & gppi0p, &   ! gamma + p -> p + pi0 channel data [size=(50,20)]
       & elg         ! Energy related to the channel?    [size=( 4,50)]


  ! Pion production threshold
  real(real64),   public, parameter, dimension(3) :: pionProdThresh = &
       & [ 0.70_real64, 0.90_real64, 0.00_real64 ]


  ! For printing in the data module
  procedure(IOHANDLER), private, pointer :: dataIO => printSDCM
  character(LEN=512),   private          :: message = ""


  ! Function for data initialization
  public  :: initializeStandardDCMData
  private :: initializeQintsData   ! Initializes data utilized by the "qints" function
  private :: initializeGammaData   ! Initializes data for photon cross section channels
  private :: printSDCM             ! Handles all printing of messages

contains

! ====================================================================

  ! Include main file for initializing sDCM data
  include "initializeStandardDCMData.f90"


  ! Include "qints" data initialization function
  include "initQintsData.f90"


  ! Include gamma cross section channel data initialization function
  include "initGammaData.f90"


  ! Include function (default) to print out warning messages
  include "printSDCM.f90"


! ====================================================================
end module standardDCMData
