
! ====================================================================
!
! This module contains all data used by the modified DCM.
!
! ====================================================================

module modifiedDCMData

  use, intrinsic:: iso_fortran_env, only: int32, real64, output_unit
  use modifiedDCMParams, only: one, degreeToRad
  use Contracts
  implicit none
  private

  integer(int32) :: j

  ! For flagging data as initialized or not
  logical, public, protected :: mDCMDataInitialized = .FALSE.

  ! For printing messages:
  integer(int32), private, parameter :: defaultMDCMDataVerbose = 4_int32
  integer(int32), public :: mDCMVerbose = defaultMDCMDataVerbose

  ! For the photon file:
  integer(int32),   private           :: decayUnit        = -1_int32
  character(len=*), public, parameter :: defaultPhotoFile = "channel1.tab"
  character(len=*), public, parameter :: defaultDecayFile = "atab.dat"

  character(len=128), private :: effectiveDecayFile = defaultDecayFile

  ! For checking that a residual is of valid size:
  integer(int32), public, parameter :: numNuclideLimits = 135
  integer(int32), public, parameter, dimension(numNuclideLimits) :: jamin = [ &
       &   1,        2,      3,      4,      5,      6,      7, &
       &   8,        9,     11,     13,     15,     17,     19, &
       &  22,       23,     25,     27,     29,     29,     32, &
       &  31,       36,     35,     40,     37,     43,     41, &
       &  48,       49,     53,     51,     57,     55,     61, &
       &  59,       64,     63,     68,     67,     72,     70, &
       &  76,       74,     80,     78,     84,     82,     88, &
       &  86,       96,     95,    100,     99,    104,    103, &
       & 109,      107,    114,    112,    118,    115,    125, &
       & 120,      130,    126,    134,    131,    140,    136, &
       & 144,      140,    148,    146,    154,    152,    162, &
       & 158,      166,    163,    172,    170,    178,    178, &
       & 184,      182,    189,    188,    195,    192,    200, &
       & 201,      206,    208,    212,    212,    218,    216, &
       & 224,      222,    229,    228,    235,    233,    241, &
       & 240,      247,    240,    253,    256,    259,    262, &
       & 266,      269,    272,    275,    278,    281,    284, &
       & 287,      290,    294,    297,    300,    303,    306, &
       & 310,      313,    316,    319,    323,    326,    329, &
       & 332,      336 ]
  integer(int32), public, parameter, dimension(numNuclideLimits) :: jamax = [ &
       &  12,       14,     17,     22,     25,     28,     31, &
       &  34,       39,     44,     47,     54,     55,     60, &
       &  61,       68,     69,     74,     75,     80,     81, &
       &  82,       85,     88,     91,    100,    103,    112, &
       & 114,      118,    119,    120,    121,    126,    127, &
       & 136,      137,    142,    143,    148,    149,    154, &
       & 157,      162,    163,    168,    169,    172,    175, &
       & 178,      181,    183,    185,    188,    191,    192, &
       & 193,      195,    198,    202,    205,    208,    211, &
       & 214,      218,    221,    224,    227,    230,    234, &
       & 237,      240,    243,    247,    250,    253,    256, &
       & 260,      263,    266,    269,    273,    276,    279, &
       & 282,      286,    289,    292,    295,    299,    302, &
       & 305,      308,    312,    315,    318,    312,    325, &
       & 328,      331,    334,    338,    339,    339,    339, &
       & 339,      339,    339,    339,    339,    339,    339, &
       & 339,      339,    339,    339,    339,    339,    339, &
       & 339,      339,    339,    339,    339,    339,    339, &
       & 339,      339,    339,    339,    339,    339,    339, &
       & 339,      339 ]

  ! For radius of projectile/target w/ [2 < A <= 10]
  ! TODO: Use r_rms module instead?
  real(real64), public, parameter, dimension(10) :: rms = [ &
       & 0.85, 2.095, 1.976, 1.671, 2.50, 2.57, 2.45, 2.519, 2.45, 2.42 ]

  ! For photon cross section information:
  integer(int32), private, parameter :: numThetaBins = 19
  real(real64),   private, parameter :: dtheta = (180.0_real64) / dble(numThetaBins - 1_int32)
  real(real64),   public,  parameter, dimension(numThetaBins) :: &
       & theta  = [ (dble(j-1)*dtheta,          j=1, numThetaBins) ]
  real(real64),   public,  parameter, dimension(numThetaBins) :: &
       & ctheta = [ (cos(theta(j))*degreeToRad, j=1, numThetaBins) ]
  ! For photon data read from the file:
  ! Note these arrays use ~180 kB.
  !> \brief Cross section data (size of [22, 50, 0:18])
  real(real64),   public, protected, dimension(:, :, :), allocatable :: xsectd
  !> \brief Center of mass energy? (size of [22, 50])
  real(real64),   public, protected, dimension(:, :), allocatable :: ecm
  !> \brief elg (size of [22, 50])
  real(real64),   public, protected, dimension(:, :), allocatable :: elg

  ! For decay data read from the file:
  !> \brief Indicates if decay channels will be simulated.
  ! TODO: Place into mDCM data class?
  logical, public :: model_decay = .true.
  !> \brief look decay data (size of 400)
  integer(int32), public, protected, dimension(:), allocatable :: look
  !> \brief cbr decay data (size of 600)
  real(real64),   public, protected, dimension(:), allocatable :: cbr
  !> \brief mode decay data (size of 5, 600)
  integer(int32), public, protected, dimension(:, :), allocatable :: mode


  ! For message handling:
  abstract interface
     subroutine IOHANDLER(msgVerb, msgType, msg)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in   ) :: msgVerb
       integer(int32),   intent(in   ) :: msgType
       character(LEN=*), intent(in   ) :: msg
     end subroutine IOHANDLER
  end interface
  type, private :: mDCMDataIO
     private
     character(LEN=512), private :: message = ""
     procedure(IOHANDLER), private, nopass, pointer :: print => printMDCM
  end type mDCMDataIO
  type(mDCMDataIO), private :: dataIO



  ! Module procedures:
  public  :: initializeModifiedDCMData
  private :: initam
  private :: readPhotonData   ! Reads the effective photon data file
  private :: readDecayData    ! Reads the effective decay data file
  private :: printMDCM

contains

! ====================================================================

#include "initializeModifiedDCMData.f90"
#include "initam.f90"
#include "readPhotonData.f90"
#include "readDecayData.f90"

#include "../printMDCM.f90"

! ====================================================================

end module modifiedDCMData
