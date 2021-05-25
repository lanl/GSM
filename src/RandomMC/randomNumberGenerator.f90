! Copyright LANS/LANL/DOE - see file COPYRIGHT_INFO

module randomNumberGenerator
  !=============================================================================
  ! Description:
  !  randomNumberGenerator.F90 -- random number generation routines
  !  NOTE: Copied from MCNP6.2 random number generator
  !=============================================================================
  !  This module contains:
  !
  !   * Constants for the RN generator, including initial RN seed for the
  !     problem & the current RN seed
  !
  !   * GSM interface routines:
  !     - random number function:           rang()
  !     - RN initialization for problem:    RN_init_problem
  !     - RN initialization for particle:   RN_init_particle
  !     - RN init for particle, special:    RN_next_particle
  !     - get info on RN parameters:        RN_query
  !     - get RN seed for n-th history:     RN_query_first
  !     - set new RN parameters:            RN_set
  !     - skip-ahead in the RN sequence:    RN_skip_ahead
  !     - Unit tests:        RN_test_basic, RN_test_skip, RN_test_mixed
  !     - not for general use:              RN_set_seed
  !
  !   * For interfacing with the rest of GSM, arguments to/from these
  !     routines will have types of I8 or I4.
  !     Any args which are to hold random seeds, multipliers,
  !     skip-distance will be type I8, so that 63 bits can be held without
  !     truncation.
  !
  ! Revisions:
  ! * 10-04-2001 - F Brown, initial MCNP version
  ! * 06-06-2002 - F Brown, mods for extended generators
  ! * 12-21-2004 - F Brown, added 3 of LeCuyer's 63-bit mult. RNGs
  ! * 01-29-2005 - J Sweezy, Modify to use MCNP modules prior to automatic
  !                io unit numbers.
  ! * 12-02-2005 - F Brown, mods for consistency with C version
  ! * 12-12-2006 - C Zeeb, added subroutine RN_next_particle
  ! * 09-25-2012 - F Brown, changed conversion of integer*8 seed to
  !                real*8 fraction, to prevent roundoff to exactly 1.0
  !                for the largest 512 integer*8 seeds when using 63-bit
  !                generators. Also changes the smallest RN to 2**-53.
  !                Always use RN_NORM=2**-53, shift seed by 53-RN_BITS
  !                in creating fraction.
  ! * 05-24-2018 - F Brown, trivial mods to allow use with MCNP6 or stand-alone
  ! * 08-22-2018 - C Juneau, modification for GSM High Energy Physics code
  !=============================================================================
  use ISO_FORTRAN_ENV, only: real64,  int32, int64,  OUTPUT_UNIT,  ERROR_UNIT
  
  implicit none
  private


  ! For filtered message printing:
  integer(int32), private, parameter :: defaultRNGVerbose = 4_int32
  integer(int32), public :: rngVerbose = defaultRNGVerbose

  ! To handle all printing, be default:
  abstract interface
     subroutine IOHANDLER(verbosity, type, text)
       use, intrinsic:: iso_fortran_env, only: int32
       implicit none
       integer(int32),   intent(in) :: verbosity
       integer(int32),   intent(in) :: type
       character(len=*), intent(in) :: text
     end subroutine IOHANDLER
  end interface
  type, private :: rngIO
     private
     character(LEN=512),   private :: message = ""
     procedure(IOHANDLER), public, nopass, pointer :: print => printRandom
  end type rngIO



  integer(int32), parameter :: iuo  = OUTPUT_UNIT
  integer(int32), parameter :: jtty = OUTPUT_UNIT
!  integer(int32), parameter :: jtty =  ERROR_UNIT

  !-----------------------------------
  ! Public functions and subroutines for this module
  !-----------------------------------
  PUBLIC :: rang
  PUBLIC :: RN_init_problem
  PUBLIC :: RN_init_particle
  PUBLIC :: RN_next_particle
  PUBLIC :: RN_set
  PUBLIC :: RN_query
  PUBLIC :: RN_query_first
  PUBLIC :: RN_update_stats
  PUBLIC :: RN_swap
  PUBLIC :: RN_test_basic
  PUBLIC :: RN_test_skip
  PUBLIC :: RN_test_mixed
  PUBLIC :: RN_set_seed  


  ! Create default parameterizations for a RNG object:
  integer(int32), public,    parameter :: defaultIndex    = 1_int32
  integer(int64), public,    parameter :: defaultMult     = 19073486328125_int64
  integer(int64), public,    parameter :: defaultAddC     = 0_int64
  integer(int32), public,    parameter :: defaultLog2Mod  = 48_int32
  integer(int64), public,    parameter :: defaultStride   = 152917_int64
  integer(int64), public,    parameter :: defaultInitSeed = 19073486328125_int64
  character(len=8), public,  parameter :: defaultName     = " GSM Std"


  !-------------------------------------
  ! Random Number Generator Object:
  !-------------------------------------
  type, PUBLIC :: RN_GEN
     private
     integer(int32),   public :: index    = defaultIndex
     integer(int64),   public :: mult     = defaultMult       ! generator (multiplier)
     integer(int64),   public :: add      = defaultAddC       ! additive constant
     integer(int32),   public :: log2mod  = defaultLog2Mod    ! log2 of modulus, must be <64
     integer(int64),   public :: stride   = defaultStride     ! stride for particle skip-ahead
     integer(int64),   public :: initseed = defaultInitSeed   ! default seed for problem
     character(len=8), public :: name     = defaultName       ! Default name for base RNG
  end type RN_GEN


  ! Create several default RNGs:
  integer(int32),      parameter, public :: numRNGenerators = 7
  type(RN_GEN),        parameter, public :: standard_generator(numRNGenerators) =  [ &
       & RN_GEN(), &   ! Creates RNG of default value
       & RN_GEN( 2, 9219741426499971445_int64, 1_int64, 63, 152917_int64, 1_int64, 'LEcuyer1' ), &
       & RN_GEN( 3, 2806196910506780709_int64, 1_int64, 63, 152917_int64, 1_int64, 'LEcuyer2' ), &
       & RN_GEN( 4, 3249286849523012805_int64, 1_int64, 63, 152917_int64, 1_int64, 'LEcuyer3' ), &
       & RN_GEN( 5, 3512401965023503517_int64, 0_int64, 63, 152917_int64, 1_int64, 'LEcuyer4' ), &
       & RN_GEN( 6, 2444805353187672469_int64, 0_int64, 63, 152917_int64, 1_int64, 'LEcuyer5' ), &
       & RN_GEN( 7, 1987591058829310733_int64, 0_int64, 63, 152917_int64, 1_int64, 'LEcuyer6' )  &
       & ]

  ! random generator actually used in the problem
  type(rngIO), public :: io
  type(RN_GEN), target, PUBLIC :: problem_generator

  !-----------------------------------------------------------------
  !   * Linear multiplicative congruential RN algorithm:
  !
  !            RN_SEED = RN_SEED*RN_MULT + RN_ADD  mod RN_MOD
  !
  !   * There are NO DEFAULT VALUES. 
  !     Must call RN_init_problem before using
  !-----------------------------------------------------------------
  integer(int32) :: RN_INDEX     !  = 1
  integer(int32) :: RN_BITS      !  = 48
  integer(int64) :: RN_MULT      !  =  19073486328125_int64
  integer(int64) :: RN_ADD       !  =               0_int64
  integer(int64) :: RN_STRIDE    !  =          152917_int64
  integer(int64) :: RN_SEED0     !  =  19073486328125_int64
  integer(int64) :: RN_MOD       !  = 281474976710656_int64
  integer(int64) :: RN_MASK      !  = 281474976710655_int64
  integer(int64) :: RN_PERIOD    !  =  70368744177664_int64
  integer(int64) :: RN_SHIFT     !  = 5

  real(real64), parameter :: RN_NORM = 2._real64**(-53)          ! 53-bits in ieee 
  !------------------------------------
  ! Private data for a single particle
  !   * There are NO DEFAULT VALUES. 
  !     Must call RN_init_problem before using
  !------------------------------------
  integer(int64) :: RN_SEED      ! = 19073486328125_int64 ! current seed
  integer(int64) :: RN_COUNT     ! = 0_int64              ! current counter
  integer(int64) :: RN_NPS       ! = 0_int64              ! current particle number



  !------------------------------------------
  ! Shared data, to collect info on RN usage
  !------------------------------------------
  integer(int64), SAVE :: RN_COUNT_TOTAL   = 0  ! total RN count all particles
  integer(int64), SAVE :: RN_COUNT_STRIDE  = 0  ! count for stride exceeded
  integer(int64), SAVE :: RN_COUNT_PERIOD  = 0  ! count for period exceeded
  integer(int64), SAVE :: RN_COUNT_MAX     = 0  ! max RN count all particles
  integer(int64), SAVE :: RN_COUNT_MAX_NPS = 1  ! part index for max count
  integer(int64), SAVE :: RN_COUNT_ADVANCES= 0  ! Used by RN_next_particle

  !---------------------------------------------------------------------
  ! Reference data:  Seeds for case of init.seed = 1,
  !                  Seed numbers for index 1-5, 123456-123460
  !---------------------------------------------------------------------
  integer(int64), parameter, dimension(10,numRNGenerators) ::  RN_CHECK = reshape(([ &
    ! ***** 1 ***** gsm standard gen *****
    &      19073486328125_int64,      29763723208841_int64,     187205367447973_int64, &
    &     131230026111313_int64,     264374031214925_int64,     260251000190209_int64, &
    &     106001385730621_int64,     232883458246025_int64,      97934850615973_int64, &
    &     163056893025873_int64, &
    ! ***** 2 *****
    & 9219741426499971446_int64,  666764808255707375_int64, 4935109208453540924_int64, &
    & 7076815037777023853_int64, 5594070487082964434_int64, 7069484152921594561_int64, &
    & 8424485724631982902_int64,   19322398608391599_int64, 8639759691969673212_int64, &
    & 8181315819375227437_int64, &
    ! ***** 3 *****
    & 2806196910506780710_int64, 6924308458965941631_int64, 7093833571386932060_int64, &
    & 4133560638274335821_int64,  678653069250352930_int64, 6431942287813238977_int64, &
    & 4489310252323546086_int64, 2001863356968247359_int64,  966581798125502748_int64, &
    & 1984113134431471885_int64, &
    ! ***** 4 *****
    & 3249286849523012806_int64, 4366192626284999775_int64, 4334967208229239068_int64, &
    & 6386614828577350285_int64, 6651454004113087106_int64, 2732760390316414145_int64, &
    & 2067727651689204870_int64, 2707840203503213343_int64, 6009142246302485212_int64, &
    & 6678916955629521741_int64, &
    ! ***** 5 *****
    & 3512401965023503517_int64, 5461769869401032777_int64, 1468184805722937541_int64, &
    & 5160872062372652241_int64, 6637647758174943277_int64,  794206257475890433_int64, &
    & 4662153896835267997_int64, 6075201270501039433_int64,  889694366662031813_int64, &
    & 7299299962545529297_int64, &
    ! ***** 6 *****
    & 2444805353187672469_int64,  316616515307798713_int64, 4805819485453690029_int64, &
    & 7073529708596135345_int64, 3727902566206144773_int64, 1142015043749161729_int64, &
    & 8632479219692570773_int64, 2795453530630165433_int64, 5678973088636679085_int64, &
    & 3491041423396061361_int64, &
    ! ***** 7 *****
    & 1987591058829310733_int64, 5032889449041854121_int64, 4423612208294109589_int64, &
    & 3020985922691845009_int64, 5159892747138367837_int64, 8387642107983542529_int64, &
    & 8488178996095934477_int64,  708540881389133737_int64, 3643160883363532437_int64, &
    & 4752976516470772881_int64  ]), shape(RN_CHECK) )
  !---------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------

  function rang()   BIND(C)
    ! GSM random number generator
    !
    ! ***************************************
    ! ***** modifies RN_SEED & RN_COUNT *****
    ! ***************************************
    use, intrinsic:: iso_fortran_env, only: int64
    use, intrinsic:: iso_C_binding, only: c_double
    implicit none
    real(c_double) ::  rang

    RN_SEED  = iand( iand( RN_MULT*RN_SEED, RN_MASK) + RN_ADD,  RN_MASK)
    rang     = max( ishft(RN_SEED,RN_SHIFT), 1_int64 )  * RN_NORM
    RN_COUNT = RN_COUNT + 1

    return
  end function rang

  !-------------------------------------------------------------------

  function RN_skip_ahead( seed, skip )
    ! advance the seed "skip" RNs:   seed*RN_MULT^n mod RN_MOD
    implicit none
    integer(int64) :: RN_skip_ahead
    integer(int64), intent(in)  :: seed, skip
    integer(int64) :: nskip, gen, g, inc, c, gp, rn, seed_old

    seed_old = seed
    ! add period till nskip>0
    nskip = skip
    do while( nskip<0_int64 )
      if( RN_PERIOD>0_int64 ) then
        nskip = nskip + RN_PERIOD
      else
        nskip = nskip + RN_MASK
        nskip = nskip + 1_int64
      endif
    enddo

    ! get gen=RN_MULT^n,  in log2(n) ops, not n ops !
    nskip = iand( nskip, RN_MASK )
    gen   = 1
    g     = RN_MULT
    inc   = 0
    c     = RN_ADD
    do while( nskip>0_int64 )
      if( btest(nskip,0) )  then
        gen = iand( gen*g, RN_MASK )
        inc = iand( inc*g, RN_MASK )
        inc = iand( inc+c, RN_MASK )
      endif
      gp    = iand( g+1,  RN_MASK )
      g     = iand( g*g,  RN_MASK )
      c     = iand( gp*c, RN_MASK )
      nskip = ishft( nskip, -1 )
    enddo
    rn = iand( gen*seed_old, RN_MASK )
    rn = iand( rn + inc, RN_MASK )
    RN_skip_ahead = rn
    return
  end function RN_skip_ahead

  !-------------------------------------------------------------------

  subroutine RN_init_problem( new_standard_gen, new_seed, &
    &                         new_stride, new_part1,  print_info )
    ! * initialize GSM random number parameters for problem,
    !   based on user input.  This routine should be called
    !   only from the main thread, if OMP threading is being used.
    !
    ! * for initial & continue runs, these args should be set:
    !     new_standard_gen - index of built-in standard RN generator,
    !                        from RAND gen=   (or dbcn(14)
    !     new_seed   - from RAND seed=        (or dbcn(1))
    !     output     - logical, print RN seed & mult if true
    !
    !     new_stride - from RAND stride=      (or dbcn(13))
    !     new_part1  - from RAND hist=        (or dbcn(8))
    !
    ! * for continue runs only, these should also be set:
    !     new_count_total   - from "rnr"   at end of previous run
    !     new_count_stride  - from nrnh(1) at end of previous run
    !     new_count_max     - from nrnh(2) at end of previous run
    !     new_count_max_nps - from nrnh(3) at end of previous run
    !
    ! * check on size of long-ints & long-int arithmetic
    ! * check the multiplier
    ! * advance the base seed for the problem
    ! * set the initial particle seed
    ! * initialize the counters for RN stats
    implicit none
    integer(int32), intent(in) :: new_standard_gen
    integer(int64), intent(in) :: new_seed
    integer(int64), intent(in) :: new_stride
    integer(int64), intent(in) :: new_part1
    integer(int32), intent(in) :: print_info
    character(len=20) :: printseed
    integer(int64)       ::  itemp1, itemp2, itemp3, itemp4

    if( new_standard_gen<1 .or. new_standard_gen>numRNGenerators ) then
      call expire( 'RN_init_problem', &
        &  ' ***** ERROR: illegal index for built-in RN generator. Using default.')
    endif
      
    ! set defaults, override if input supplied: seed, mult, stride
    RN_INDEX   = new_standard_gen
    RN_MULT    = standard_generator(RN_INDEX)%mult
    RN_ADD     = standard_generator(RN_INDEX)%add
    RN_STRIDE  = standard_generator(RN_INDEX)%stride
    RN_SEED0   = standard_generator(RN_INDEX)%initseed
    RN_BITS    = standard_generator(RN_INDEX)%log2mod
    RN_MOD     = ishft( 1_int64,       RN_BITS )
    RN_MASK    = ishft( not(0_int64),  RN_BITS-64 )
    RN_SHIFT   = 53 - RN_BITS

    ! Set period
    if( RN_ADD==0_int64) then
      RN_PERIOD  = ishft( 1_int64, RN_BITS-2 )
    else
      RN_PERIOD  = ishft( 1_int64, RN_BITS )
    endif

    ! Set seed (if valid one is supplied)
    if( new_seed>0_int64 ) then
      RN_SEED0  = new_seed
    endif

    ! Set stride for each event
    if( new_stride>0_int64 ) then
      RN_STRIDE = new_stride
    endif

    ! Set initial information for RNG
    RN_COUNT_TOTAL   = 0
    RN_COUNT_STRIDE  = 0
    RN_COUNT_PERIOD  = 0
    RN_COUNT_MAX     = 0
    RN_COUNT_MAX_NPS = 1
    RN_COUNT_ADVANCES = 0

    ! Print information to user output and terminal
    if( print_info /= 0 ) then
      write(printseed,'(i20)') RN_SEED0
      write( iuo,1) RN_INDEX, RN_SEED0, RN_MULT, RN_ADD, RN_BITS, RN_STRIDE
      write(jtty,2) RN_INDEX, trim(adjustl(standard_generator(RN_INDEX)%name)), adjustl(printseed)
1     format( &
        & /,' ***************************************************', &
        & /,' * Random Number Generator  = ',i20,             ' *', &
        & /,' * Random Number Seed       = ',i20,             ' *', &
        & /,' * Random Number Multiplier = ',i20,             ' *', &
        & /,' * Random Number Adder      = ',i20,             ' *', &
        & /,' * Random Number Bits Used  = ',i20,             ' *', &
        & /,' * Random Number Stride     = ',i20,             ' *', &
        & /,' ***************************************************',/)
2     format(' Using random number generator ',i2, &
        &    ', named "', A, '", initial seed = ',a20)
    endif

    ! double-check on number of bits in a long int
    if( bit_size(RN_SEED)<64 ) then
      call expire( 'RN_init_problem', &
        &  ' ***** ERROR: <64 bits in long-int, can-t generate RN-s')
    endif
    itemp1 = 5_int64**25
    itemp2 = 5_int64**19
    itemp3 = ishft(2_int64**62-1_int64,1) + 1_int64
    itemp4 = itemp1*itemp2
    if( iand(itemp4,itemp3)/=8443747864978395601_int64 ) then
      call expire( 'RN_init_problem', &
        &  ' ***** ERROR: can-t do 64-bit integer ops for RN-s')
    endif

    if( new_part1>1_int64 ) then
      ! advance the problem seed to what it would be [(new_part-1)*stride] random numbers from now, for part1
      RN_SEED0 = RN_skip_ahead( RN_SEED0, (new_part1-1_int64)*RN_STRIDE )
      itemp1   = RN_skip_ahead( RN_SEED0, RN_STRIDE )
      if( print_info /= 0 ) then
        write(printseed,'(i20)') itemp1
        write( iuo,3) new_part1,  RN_SEED0, itemp1
        write(jtty,4) new_part1,  adjustl(printseed)
3       format( &
          & /,' ***************************************************', &
          & /,' * Random Number Seed will be advanced to that for *', &
          & /,' * previous particle number = ',i20,             ' *', &
          & /,' * New RN Seed for problem  = ',i20,             ' *', &
          & /,' * Next Random Number Seed  = ',i20,             ' *', &
          & /,' ***************************************************',/)
4       format(' comment. advancing random number to particle ',i12, &
          &    ', initial seed = ',a20)
      endif
    endif

    ! assign the data structure for the problem random generator
    problem_generator = RN_GEN( RN_INDEX, RN_MULT, RN_ADD, RN_BITS, RN_STRIDE, &
 &                              RN_SEED0, standard_generator( RN_INDEX )%name )

    continue 
    ! set the initial particle seed
    RN_SEED  = RN_SEED0
    RN_COUNT = 0
    RN_NPS   = 0

    return
  end subroutine RN_init_problem

  !-------------------------------------------------------------------

  subroutine RN_init_particle( nps )
    ! initialize GSM random number parameters for event "nps"
    !
    !     * generate a new particle seed from the base seed
    !       & particle index
    !     * set the RN count to zero
    implicit none
    integer(int64), intent(in) :: nps

    RN_SEED  = RN_skip_ahead( RN_SEED0, nps*RN_STRIDE )
    RN_COUNT = 0
    RN_NPS   = nps

    return
  end subroutine RN_init_particle

  !-------------------------------------------------------------------

  subroutine RN_next_particle( nps, skip, np_run )
    ! advance the GSM random number parameters to the next event
    !
    !     * generate a new particle seed from the base seed
    !       & particle index
    !     * set the RN count to zero
    implicit none
    integer(int64), intent(in) :: nps
    integer(int64), intent(in) :: skip
    integer(int64), intent(in) :: np_run

    ! Check for stride being exceeded:
    if ( RN_COUNT > RN_STRIDE ) then
       RN_COUNT_STRIDE = RN_COUNT_STRIDE + int( RN_COUNT / RN_STRIDE )
       write(io%message, 1000) np_run, int( RN_COUNT / RN_STRIDE )
       call io%print(3, 3, io%message)
    end if

    !$OMP CRITICAL (RN_NEXT_PART)
    !>>> RN_COUTN_ADVANCES and RN_SEED need to be threadprivate!
    RN_COUNT_ADVANCES = np_run + skip
    RN_SEED  = RN_skip_ahead( RN_SEED0, RN_COUNT_ADVANCES*RN_STRIDE )
    !$OMP END CRITICAL (RN_NEXT_PART)
    RN_COUNT_TOTAL = RN_COUNT_TOTAL + RN_COUNT
    RN_COUNT = 0
    RN_NPS   = nps
    return
1000 format("Event ", i9, " exceeded the RN stride ", i3, " time(s).")
! 2000 format("The RN period has been exceeded ", i5, " time(s) by the ", &
!           & " end of event ", i9, ".")
  end subroutine RN_next_particle

  !-------------------------------------------------------------------

  subroutine RN_set(  key,  newValue )
    implicit none
    character(len=*), intent(in) :: key
    integer(int64),   intent(in) :: newValue
    character(len=20) :: printseed
    integer(int64) :: itemp1

    if( key == "stride"        ) then
      if( newValue>0_int64 ) then
        RN_STRIDE        = newValue
      endif
    endif
    if( key == "count_total"   )  RN_COUNT_TOTAL   = newValue
    if( key == "count_stride"  )  RN_COUNT_STRIDE  = newValue
    if( key == "count_period"  )  RN_COUNT_PERIOD  = newValue
    if( key == "count_max"     )  RN_COUNT_MAX     = newValue
    if( key == "count_max_nps" )  RN_COUNT_MAX_NPS = newValue
    if( key == "rn_srcnum"     )  RN_NPS           = newValue
    if( key == "seed"          )  then
      if( newValue>0_int64 ) then
        RN_SEED0 = newValue
        RN_SEED  = RN_SEED0
        RN_COUNT = 0
        RN_NPS   = 0
      endif
    endif
    if( key == "part1" ) then
      if( newValue>1_int64 ) then
        ! advance the problem seed to that for part1
        RN_SEED0 = RN_skip_ahead( RN_SEED0, (newValue-1_int64)*RN_STRIDE )
        itemp1   = RN_skip_ahead( RN_SEED0, RN_STRIDE )
        write(printseed,'(i20)') itemp1
        write( iuo,3) newValue,  RN_SEED0, itemp1
        write(jtty,4) newValue,  adjustl(printseed)
3       format( &
          & /,' ***************************************************', &
          & /,' * Random Number Seed will be advanced to that for *', &
          & /,' * previous particle number = ',i20,             ' *', &
          & /,' * New RN Seed for problem  = ',i20,             ' *', &
          & /,' * Next Random Number Seed  = ',i20,             ' *', &
          & /,' ***************************************************',/)
4       format(' comment. advancing random number to particle ',i12, &
          &    ', initial seed = ',a20)
        RN_SEED  = RN_SEED0
        RN_COUNT = 0
        RN_NPS   = 0
      endif
    endif
    return
  end subroutine RN_set

  !-------------------------------------------------------------------

  function RN_query( key )
    implicit none
    integer(int64)                  :: RN_query
    character(len=*), intent(in) :: key
    RN_query = 0_int64
    if( key == "seed"           )  RN_query = RN_SEED
    if( key == "stride"         )  RN_query = RN_STRIDE
    if( key == "mult"           )  RN_query = RN_MULT
    if( key == "add"            )  RN_query = RN_ADD
    if( key == "count"          )  RN_query = RN_COUNT
    if( key == "period"         )  RN_query = RN_PERIOD
    if( key == "count_total"    )  RN_query = RN_COUNT_TOTAL
    if( key == "count_stride"   )  RN_query = RN_COUNT_STRIDE
    if( key == "count_period"   )  RN_query = RN_COUNT_PERIOD
    if( key == "count_max"      )  RN_query = RN_COUNT_MAX
    if( key == "count_max_nps"  )  RN_query = RN_COUNT_MAX_NPS
    if( key == "count_advances" )  RN_query = RN_COUNT_ADVANCES
    if( key == "first"          )  RN_query = RN_SEED0
    return
  end function RN_query

  !-------------------------------------------------------------------

  subroutine RN_set_seed( ISEED )
    implicit none
    integer(int64), intent(in) :: ISEED
    RN_SEED = ISEED
    return
  end subroutine RN_set_seed

  !-------------------------------------------------------------------

  function RN_query_first( nps )
    implicit none
    integer(int64)                  :: RN_query_first
    integer(int64),      intent(in) :: nps
    RN_query_first = RN_skip_ahead( RN_SEED0, nps*RN_STRIDE )
    return
  end function RN_query_first

  !-------------------------------------------------------------------

  subroutine RN_update_stats()
    ! update overall RN count info
    implicit none

    !$OMP CRITICAL (RN_STATS)

    RN_COUNT_TOTAL = RN_COUNT_TOTAL + RN_COUNT

    if( RN_COUNT>RN_COUNT_MAX ) then
      RN_COUNT_MAX     = RN_COUNT
      RN_COUNT_MAX_NPS = RN_NPS
    elseif ( RN_COUNT==RN_COUNT_MAX) then
      RN_COUNT_MAX_NPS = min(RN_NPS, RN_COUNT_MAX_NPS)
    endif

    if( RN_COUNT>RN_STRIDE ) then
      RN_COUNT_STRIDE = RN_COUNT_STRIDE + 1
    endif

    !$OMP END CRITICAL (RN_STATS)

    RN_COUNT = 0
    RN_NPS   = 0

    return
  end subroutine RN_update_stats

  !-------------------------------------------------------------------

  subroutine RN_swap( RNG, SEED, CNT )
    implicit none

    type(RN_GEN), intent(in) :: RNG
    integer(int64), intent(inout)  :: SEED, CNT

    integer(int64) :: swap

    ! set the parameters to the new random number generator
    RN_INDEX  = RNG%index
    RN_MULT   = RNG%mult
    RN_ADD    = RNG%add
    RN_BITS   = RNG%log2mod
    RN_STRIDE = RNG%stride
    RN_SEED0  = RNG%initseed

    RN_MOD     = ishft( 1_int64,       RN_BITS )
    RN_MASK    = ishft( not(0_int64),  RN_BITS-64 )
    RN_SHIFT   = 53 - RN_BITS

    if( RN_ADD==0_int64) then
      RN_PERIOD  = ishft( 1_int64, RN_BITS-2 )
    else
      RN_PERIOD  = ishft( 1_int64, RN_BITS )
    endif

    ! swap the seed and count (nps must stay the same in swaps)
    swap    = RN_SEED
    RN_SEED = SEED
    SEED    = swap

    swap     = RN_COUNT
    RN_COUNT = CNT
    CNT      = swap

  end subroutine RN_swap

  !-------------------------------------------------------------------
  subroutine expire( c1,c2 )
    implicit none
    character(len=*), intent(in) :: c1, c2
    write(*,*) ' ********** error: ',trim(c1)
    write(*,*) ' ********** error: ',trim(c2)
    stop       ' ********** error **********'
  end subroutine expire
  !-------------------------------------------------------------------
  !###################################################################
  !#
  !#  Unit tests
  !#
  !###################################################################

  subroutine RN_test_basic( new_gen )
    ! test routine for basic random number generator
    implicit none
    integer(int32), intent(in) :: new_gen
    real(real64)    :: s
    integer(int64)  :: seeds(10)
    integer(int32)  :: i, j

    write(jtty,"(/,a)")  " ***** random number - basic test *****"

    ! set the seed
    call RN_init_problem( new_gen, 1_int64, 0_int64, 0_int64, 0 )

    ! get the first 5 seeds, then skip a few, get 5 more - directly
    s = 0.0_real64
    do  i = 1,5
      s = s + rang()
      seeds(i) = RN_query( "seed" )
    enddo
    do  i = 6,123455
      s = s + rang()
    enddo
    do  i = 6,10
      s = s + rang()
      seeds(i) = RN_query( "seed" )
    enddo

    ! compare
    do  i = 1,10
      j = i
      if( i>5  ) j = i + 123450
      write(jtty,"(1x,i6,a,i20,a,i20)") &
        &  j, "  reference: ", RN_CHECK(i,new_gen), "  computed: ", seeds(i)
      if( seeds(i)/=RN_CHECK(i,new_gen) ) then
        write(jtty,"(a)")  " ***** basic_test of RN generator failed:"
      endif
    enddo
    return
  end subroutine RN_test_basic

  !-------------------------------------------------------------------

  subroutine RN_test_skip( new_gen )
    ! test routine for basic random number generation & skip-ahead
    implicit none
    integer(int32), intent(in) :: new_gen
    integer(int64) :: seeds(10)
    integer(int32) :: i, j

    ! set the seed
    call RN_init_problem( new_gen, 1_int64, 0_int64, 0_int64, 0 )

    ! use the skip-ahead function to get first 5 seeds, then 5 more
    do i = 1,10
      j = i
      if( i>5 )  j = i + 123450
      seeds(i) = RN_skip_ahead( 1_int64, int(j,int64) )
    enddo

    ! compare
    write(jtty,"(/,a)")  " ***** random number - skip test *****"
    do i = 1,10
      j = i
      if( i>5  ) j = i + 123450
      write(jtty,"(1x,i6,a,i20,a,i20)") &
        &  j, "  reference: ", RN_CHECK(i,new_gen),  "  computed: ", seeds(i)
      if( seeds(i)/=RN_CHECK(i,new_gen) ) then
        write(jtty,"(a)")  " ***** skip_test of RN generator failed:"
      endif
    enddo
    return
  end subroutine RN_test_skip

  !-------------------------------------------------------------------

  subroutine RN_test_mixed( new_gen )
    ! test routine -- print RN's 1-5 & 123456-123460,
    !                 with reference vals
    implicit none
    integer(int32), intent(in) :: new_gen
    integer(int64) :: r
    integer(int32) :: i, j

    write(jtty,"(/,a)")  " ***** random number - mixed test *****"
    ! set the seed & set the stride to 1
    call RN_init_problem( new_gen, 1_int64, 1_int64, 0_int64, 0 )

    write(jtty,"(a,i20,z20)") " RN_MULT   = ", RN_MULT, RN_MULT
    write(jtty,"(a,i20,z20)") " RN_ADD    = ", RN_ADD,  RN_ADD
    write(jtty,"(a,i20,z20)") " RN_MOD    = ", RN_MOD,  RN_MOD
    write(jtty,"(a,i20,z20)") " RN_MASK   = ", RN_MASK, RN_MASK
    write(jtty,"(a,i20)")     " RN_BITS   = ", RN_BITS
    write(jtty,"(a,i20)")     " RN_PERIOD = ", RN_PERIOD
    write(jtty,"(a,es21.13)") " RN_NORM   = ", RN_NORM
    write(jtty,"(a)")  " "
    do i = 1,10
      j = i
      if( i>5  ) j = i + 123450
      call RN_init_particle( int(j,int64) )
      r = RN_query( "seed" )
      write(jtty,"(1x,i6,a,i20,a,i20)") &
        &  j, "  reference: ", RN_CHECK(i,new_gen),"  computed: ", r
      if( r/=RN_CHECK(i,new_gen) ) then
        write(jtty,"(a)")  " ***** mixed test of RN generator failed:"
      endif
    enddo
    return
  end subroutine RN_test_mixed

  include "printRandom.f90"

  !-------------------------------------------------------------------
end module randomNumberGenerator
