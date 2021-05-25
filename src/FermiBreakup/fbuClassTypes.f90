
! ====================================================================
!
! This file contains the various data types utilized by the Fermi
! Break-Up object for its simulations.
! 
! These data types primarily consist of:
!    -I/O [handles all messages and their printing]
!    -Options [handles object behavior and numerics]
!    -Progeny [contains all information tracked for FBU progeny]
!    -Results [contains all results for the FBU object]
!    -Calculation [contains all interim results not needed for keeping
!
!
! Written by CMJ, XCP-3 (3/2019)
!
! ====================================================================

  type, private:: fbuIO
     private
     character(LEN=512),   private :: message = ""
     procedure(IOHANDLER), private, nopass, pointer :: print => printFermi
  end type fbuIO


  type, public:: fermiBreakUpOptions
     private
     integer(int32), public  :: recNumNucleons = defaultRecNumNucleons
     real(real64),   public  :: akpScalingFlag = defaultAkpScalingFlag
  end type fermiBreakUpOptions


  type, public :: fermiBreakUpProgeny
     private
     real(real64),   public  :: numBaryons = 0.0_real64   ! A Number (num. nucleons)
     real(real64),   public  :: numProtons = 0.0_real64   ! Z Number (num. protons)
     real(real64),   public  :: kinEnergy  = 0.0_real64   ! Kinetic energy of fragment [MeV]
     real(real64),   public  :: restMass   = 0.0_real64   ! Rest mass of fragment [MeV/c**2]
     real(real64),   public  :: linearXMom = 0.0_real64   ! X-component of linear momentum [MeV/c]
     real(real64),   public  :: linearYMom = 0.0_real64   ! Y-component of linear momentum [MeV/c]
     real(real64),   public  :: linearZMom = 0.0_real64   ! Z-component of linear momentum [MeV/c]
  end type fermiBreakUpProgeny


  type, public  :: fermiBreakUpNucleus
     private
     real(real64),   public  :: numBaryons = 0.0_real64   ! A Number (num. nucleons)
     real(real64),   public  :: numProtons = 0.0_real64   ! Z Number (num. protons)
     real(real64),   public  :: kinEnergy  = 0.0_real64   ! Kinetic energy of fragment [GeV]
     real(real64),   public  :: linearXMom = 0.0_real64   ! X-component of linear momentum [GeV/c]
     real(real64),   public  :: linearYMom = 0.0_real64   ! Y-component of linear momentum [GeV/c]
     real(real64),   public  :: linearZMom = 0.0_real64   ! Z-component of linear momentum [GeV/c]
  end type fermiBreakUpNucleus



  type, public :: fermiBreakUpResults
     private
     ! Regarding progeny:
     integer(int32), public  :: numProgeny = 0_int32
     integer(int32), public  :: maxProgeny = 0_int32
     type(fermiBreakUpProgeny), public, dimension(:), pointer :: progenyBnk => NULL()

     ! Misc. information:
     logical,        private :: constructed = .FALSE.
     integer(int32), public  :: simState = 0_int32

  end type fermiBreakUpResults


  type, private :: fbuInterimCalc
     private
     real(real64),   private, dimension(mp)    :: iaz = 0.0_real64
     real(real64),   private, dimension(mf_wnq):: wnq = 0.0_real64
  end type fbuInterimCalc


  type, private :: fbuCrackCalculation
     private
     integer(int32), private, dimension(mpraz) :: mpa = 0_int32
     integer(int32), private, dimension(mpraz) :: mpz = 0_int32
     integer(int32), private, dimension(mpraz) :: mra = 0_int32
     integer(int32), private, dimension(mpraz) :: mrz = 0_int32
  end type fbuCrackCalculation
