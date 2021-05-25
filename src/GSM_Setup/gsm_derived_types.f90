! Copyright LANS/LANL/DOE - see file COPYRIGHT_INFO

  module gsm_derived_types

    use, intrinsic :: iso_fortran_env, only: real64, int32
    use gsm_params, only: zro


    implicit NONE

    ! For use with parz/spt variable.
    integer(int32), parameter,public :: &
         & bnk_size = 150, &
         & parz_size = bnk_size, &
         & spt_size = bnk_size

    !----------------------------------------------------------- Residual/Compound Nucleus ----------
    type, public :: nucleus
       ! Particle A, Z
       real(real64)   :: ap    = zro     ! Mass Number (A)
       real(real64)   :: zp    = zro     ! Atomic Number (Z)
       ! Particle Energy
       real(real64)   :: up    = zro     ! Kinetic energy [MeV]
       real(real64)   :: ue    = zro     ! Thermal energy [MeV]
       real(real64)   :: trec  = zro     ! Recoil energy  [MeV]
       real(real64)   :: rstMass = zro   ! Rest Mass      [GeV/c**2]
       ! Lineary Momentum
       real(real64)   :: pnx   = zro     ! X-component of linear momentum [GeV/c]
       real(real64)   :: pny   = zro     ! Y-component of linear momentum [GeV/c]
       real(real64)   :: pnz   = zro     ! Z-component of linear momentum [GeV/c]
       ! Angular Momentum
       real(real64)   :: elx   = zro     ! X-component of angular momentum [GeV * "time"]
       real(real64)   :: ely   = zro     ! Y-component of angular momentum [GeV * "time"]
       real(real64)   :: elz   = zro     ! Z-component of angular momentum [GeV * "time"]
       integer(int32) :: ln    = zro     ! Quantum angular momentum number
       ! Fission barrier Height
       real(real64)   :: bf0   = zro     ! Fission barrier [MeV]
       ! 
       integer(int32) :: state = 0_int32 ! State of compound
                                         ! ( excited        = 0,
                                         !   compound       = 1,
                                         !   stable         = 2,
                                         !   fermi broke-up < 0 )
    end type nucleus


  end module gsm_derived_types

