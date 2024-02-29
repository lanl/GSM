
  subroutine initam( )

! ====================================================================
!
! Initializes several common blocks and calls to setup constants/decay table
!
! ====================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64

    implicit none

! ====================================================================
! Simulation varibles (may change w/ simulation)
! enbou
! yesela
! keypla
! qsee, qvsee

! Simulation OPTIONS (static, but can be configured
! keyhh
! inside
! icms
! ivalon
! wtime
! iorhei

! enbou appears to be related to a Coulomb barrier possibly; should
! be made a simulation option and better described.
    real(real64) :: enbou
    common/comenb/ enbou
    integer(int32) :: icms
    common/comfr/  icms
    logical yesela
    common/yesela/ yesela
    logical keyhh
    common/keyhh/ keyhh
    logical keypla
    common/keypla/ keypla
    logical qsee,qvsee
    common/comqse/ qsee,qvsee
    integer(int32) :: inside
    common/cinsid/ inside
    integer(int32) :: ivalon
    common/cvalon/ ivalon
    logical wtime
    common/comwti/ wtime
    integer(int32) :: iorhei
    common/iorhei/ iorhei

! ====================================================================

    yesela=.false.
    keyhh=.false.
    keypla=.false.
    qsee=.false.
    qvsee=.false.
!    WTIME=.TRUE.
    wtime=.false.
    enbou= 4.4
    inside=0
    ivalon=0
    iorhei=1

    ! IF LAB. FRAME
    ICMS=0
    ! IF  C.M. FRAME
    icms=1
    return
! ====================================================================
  end subroutine initam
