
  subroutine setALJ(calcVars, n, ip, p, exn)
!  subroutine setALJ(preeqObj, calcVars, n, ip, p, exn)

! ======================================================================
!
! This routine establishes the alj(:) array in the blalj common block.
! Written to simplifiy peqemt routine.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use preequilibriumParams, only: one, two, thr, four
    use preequilibriumData, only: alj0, alj0D1Limit, alj0D2Limit

    implicit none
!    class(Preequilibrium), intent(inout) :: preeqObj   ! UNUSED!
    type(preequilibriumCalculation), intent(inout) :: calcVars
    integer(int32),        intent(in   ) :: n   ! Used for flags on alj
    integer(int32),        intent(in   ) :: ip  ! Used for flags on alj
    real(real64),          intent(in   ) :: p   ! Used to determine alj(:) values
    real(real64),          intent(in   ) :: exn ! Used to determine alj(:) values

! ======================================================================

    if (n <= alj0D2Limit .and. ip <= alj0D1Limit) then
       ! Use data
       calcVars%alj(1) = alj0(ip,n,1)
       calcVars%alj(2) = calcVars%alj(1)
       calcVars%alj(3) = alj0(ip,n,2)
       calcVars%alj(4) = alj0(ip,n,3)
       calcVars%alj(5) = calcVars%alj(4)
       calcVars%alj(6) = alj0(ip,n,4)
    else
       ! Outside data table
       calcVars%alj(1) = p*(exn - one)
       calcVars%alj(2) = calcVars%alj(1)
       calcVars%alj(3) = calcVars%alj(2)*(p - one)*(exn - two)/two
       calcVars%alj(4) = calcVars%alj(3)*(p - two)*(exn - thr)/6.d0
       calcVars%alj(5) = calcVars%alj(4)
       calcVars%alj(6) = calcVars%alj(5)*(p - thr)*(exn - four)/12.d0
    endif
! LMK 07/2012
    calcVars%alj(7) = calcVars%alj(6)*(p-4.d0)*(p-5.d0)*(exn-5.d0)*(exn-6.d0)/600.d0      ! He-6
    calcVars%alj(8) = calcVars%alj(7)*(p-6.d0)*(p-7.d0)*(exn-7.d0)*(exn-8.d0)/2352.d0     ! He-8
    calcVars%alj(9) = calcVars%alj(7)                                           ! Li-6
    calcVars%alj(10)= calcVars%alj(9)*(p-6.d0)*(exn-7.d0)/42.d0                       ! Li-7
    calcVars%alj(11)= calcVars%alj(8)                                           ! Li-8
    calcVars%alj(12)= calcVars%alj(11)*(p-8.d0)*(exn-9.d0)/72.d0                      ! Li-9
    calcVars%alj(13)= calcVars%alj(10)                                          ! Be-7
    calcVars%alj(14)= calcVars%alj(12)                                          ! Be-9
    calcVars%alj(15)= calcVars%alj(14)*(p-9.d0)*(exn-10.d0)/90.d0                     ! Be-10
    calcVars%alj(16)= calcVars%alj(15)*(p-10.d0)*(exn-11.d0)/110.d0                    ! Be-11
    calcVars%alj(17)= calcVars%alj(16)*(p-11.d0)*(exn-12.d0)/132.d0                   ! Be-12
    calcVars%alj(18)= calcVars%alj(11)                                          ! B-8
    calcVars%alj(19)= calcVars%alj(15)                                          ! B-10
    calcVars%alj(20)= calcVars%alj(16)                                          ! B-11
    calcVars%alj(21)= calcVars%alj(17)                                          ! B-12
    calcVars%alj(22)= calcVars%alj(21)*(p-12.d0)*(exn-13.d0)/156.d0                   ! B-13
    calcVars%alj(23)= calcVars%alj(19)                                          ! C-10
    calcVars%alj(24)= calcVars%alj(20)                                          ! C-11
    calcVars%alj(25)= calcVars%alj(21)                                          ! C-12
    calcVars%alj(26)= calcVars%alj(22)                                          ! C-13
    calcVars%alj(27)= calcVars%alj(26)*(p-13.d0)*(exn-14.d0)/182.d0                   ! C-14
    calcVars%alj(28)= calcVars%alj(27)*(p-14.d0)*(exn-15.d0)/210.d0                   ! C-15
    calcVars%alj(29)= calcVars%alj(28)*(p-15.d0)*(exn-16.d0)/240.d0                   ! C-16
! LMK 12/2012
    calcVars%alj(30)= calcVars%alj(17)                                        ! N-12
    calcVars%alj(31)= calcVars%alj(22)                                        ! N-13
    calcVars%alj(32)= calcVars%alj(27)                                        ! N-14
    calcVars%alj(33)= calcVars%alj(28)                                        ! N-15
    calcVars%alj(34)= calcVars%alj(29)                                        ! N-16
    calcVars%alj(35)= calcVars%alj(34)*(p-16.d0)*(exn-17.d0)/272.d0           ! N-17
    calcVars%alj(36)= calcVars%alj(27)                                        ! O-14
    calcVars%alj(37)= calcVars%alj(28)                                        ! O-15
    calcVars%alj(38)= calcVars%alj(29)                                        ! O-16
    calcVars%alj(39)= calcVars%alj(35)                                        ! O-17
    calcVars%alj(40)= calcVars%alj(39)*(p-17.d0)*(exn-18.d0)/306.d0           ! O-18
    calcVars%alj(41)= calcVars%alj(40)*(p-18.d0)*(exn-19.d0)/342.d0           ! O-19
    calcVars%alj(42)= calcVars%alj(41)*(p-19.d0)*(exn-20.d0)/380.d0           ! O-20
    calcVars%alj(43)= calcVars%alj(35)                                        ! F-17
    calcVars%alj(44)= calcVars%alj(40)                                        ! F-18
    calcVars%alj(45)= calcVars%alj(41)                                        ! F-19
    calcVars%alj(46)= calcVars%alj(42)                                        ! F-20
    calcVars%alj(47)= calcVars%alj(46)*(p-20.d0)*(exn-21.d0)/420.d0           ! F-21
    calcVars%alj(48)= calcVars%alj(40)                                        ! Ne-18
    calcVars%alj(49)= calcVars%alj(41)                                        ! Ne-19
    calcVars%alj(50)= calcVars%alj(42)                                        ! Ne-20
    calcVars%alj(51)= calcVars%alj(47)                                        ! Ne-21
    calcVars%alj(52)= calcVars%alj(51)*(p-21.d0)*(exn-22.d0)/462.d0           ! Ne-22
    calcVars%alj(53)= calcVars%alj(52)*(p-22.d0)*(exn-23.d0)/506.d0           ! Ne-23
    calcVars%alj(54)= calcVars%alj(53)*(p-23.d0)*(exn-24.d0)/552.d0           ! Ne-24
    calcVars%alj(55)= calcVars%alj(47)                                        ! Na-21
    calcVars%alj(56)= calcVars%alj(52)                                        ! Na-22
    calcVars%alj(57)= calcVars%alj(53)                                        ! Na-23
    calcVars%alj(58)= calcVars%alj(54)                                        ! Na-24
    calcVars%alj(59)= calcVars%alj(58)*(p-24.d0)*(exn-25.d0)/600.d0           ! Na-25
    calcVars%alj(60)= calcVars%alj(52)                                        ! Mg-22
    calcVars%alj(61)= calcVars%alj(53)                                        ! Mg-23
    calcVars%alj(62)= calcVars%alj(54)                                        ! Mg-24
    calcVars%alj(63)= calcVars%alj(59)                                        ! Mg-25
    calcVars%alj(64)= calcVars%alj(63)*(p-25.d0)*(exn-26.d0)/650.d0           ! Mg-26
    calcVars%alj(65)= calcVars%alj(64)*(p-26.d0)*(exn-27.d0)/702.d0           ! Mg-27
    calcVars%alj(66)= calcVars%alj(65)*(p-27.d0)*(exn-28.d0)/756.d0           ! Mg-28

    return
! ======================================================================
  end subroutine setALJ
