
  subroutine opandist (gsmObj, event, fisevent, nn)

! ======================================================================
!
!      Accumulate distribution of fission fragments' opening angles.
!      nn is total number of produced neutrons.
!      bf12 contains the fragment velocity vectors.
!
!   Written by K. K. Gudima, December, 2004.
!   Edited by A. J. Sierk, LANL T-16, January, 2005.
!   Converted to separate SR by A. J. Sierk, LANL T-16, March, 2006.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: one, radiantodegree

    implicit none
    class(GSM),       intent(inout) :: gsmObj
    class(EventData), intent(in   ) :: event
    logical,          intent(in   ) :: fisevent
    integer(int32),   intent(in   ) :: nn

    integer(int32) :: it12, jk
    real(real64)   :: b1, b2, ct12, st12, t12, temp

! ======================================================================

    real(real64) :: opan, dth12
    common /fisopa/ opan(7,185), dth12

! ======================================================================

    if (.not.fisevent) return
    b1 = sqrt(event%fissFrag(1)%linearMomFrac(1)**2 + &
         &    event%fissFrag(1)%linearMomFrac(2)**2 + &
         &    event%fissFrag(1)%linearMomFrac(3)**2)
    b2 = sqrt(event%fissFrag(2)%linearMomFrac(1)**2 + &
         &    event%fissFrag(2)%linearMomFrac(2)**2 + &
         &    event%fissFrag(2)%linearMomFrac(3)**2)
    temp = b1*b2
    if ( temp < div0Lim .and.temp > -div0Lim ) then
       temp = div0Lim
       write(gsmObj%io%message,1000) "44"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    ! Angle between the two fission fragments?
    ct12 = (event%fissFrag(1)%linearMomFrac(1) * &
         &    event%fissFrag(2)%linearMomFrac(1) + &
         & event%fissFrag(1)%linearMomFrac(2) * &
         &    event%fissFrag(2)%linearMomFrac(2) + &
         & event%fissFrag(1)%linearMomFrac(3) * &
         &    event%fissFrag(2)%linearMomFrac(3)) / (temp)
    st12 = sqrt(abs(one - ct12**2))
    t12 = atan2(st12,ct12)*radiantodegree
    if ( dth12 < div0Lim .and.dth12 > -div0Lim ) then
       dth12 = div0Lim
       write(gsmObj%io%message,1000) "52"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    it12 = int(t12/dth12) + 1
    opan(1,it12) = opan(1,it12) + one
    opan(1,183)  = opan(1,183)  + t12
    opan(1,184)  = opan(1,184)  + t12*t12
    opan(1,185)  = opan(1,185)  + one
    if (nn <= 5) then
!   n multiplictiy 0-5
       jk = 2
    elseif (nn <= 8) then
!   n multiplictiy 6-8
       jk = 3
    elseif (nn <= 12) then
!   n multiplictiy 9-12
       jk = 4
    elseif (nn <= 15) then
!   n multiplictiy 13-15
       jk = 5
    elseif (nn <= 19) then
!   n multiplictiy 16-19
       jk = 6
    elseif (nn >= 20) then
!   n multiplicity 20-.....
       jk = 7
    endif
    opan(jk,it12) = opan(jk,it12) + one
    opan(jk,183)  = opan(jk,183)  + t12
    opan(jk,184)  = opan(jk,184)  + t12*t12
    opan(jk,185)  = opan(jk,185)  + one
    return

! ======================================================================
1000 format("Divide by zero error prevented in 'opandist.f90', line(s) ", A)
! ======================================================================
  end subroutine opandist
