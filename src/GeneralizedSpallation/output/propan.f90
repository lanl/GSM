
  subroutine propan (gsmObj, fnfis)

! ======================================================================
!
!   This subroutine prints out the distribution of fission-fragment
!   opening angles (in the lab. syst.) in different bins of neutron
!   multiplicity.
!
!   Written by K. K. Gudima, December, 2004.
!   Edited by A. J. Sierk, LANL T-16, January, 2005.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro, one

    implicit none
    class(GSM),   intent(inout) :: gsmObj
    real(real64), intent(in   ) :: fnfis

    integer(int32) :: i, j
    real(real64)   :: av, fct0, s, temp, tet1, tet2
    real(real64), dimension(7) :: fct1=zro, pro=zro

! ======================================================================

    real(real64) :: opan, dth12
    common /fisopa/  opan(7,185), dth12

! ======================================================================

    temp = fnfis
    if ( temp < div0Lim .and. temp > -div0Lim ) then
       temp = div0Lim
       write(gsmObj%io%message, 3000) "40, 68"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if

    write (31, 1000)
    fct0 = one/(temp*dth12)
    do i=1,180
       tet1 = dble(i-1)*dth12
       tet2 = tet1 + dth12
       s = zro
       do j=1,7
          pro(j) = opan(j,i)*fct0
          s = s + pro(j)
       end do
       if (s > zro) write (31, 1100) tet1, tet2, pro
    end do
    do j=1,7
       pro(j) = zro
       fct1(j) = zro
       if (opan(j,185) > zro) then
          fct1(j) = one/opan(j,185)
          pro(j) = opan(j,183)*fct1(j)
       endif
    end do
    write (31, 1200) pro
    do j=1,7
       if (pro(j) > zro) then
          av = pro(j)
          pro(j) = sqrt(abs(opan(j,184)*fct1(j) - av**2))
       endif
    end do
    write (31, 1300) pro
    do j=1,7
       pro(j) = opan(j,185)/temp
    end do
    write (31, 1400) pro
    return

! ======================================================================
1000 format (/3x,'Distribution of fission-fragment opening angles ', &
          & '[1/deg.] (lab.sys.) in different'/29x,' bins ', &
          & 'of neutron multiplicity:'/1x,'theta(deg.)',1x, &
          & 'All events ',' n = 0-5  ',' n = 6-8 ',' n = 9-12 ', &
          & ' n = 13-15 ','n = 16-19 ','  n > 20')
1100 format (1x,f4.0,' - ',f4.0,1x,7(1x,1pe9.3))
1200 format (4x,'<thet> = ',7(1x,1pe9.3))
1300 format ('St dv thet = ',7(1x,1pe9.3))
1400 format (5x,'norm. = ',7(1x,1pe9.3)/)
3000 format("Divide by zero error prevented in 'propan.f90', line(s) ", A)
! ======================================================================
  end subroutine propan
