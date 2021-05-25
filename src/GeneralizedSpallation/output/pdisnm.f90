
  subroutine pdisnm (gsmObj, fn)

! ======================================================================
!
!   This subroutine prints out the multiplicity distribution of neutrons
!   produced:
!              1) after all stages,
!              2) in cascade stage,
!              3) in preeq. stage,
!              4) in evaporation without fission,
!              5) in prefission stage,
!              6) in postfission stage.
!
!   Written by K. K. Gudima, December, 2004.
!   Modified by A. J. Sierk, LANL T-16, January, 2005.
!   Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro, one

    implicit none
    class(GSM),   intent(inout) :: gsmObj
    real(real64), intent(in   ) :: fn

    integer(int32) :: i, j, nn
    real(real64)   :: av, fct0a, s
    real(real64), dimension(6) :: fct1, pro

! ======================================================================

    real(real64) :: disnm
    common /disnmu/ disnm(6,155)

! ======================================================================

    write (31, 1000)
    fct0a = one/fn
    do i=1,152
       nn = i - 1
       s = zro
       do j=1,6
          pro(j) = disnm(j,i)*fct0a
          s = s + pro(j)
       end do
       if (s > zro) write (31, 1100) nn, pro
    end do
!   End of i loop ^
    do j=1,6
       pro(j) = zro
       fct1(j) = zro
       if (disnm(j,155) > zro) then
          fct1(j) = one/disnm(j,155)
          pro(j) = disnm(j,153)*fct1(j)
       endif
    end do
    write (31, 1200) pro
    do j=1,6
       if (pro(j) > zro) then
          av = pro(j)
          pro(j) = sqrt(abs(disnm(j,154)*fct1(j) - av**2))
       endif
    end do
    write (31, 1300) pro
    do j=1,6
       pro(j) = disnm(j,155)*fct0a
    end do
    write (31, 1400) pro
    return

! ======================================================================
1000 format (/22x,'Neutron-multiplicity probability:'// &
         & 5x,'Nn',6x,'Total',6x,'Cascade',4x,'Preequil.',3x, &
         & 'Evap. res.',2x,'Pre-fiss.',3x,'Post-fiss.')
1100 format (4x,i3,2x,6(2x,1pe9.3,1x))
1200 format (4x,'<n> =',6(2x,1pe9.3,1x))
1300 format ('St dv n =',6(2x,1pe9.3,1x))
1400 format (3x,'norm =',6(2x,1pe9.3,1x)/)
! ======================================================================
  end subroutine pdisnm
