
  subroutine prrdis (gsmObj, anucl, znucl, mb0, fn)

! ======================================================================
!
!   This subroutine prints out the information on residual nuclei after
!   the cascade, preequilibrium, evaporation and fission stages of
!   the reaction.
!
!   Written by K. K. Gudima, December, 2004.
!   Edited by A. J. Sierk, LANL T-16, January, 2005.
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro, one, fiv

    implicit none
    class(GSM),     intent(inout) :: gsmObj
    real(real64),   intent(in   ) :: anucl
    real(real64),   intent(in   ) :: znucl
    integer(int32), intent(in   ) :: mb0
    real(real64),   intent(in   ) :: fn

    integer(int32) :: i, iamax, ip, ipz, izmax, j
    real(real64)   :: abeg, am1, am2, av, dam, den1a, ex1, ex2, fct0a, &
         & pm1, pm2, s, temp, zbeg

    integer(int32), dimension(250) :: iar=0_int32, izr=0_int32
    real(real64),   dimension(  5) :: den2=zro, pro=zro

! ======================================================================

    real(real64) :: rdis, dex, dpm
    common /resdis/  rdis(5,5,250), dex, dpm

! ======================================================================

    abeg = anucl + fiv
    zbeg = znucl + fiv
    izmax = nint(zbeg)
    izmax = min(izmax, 247)
    iamax = nint(abeg) + mb0
    iamax = min(iamax, 247)
    ip = nint(abeg)
    ipz = nint(zbeg)
    temp = fn
    if ( temp < div0Lim .and. temp > -div0Lim ) then
       temp = div0Lim
       write(gsmObj%io%message, 5000) "50"
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    fct0a = one/temp
    do i=1,ip
       iar(i) = ip - i + 1
    end do
    do i=1,ipz
       izr(i) = ipz - i + 1
    end do
!   Print out the A-distribution of residual nuclei
    write (31, 1000)
    do i=1,ip
       s = zro
       if (i <= 247) then
          do j=1,5
             pro(j) = rdis(j,1,i)*fct0a
             s = s + pro(j)
          end do
          if (s > 1.d-10) write (31, 1100) iar(i), pro
       endif
    end do
    do j=1,5
       den2(j) = zro
       if (rdis(j,1,250) > zro) then
          den2(j) = one/rdis(j,1,250)
          pro(j) = rdis(j,1,248)*den2(j)
       else
          pro(j) = zro
       endif
    end do
    write (31, 1200) pro
    do j=1,5
       if (den2(j) > zro) then
          av = pro(j)
          pro(j) = sqrt(abs(rdis(j,1,249)*den2(j) - av**2))
       else
          pro(j) = zro
       endif
    end do
    write (31, 1300) pro
    do j=1,5
       pro(j) = rdis(j,1,250)*fct0a
    end do
    write (31, 1400) pro
!
!   Print out the Z-distribution of residual nuclei:
    write (31, 1500)
    do i=1,ipz
       s = zro
       if (i <= 247) then
          do j=1,5
             pro(j) = rdis(j,2,i)*fct0a
             s = s + pro(j)
          end do
          if (s > 1.d-10) write (31, 1600) izr(i), pro
       endif
    end do
!  Average value of Z:
    do j=1,5
       den2(j) = zro
       if (rdis(j,2,250) > zro) then
          den2(j) = one/rdis(j,2,250)
          pro(j) = rdis(j,2,248)*den2(j)
       else
          pro(j) = zro
       endif
    end do
    write (31, 1700) pro
!  Width of Z distribution:
    do j=1,5
       if (den2(j) > zro) then
          av  = pro(j)
          pro(j) = sqrt(abs(rdis(j,2,249)*den2(j) - av**2))
       else
          pro(j) = zro
       endif
    end do
    write (31, 1800) pro
!  Normalization of distribution; should be 1.0000:
    do j=1,5
       pro(j) = rdis(j,2,250)*fct0a
    end do
    write (31, 1900) pro
!
!   Print out the E*-distribution of residual nuclei:
    write (31, 2000)
    do i=1,247
       ex1 = dble(i-1)*dex
       ex2 = ex1 + dex
       s = zro
       temp = fn*dex
       if ( temp < div0Lim .and. temp > -div0Lim ) then
          temp = div0Lim
          write(gsmObj%io%message, 5000) "145"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       den1a = one/(temp)
       do j=1,5
          pro(j) = rdis(j,3,i)*den1a
          s = s + pro(j)
       end do
       if (s > 1.d-10) write (31, 2100) ex1, ex2, pro
    end do
    do j=1,5
       den2(j) = zro
       if (rdis(j,3,250) > zro) then
          den2(j) = one/rdis(j,3,250)
          pro(j) = rdis(j,3,248)*den2(j)
       else
          pro(j) = zro
       endif
    end do
    write (31, 2200) pro
    do j=1,5
       if (den2(j) > zro) then
          av = pro(j)
          pro(j) = sqrt(abs(rdis(j,3,249)*den2(j) - av**2))
       else
          pro(j) = zro
       endif
    end do
    write (31, 2300) pro
    do j=1,5
       pro(j) = rdis(j,3,250)*fct0a
    end do
    write (31, 2400) pro
!
!   Print out the P-distribution of residual nuclei:
    write (31, 2500)
    do i=1,247
       pm1 = dble(i-1)*dpm
       pm2 = pm1 + dpm
       s = zro
       temp = fn*dpm
       if ( temp < div0Lim .and. temp > -div0Lim ) then
          temp = div0Lim
          write(gsmObj%io%message, 5000) "184"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       den1a = one/(temp)
       do j=1,5
          pro(j) = rdis(j,4,i)*den1a
          s = s + pro(j)
       end do
       if (s > 1.d-10) write (31, 2600) pm1, pm2, pro
    end do
    do j=1,5
       if (rdis(j,4,250) > zro) then
          den2(j) = one/rdis(j,4,250)
          pro(j) = rdis(j,4,248)*den2(j)
       else
          pro(j) = zro
       endif
    end do
    write (31, 2700) pro
    do j=1,5
       if (den2(j) > zro) then
          av  = pro(j)
          pro(j) = sqrt(abs(rdis(j,4,249)*den2(j) - av**2))
       else
          pro(j) = zro
       endif
    end do
    write (31, 2800) pro
    do j=1,5
       pro(j) = rdis(j,4,250)*fct0a
    end do
    write (31, 2900) pro
!
!   Print out the L-distribution of residual nuclei:
    write (31, 3000)
    dam = one
    do i=1,247
       am1 = dble(i-1)*dam
       am2 = am1 + dam
       s = zro
       temp = fn*dam
       if ( temp < div0Lim .and. temp > -div0Lim ) then
          temp = div0Lim
          write(gsmObj%io%message, 5000) "226"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       den1a = one/(temp)
       do j=1,5
          pro(j) = rdis(j,5,i)*den1a
          s = s + pro(j)
       end do
       if (s > 1.d-10) write (31, 3100) am1, am2, pro
    end do
    do j=1,5
       pro(j) = zro
       if (rdis(j,5,250) > zro) then
          den2(j) = one/rdis(j,5,250)
          pro(j) = rdis(j,5,248)*den2(j)
       endif
    end do
    write (31, 3200) pro
    do j=1,5
       if (den2(j) > zro) then
          av  = pro(j)
          pro(j) = sqrt(abs(rdis(j,5,249)*den2(j) - av**2))
       else
          pro(j) = zro
       endif
    end do
    write (31, 3300) pro
    do j=1,5
       pro(j) = rdis(j,5,250)*fct0a
    end do
    write (31, 3400) pro
    return

! ======================================================================
1000 format (/18x,' Mass distributions of nuclei:'/38x, &
          & 'at start of evap, which:  just prior to' &
          & /12x,'after cascade',' after preeq   evap. only', &
          & 4x,'fission',6x,'fission')
1100 format (3x,'A =',i4,5(4x,1pe9.3))
1200 format (4x,'<A> =',1x,5(4x,1pe9.3))
1300 format ('St Dev A =',1x,5(4x,1pe9.3))
1400 format (3x,'norm =',1x,5(4x,1pe9.3))
1500 format (/18x,' Charge distributions of nuclei:'/38x, &
          & 'at start of evap, which:  just prior to' &
          & /12x,'after cascade',' after preeq   evap. only', &
          & 4x,'fission',6x,'fission')
1600 format (3x,'Z =',i4,5(4x,1pe9.3))
1700 format (4x,'<Z> =',1x,5(4x,1pe9.3))
1800 format ('St Dev Z =',1x,5(4x,1pe9.3))
1900 format (3x,'norm =',1x,5(4x,1pe9.3))
2000 format (/10x,' Excitation energy distributions [1/MeV]', &
          & ' of nuclei:'/40x, &
          & 'at start of evap, which:  just prior to' &
          & /4x,'E*(MeV)',3x,'after cascade', &
          & ' after preeq   evap. only',4x, &
          & 'fission',6x, 'fission')
2100 format (1x,f5.0,'-',f5.0,5(4x,1pe9.3))
2200 format (5x,'<E*> =',1x,5(4x,1pe9.3))
2300 format ('St Dev E* =',1x,5(4x,1pe9.3))
2400 format (5x,'norm =',1x,5(4x,1pe9.3))
2500 format (/10x,' Linear momentum distributions [1/MeV/c]', &
          & ' of nuclei:'/40x, &
          & 'at start of evap, which:  just prior to' &
          & /4x,'P(MeV/c)',2x,'after cascade', &
          & ' after preeq   evap. only',4x, &
          & 'fission',6x, 'fission')
2600 format (1x,f5.0,'-',f5.0,5(4x,1pe9.3))
2700 format (5x,' <P> =',1x,5(4x,1pe9.3))
2800 format (' St Dev P =',1x,5(4x,1pe9.3))
2900 format (5x,'norm =',1x,5(4x,1pe9.3))
3000 format (/14x,' Angular momentum distributions [1/hbar]', &
          & ' of nuclei:'/39x, &
          & 'at start of evap, which:  just prior to' &
          & /4x,'  L',5x,'after cascade', &
          & ' after preeq   evap. only',4x, &
          & 'fission',6x, 'fission')
3100 format (1x,f4.0,'-',f4.0,5(4x,1pe9.3))
3200 format (3x,' <L> =',1x,5(4x,1pe9.3))
3300 format ('St Dev L =',1x,5(4x,1pe9.3))
3400 format (3x,'norm =',1x,5(4x,1pe9.3))
5000 format("Divide by zero error prevented in 'prrdis.f90', line(s) ", A)
! ======================================================================
  end subroutine prrdis
