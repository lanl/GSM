
   subroutine razval (fbuObj, an, zn, up, nf, np, calc)
 
! ======================================================================
!
!    Called by: RASTAR
!
!    Calls: CRACK, TCUL, WECHAN, WEPT
!
!    Last change: 13-Aug-2003 BY NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Modified by R E Prael, 03 May 2006.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class)
!
! ======================================================================

     use, intrinsic :: iso_fortran_env, only: int32, real64
     use fermiBreakupParams, only : zro
     use fermiBreakUpData, only : mp, mf

     implicit none
     class(FermiBreakup),  intent(inout) :: fbuObj
     real(real64),         intent(in   ) :: an
     real(real64),         intent(in   ) :: zn
     real(real64),         intent(in   ) :: up
     integer(int32),       intent(  out) :: nf
     integer(int32),       intent(  out) :: np
     type(fbuInterimCalc), intent(  out) :: calc

     integer(int32) :: i, ia, iz, jz, k, ka, kp, mpt, n, ni, nk, nt
     integer(int32), dimension(fbuObj%options%recNumNucleons) :: ntv
     real(real64) :: daz, rmq, srq, temp, tq
 
! ====================================================================== 

     ia = nint(an)
     mpt = mp - ia
     iz = nint(zn)
     daz = up + fbuObj%wept (ia, iz)
     nf = 0
     np = 0
     srq = zro
     ka = ia - 1
     if (ka > 2) then
        do k = 2,ka
           kp = k
           nt = np
           nk = nt - kp
           call fbuObj%crack (kp, ia, iz, np, calc%iaz)
           if (np > nt) then
10            continue
              nk = nk + kp
              if (nk > mpt) then
                 ! iaz(X) array will be exceeded; warn client and approximate
                 write (fbuObj%io%message, 100) "iaz", ia, iz, kp, nf, nt, nk, np
                 call fbuObj%io%print(3, 3, fbuObj%io%message)
                 nk = mpt
              endif
20            do i = 1,k
                 ntv(i) = calc%iaz(nk+i)
              end do
              tq = fbuObj%tcul (kp, ntv)
              rmq = fbuObj%wechan (ia, kp, ntv, daz, tq)
!    REP 03 May 2006:
              if (rmq > zro) then
                 nf = nf + 1
                 if (nf > mf) then
                    ! wnq(X) array will be exceeded; warn client and approximate
                    write (fbuObj%io%message, 100) "wnq", ia, iz, kp, nf, nt, nk, np
                    call fbuObj%io%print(3, 3, fbuObj%io%message)
                    nf = mf
                 endif
                 srq = srq + rmq
                 calc%wnq(nf) = rmq
                 if ((np - nk) > kp) go to 10
              else
                 np = np - kp
                 ni = nk + 1
                 if (np >= ni) then
                    do n = ni,np
                       calc%iaz(n) = calc%iaz(n+kp)
                    end do
                    go to 20
                 else
!   REP 03 May 2006:
                    if ((np <= nt) .and. (kp > 2) .and. srq > zro) then
                       if (nf >= 1) then
                          do n = 1,nf
                             temp = srq
                             if (temp < div0Lim .and. temp > -div0Lim) then
                                temp = div0Lim
                                write(fbuObj%io%message,1000) "882"
                                call fbuObj%io%print(4, 3, fbuObj%io%message)
                             end if
                             calc%wnq(n) = calc%wnq(n)/temp
                          end do
                       endif
                       return
                    endif
                 endif
              endif
           endif
        end do
     endif
     do i = 1,ia
        jz = 1
        if (i > iz) jz = 0
        ntv(i) = 100 + jz
     end do
     tq = fbuObj%tcul (ia, ntv)
     rmq = fbuObj%wechan (ia, ia, ntv, daz, tq)
     if (rmq >= zro) then
        nf = nf + 1
        calc%wnq(nf) = rmq
        srq = srq + rmq
        do i = 1,ia
           calc%iaz(np+i) = ntv(i)
        end do
        np = np + ia
     endif
     if (nf >= 1) then
        temp = srq
        if (temp < div0Lim .and. temp > -div0Lim) then
           temp = div0Lim
           write(fbuObj%io%message,1000) "921"
           call fbuObj%io%print(4, 3, fbuObj%io%message)
        end if
        do n = 1,nf
           calc%wnq(n) = calc%wnq(n)/temp
        end do
     endif

     return
! ======================================================================
100  format ("The Fermi Break-up array '", A, "' will be ", &
          & "exceeded. ", /, 9x, "Approximating with the largest available value ", &
          & "instead.", /, 9x, "A=", i3, ", Z=", i3, ", kp=", &
          & i5, ", nf=", i5, ", nt=", i5, ", nk=", i5, ", np=", i5)
1000 format("Divide by zero error prevented in 'razval.f90', line ", A)
! ======================================================================
   end subroutine razval
