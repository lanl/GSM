
  subroutine crack (fbuObj, k, ia, iz, np, iaz)

! ======================================================================
!
!    Called by: RAZVAL
!
!    Calls: DIVAZ
!
!    Last change: 13-Aug-2003 BY NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32

    implicit none
    class(FermiBreakup), intent(inout) :: fbuObj
    integer(int32),      intent(in   ) :: k
    integer(int32),      intent(in   ) :: ia
    integer(int32),      intent(in   ) :: iz
    integer(int32),      intent(inout) :: np
    real(real64),        intent(  out), dimension(mp) :: iaz

    integer(int32) :: i, ka, km, kp, kz, la, lz, m, mi, mk, n, n1, n2, &
         & nc, ni, nk, nm, nn, nq, nr, nw
    integer(int32), dimension(fbuObj%options%recNumNucleons) :: maz

    type(fbuCrackCalculation) :: crackCalc

! ======================================================================

    kp = k
    nw = 0
    call fbuObj%divaz (ia, iz, ia, iz, kp, nr, crackCalc)
    if (nr < 1) return
    nk = np - kp

    do n = 1,nr
       nk = nk + kp
       iaz(nk+1) = 100*crackCalc%mpa(n) + crackCalc%mpz(n)
       iaz(nk+2) = 100*crackCalc%mra(n) + crackCalc%mrz(n)
       nw = nw + 1
    end do
    nn = kp*nw + np
    if (kp < 3) then
       np = nn
    else
       nq = nw
       nw = 0
       nc = 0
       mk = kp - 1
       do m = 2,mk
          km = kp - m + 1
          nk = np - kp
          mi = m - 1
10        continue
          nk = nk + kp
          nm = nk + m
20        continue
          do i = 1,mi
             maz(i) = iaz(nk+i)
          end do
          la = iaz(nm-1)/100
          lz = iaz(nm-1) - la*100
          ka = iaz(nm)/100
          kz = iaz(nm) - ka*100
          call fbuObj%divaz (la, lz, ka, kz, km, nr, crackCalc)
          if (nr < 1) go to 30
          iaz(nm) = 100*crackCalc%mpa(1) + crackCalc%mpz(1)
          iaz(nm+1) = 100*crackCalc%mra(1) + crackCalc%mrz(1)
          nw = nw + 1
          if (nr >= 2) then
             do n = 2,nr
                do i = 1,mi
                   iaz(nn+i) = maz(i)
                end do
                nm = nn + m
                iaz(nm) = 100*crackCalc%mpa(n) + crackCalc%mpz(n)
                iaz(nm+1) = 100*crackCalc%mra(n) + crackCalc%mrz(n)
                nn = nn + kp
                nw = nw + 1
             end do
          endif
          nc = nc + 1
          if (nc < nq) go to 10
          go to 40
30        nq = nq - 1
          nn = nn - kp
          if (nq < 1) then
             np = nn
             return
          endif
          n1 = nk + 1
          n2 = nn
          do ni = n1,n2
             iaz(ni) = iaz(ni+kp)
          end do
          if (nc < nq) go to 20
40        nq = nw
          nc = 0
          nw = 0
       end do
       np = nn
    endif

    return
! ======================================================================
  end subroutine crack
