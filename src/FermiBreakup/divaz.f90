
  subroutine divaz (fbuObj, la, lz, ia, iz, k, nr, crackCalc)

! ======================================================================
!
!    Called by: RAZVAL
!
!     Last change: 13-Aug-2003 BY NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fermiBreakUpData, only : ms, mpraz, min_fermi_AData

    implicit none
    class(FermiBreakup), intent(inout) :: fbuObj
    integer(int32),      intent(in   ) :: la
    integer(int32),      intent(in   ) :: lz
    integer(int32),      intent(in   ) :: ia
    integer(int32),      intent(in   ) :: iz
    integer(int32),      intent(in   ) :: k
    integer(int32),      intent(  out) :: nr
    type(fbuCrackCalculation), intent(  out ) :: crackCalc

    integer(int32) :: ip, j, j1, j2, jp, jz, m, map, mar, mzp, mzr, &
         & n, nm, nmk, nmn, nmx, npr
    real(real64)   :: temp

! ======================================================================

    npr = mpraz
    m = 0
    nmx = la
    nmk = ia - k + 1
    nmx = min (nmx, nmk)
    temp = k
    if (temp < div0Lim .and. temp > -div0Lim) then
       temp = div0Lim
       write(fbuObj%io%message,1000) "732"
       call fbuObj%io%print(4, 3, fbuObj%io%message)
    end if
    nmn = (ia - 1)/temp + 1
    nmn = min(nmx, nmn)
    nm = nmx - nmn + 1
    do n = 1,nm
       map = nmx - n + 1
       mar = ia - map
        j1 = iz - mar
        j1 = max (j1,0)
        jz = 0
        if (map == la) jz = lz
        j1 = max(j1, jz)
        j2 = iz
        if (k < 3 .and. map == mar) j2 = iz/2
        j2 = min (map, j2)
        j1 = j1 + 1
        j2 = j2 + 1
        if (j2 >= j1) then
           do j = j1,j2
              mzp = j - 1
              mzr = iz - mzp
              jp = mzp + 1
              ip = map - mzp + 1
!  KKG 11/02/04:
              if (ip <= 0 .or. jp <= 0 .or. ip > min_fermi_AData .or. jp > (min_fermi_AData-1)) go to 10
              if (ms(ip,jp) >= 1) then
                 if (k <= 2) then
                    jp = mzr + 1
                    ip = mar - mzr + 1
                    !  KKG 11/02/04:
                    if (ip <= 0 .or. jp <= 0 .or. ip > 12 .or. jp > 11)    &
                         go to 10
                 endif
                 m = m + 1
                 if (m > mpraz) then
                    write(fbuObj%io%message, 100) npr, n
                    call fbuObj%io%print(3, 3, fbuObj%io%message)
                    go to 10
                 endif
                 crackCalc%mpa(m) = map
                 crackCalc%mpz(m) = mzp
                 crackCalc%mra(m) = mar
                 crackCalc%mrz(m) = mzr
              endif
10            continue
           end do
        endif
     end do
     nr = m

     return
! ======================================================================
100  format ("The internal Fermi Break-up arrays will be ", &
          & "exceeded in 'divaz' routine. The fragment will be ignored.", &
          & /, 9x, i3, " particles produced, n=", i3)
1000 format("Divide by zero error prevented in 'divaz.f90', line ", A)
! ======================================================================
   end subroutine divaz
