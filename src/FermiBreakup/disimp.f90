
  subroutine disimp (fbuObj, k, amk, pn0, tn)

! ======================================================================
!
!   Called by: RASTAR
!
!   Calls: RANG
!
!    Modified to incorporate ISOTR subroutine inline.
!
!    "Last" change: 13-Aug-2003 BY NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, July 2018 (creation of FermiBreakup class).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fermiBreakupParams, only : zro, one, two, twpi

    implicit none
    class(FermiBreakup), intent(inout) :: fbuObj
    integer(int32),      intent(in   ) :: k
    real(real64),        intent(in   ) :: amk(fbuObj%options%recNumNucleons)
    real(real64),        intent(  out) :: pn0(3, fbuObj%options%recNumNucleons)
    real(real64),        intent(in   ) :: tn

    integer(int32) :: i, kk, kl, l, lk
    real(real64)   :: amp, amr, csi, ct, fi, fk, pex, pmc, rnksi, smk, &
         & st, temp, temp2, tkm, tkn, tpr
    real(real64), dimension(3) :: anl = zro, pnc = zro, vrs = zro

! ======================================================================

    smk = zro
    do i = 1,k
       smk = smk + amk(i)
    end do
    do i = 1,3
       vrs(i) = zro
    end do
    tkn = tn
    kl = k - 1
    do l = 1,kl
       lk = k - l + 1
       amp = amk(lk)
       smk = smk - amp
       amr = smk
       tpr = tkn
       if (lk >= 3) then
10        csi = fbuObj%rang()
          kk = lk - 1
          if (kk >= 3) then
             pex = one/(1.5d0*dble(kk) - two)
             csi = csi**pex
          endif
          temp = csi*(one - csi)
          if (temp < 0.0d0) then
             temp = 0.01d0
             write(fbuObj%io%message,1100) "63"
             call fbuObj%io%print(4, 3, fbuObj%io%message)
          end if
          fk = two*sqrt(temp)
          if (fbuObj%rang() > fk) go to 10
          rnksi = csi
          tkm = tkn*rnksi
          tpr = tkn - tkm
       endif
       temp = amp + amr
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(fbuObj%io%message,1000) "74"
          call fbuObj%io%print(4, 3, fbuObj%io%message)
       end if
       temp2 = two*((amp*amr)/(temp))*tpr
       if (temp2 < 0.0d0) then
          temp2 = 0.01d0
          write(fbuObj%io%message,1100) "79"
          call fbuObj%io%print(4, 3, fbuObj%io%message)
       end if
       pmc = sqrt(temp2)
       ct = one - two*fbuObj%rang()
       st = sqrt(abs(one - ct**2))
       anl(3) = ct
       fi = twpi*fbuObj%rang()
       anl(1) = cos(fi)*st
       anl(2) = Sin(fi)*st
       do i = 1,3
          pnc(i) = pmc*anl(i)
          pn0(i,lk) = vrs(i)*amp + pnc(i)
          temp = amr
          if (temp < div0Lim .and. temp > -div0Lim) then
             temp = div0Lim
             write(fbuObj%io%message,1000) "95"
             call fbuObj%io%print(4, 3, fbuObj%io%message)
          end if
          vrs(i) = vrs(i) - pnc(i)/temp
       end do
       tkn = tkm
    end do
    do i = 1,3
       pn0(i,1) = vrs(i)*amr
    end do

    return
! ======================================================================
1000 format ("Divide by zero error prevented in 'disimp.f90', line(s) ", A)
1100 format ("Square root error prevented in 'disimp.f90', line(s) ", A)
! ======================================================================
  end subroutine disimp
