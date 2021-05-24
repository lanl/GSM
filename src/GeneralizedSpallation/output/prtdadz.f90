
  subroutine prtdadz (gsmObj, anucl, znucl, mb0, sigin, fn)

! ======================================================================
!
!   This subroutine prints out average kinetic energies of
!   fragments and fragment yields emitted in forward/backward direction.
!
!   Written by K. K. Gudima, November, 2004.
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
    real(real64),   intent(in   ) :: sigin
    real(real64),   intent(in   ) :: fn

    integer(int32) :: i, ia, iaa, iamax, ip, ipt, ipz, iz, iz1, iz2, &
         & iz3, izmax, k, lpr, lprm, npt
    real(real64)   :: abeg, atet, atk, atks, datet, datk, datks, dfbs, &
         & dtets, dvzs, dy1, dy2, dy3, dys1, dys2, dys3, dyy, dyyb, &
         & dyybs, dyyf, dyyfs, dyys, faca, fbs, ppt, ss, tets, vzs, &
         & y1, y2, y3, ys1, ys2, ys3, yy, yyb, yybs, yyf, yyfs, yys, &
         & zbeg, tempfn

    integer(int32), dimension(351) :: iar=0_int32, izr=0_int32
    real(real64),   dimension(351) :: ar=zro, avz=zro, &
         & davz=zro, fb=zro, dfb=zro
    real(real64),   dimension(151) :: dasum=zro
    real(real64),   dimension( nint(znucl+12) ) :: z

    logical :: printit, prtz1, prtz2, prtz3
    logical, dimension(151) :: prtz = .FALSE.

! ======================================================================

    real(real64) :: spec, ang, chan, dadz
    common /result/  spec(9,200), ang(9,20), chan(218), dadz(351,151)
    real(real64) :: dadz4
    common /result1/ dadz4(4,351,151)

! ======================================================================

    abeg = anucl + fiv
    zbeg = znucl + fiv
    izmax = nint(zbeg)
    izmax = min(izmax, 148)
    iamax = nint(abeg) + mb0
    iamax = min(iamax, 348)
    lprm = (izmax + 3)/3
    lprm = max(lprm, 1)
    tempfn = fn
    if ( tempfn < div0Lim .and. tempfn > -div0Lim ) then
       tempfn = div0Lim
       write(gsmObj%io%message,5000) tempfn
       call gsmObj%io%print(4, 3, gsmObj%io%message)
    end if
    faca = sigin/fn
    ip = nint(abeg)
    ipz = nint(zbeg)
    do i = 1,ip
       iar(i) = ip - i + 1
    end do
    do i = 1,ipz
       izr(i) = ipz - i + 1
    end do
!
!   Print out forward isotope formation cross sections:
!
    do iz = 1,izmax + 1
       z(iz) = max (-one, zbeg - dble(iz-1))
       dasum(iz) = zro
       do ia = 1,iamax
          ar(ia) = max(abeg - dble(ia-1), zro)
          dasum(iz) = dasum(iz) + dadz4(3,ia,iz)
       end do
       prtz(iz) = z(iz) >= zro .and. dasum(iz) > zro
    end do
    write (31, 1000)
    do lpr = 1,lprm
       iz1 = 3*lpr - 2
       iz2 = 3*lpr - 1
       iz3 = 3*lpr
       prtz1 = prtz(iz3)
       prtz2 = .not.prtz(iz3) .and. prtz(iz2)
       prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. prtz(iz1)
       if (prtz1) then
          write (31, 1100) z(iz1), z(iz2), z(iz3)
       elseif (prtz2) then
          write (31, 1200) z(iz1), z(iz2)
       elseif (prtz3) then
          write (31, 1300) z(iz1)
       endif
       ppt = zro
       do k = 1,iamax
          y1 = dadz4(3,k,iz1)*faca
          dy1 = sqrt(abs(dadz4(3,k,iz1)))*faca
          y2 = dadz4(3,k,iz2)*faca
          dy2 = sqrt(abs(dadz4(3,k,iz2)))*faca
          y3 = dadz4(3,k,iz3)*faca
          dy3 = sqrt(abs(dadz4(3,k,iz3)))*faca
          printit = y1 > zro .or. y2 > zro .or. y3 > zro
          if (printit) then
!  KKG 11/15/04
             if (ar(k) > zro .and. ar(k) >= z(iz3)) then
                iaa = nint(ar(k))
                prtz1 = prtz(iz3)
                prtz2 = .not.prtz(iz3) .and. prtz(iz2)
                prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. &
                     & prtz(iz1)
                if (prtz1) then
                   ppt = ppt + one
                   write (31, 1400) iaa, y1, dy1, y2, dy2, y3, dy3
                elseif (prtz2) then
                   ppt = ppt + one
                   write (31, 1500) iaa, y1, dy1, y2, dy2
                elseif (prtz3) then
                   ppt = ppt + one
                   write (31, 1600) iaa, y1, dy1
                endif
             endif
          endif
       end do
       if (ppt > zro) then
          ipt = nint(ppt)
          ys1  = dasum(iz1)*faca
          dys1 = sqrt(abs(dasum(iz1)))*faca
          ys2  = dasum(iz2)*faca
          dys2 = sqrt(abs(dasum(iz2)))*faca
          ys3  = dasum(iz3)*faca
          dys3 = sqrt(abs(dasum(iz3)))*faca
          if (ys3 > zro) then
             write (31, 1700) ipt, ys1, dys1, ys2, dys2, ys3, dys3
          elseif (ys2 > zro) then
             write (31, 1800) ipt, ys1, dys1, ys2, dys2
          elseif (ys1 > zro) then
             write (31, 1900) ipt, ys1, dys1
          endif
       endif
    end do
    write (31, 2000)
    ss = faca
    write (31, 2100)
    atks = zro
    datks= zro
    yys = zro
    npt = 0
    do i = 1,ip
       yy = dadz4(3,i,150)*ss
       dyy = sqrt(abs(dadz4(3,i,150)))*ss
       yys = yys + dadz4(3,i,150)
       if (dadz4(3,i,150) > zro) then
          atk   = dadz4(3,i,151)/dadz4(3,i,150)
          atks  = atks + dadz4(3,i,151)
          datk  = sqrt(abs(dadz4(3,i,149)/dadz4(3,i,150) - atk**2))
          datks = datks + dadz4(3,i,149)
       else
          atk  = zro
          datk = zro
       endif
       printit = yy > zro
       if (printit) then
          npt = npt + 1
          write (31, 2200) iar(i), yy, dyy, atk, datk
       endif
    end do
    dyys = sqrt(abs(yys))*ss
    if (yys > zro) then
       atks  = atks/yys
       datks = sqrt(abs(datks/yys - atks**2))
    endif
    yys = yys*ss
    write (31, 2300) npt, yys, dyys, atks, datks
    write (31, 2400)
    yys = zro
    atks = zro
    datks= zro
    npt = 0
    do i = 1,ipz+1
       yy = dadz4(3,351,i)*ss
       dyy = sqrt(abs(dadz4(3,351,i)))*ss
       yys = yys + dadz4(3,351,i)
       if (dadz4(3,351,i) > zro) then
          atk   = dadz4(3,350,i)/dadz4(3,351,i)
          atks  = atks + dadz4(3,350,i)
          datk  = sqrt(abs(dadz4(3,349,i)/dadz4(3,350,i) - atk**2))
          datks = datks + dadz4(3,349,i)
       else
          atk  = zro
          datk = zro
       endif
       printit = yy > zro
       if (printit) then
          npt = npt + 1
          write (31, 2500) izr(i), yy, dyy, atk, datk
       endif
    end do
    dyys = sqrt(abs(yys))*ss
    if (yys > zro) then
       atks  = atks/yys
       datks = sqrt(abs(datks/yys - atks**2))
    endif
    yys = yys*ss
    write (31, 2600) npt, yys, dyys, atks, datks
!
!   Print out backward isotope formation cross sections:
!
    do iz = 1,izmax + 1
       z(iz) = max (-one, zbeg - dble(iz-1))
       dasum(iz) = zro
       do ia = 1,iamax
          ar(ia) = max(abeg - dble(ia-1), zro)
          dasum(iz) = dasum(iz) + dadz4(4,ia,iz)
       end do
       prtz(iz) = z(iz) >= zro .and. dasum(iz) > zro
    end do
    write (31, 2700)
    do lpr = 1,lprm
       iz1 = 3*lpr - 2
       iz2 = 3*lpr - 1
       iz3 = 3*lpr
       prtz1 = prtz(iz3)
       prtz2 = .not.prtz(iz3) .and. prtz(iz2)
       prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. prtz(iz1)
       if (prtz1) then
          write (31, 1100) z(iz1), z(iz2), z(iz3)
       elseif (prtz2) then
          write (31, 1200) z(iz1), z(iz2)
       elseif (prtz3) then
          write (31, 1300) z(iz1)
       endif
       ppt = zro
       do k = 1,iamax
          y1 = dadz4(4,k,iz1)*faca
          dy1 = sqrt(abs(dadz4(4,k,iz1)))*faca
          y2 = dadz4(4,k,iz2)*faca
          dy2 = sqrt(abs(dadz4(4,k,iz2)))*faca
          y3 = dadz4(4,k,iz3)*faca
          dy3 = sqrt(abs(dadz4(4,k,iz3)))*faca
          printit = y1 > zro .or. y2 > zro .or. y3 > zro
          if (printit) then
!  KKG 11/15/04
             if (ar(k) > zro .and. ar(k) >= z(iz3)) then
                iaa = nint(ar(k))
                prtz1 = prtz(iz3)
                prtz2 = .not.prtz(iz3) .and. prtz(iz2)
                prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. &
                     & prtz(iz1)
                if (prtz1) then
                   ppt = ppt + one
                   write (31, 1400) iaa, y1, dy1, y2, dy2, y3, dy3
                elseif (prtz2) then
                   ppt = ppt + one
                   write (31, 1500) iaa, y1, dy1, y2, dy2
                elseif (prtz3) then
                   ppt = ppt + one
                   write (31, 1600) iaa, y1, dy1
                endif
             endif
          endif
       end do
!   End of k loop ^
       if (ppt > zro) then
          ipt = nint(ppt)
          ys1  = dasum(iz1)*faca
          dys1 = sqrt(abs(dasum(iz1)))*faca
          ys2  = dasum(iz2)*faca
          dys2 = sqrt(abs(dasum(iz2)))*faca
          ys3  = dasum(iz3)*faca
          dys3 = sqrt(abs(dasum(iz3)))*faca
          if (ys3 > zro) then
             write (31, 1700) ipt, ys1, dys1, ys2, dys2, ys3, dys3
          elseif (ys2 > zro) then
             write (31, 1800) ipt, ys1, dys1, ys2, dys2
          elseif (ys1 > zro) then
             write (31, 1900) ipt, ys1, dys1
          endif
       endif
    end do
!   End of lpr loop ^
    write (31, 2800)
    ss = faca
    write (31, 2900)
    atks = zro
    datks= zro
    yys = zro
    npt = 0
    do i = 1,ip
       yy = dadz4(4,i,150)*ss
       dyy = sqrt(abs(dadz4(4,i,150)))*ss
       yys = yys + dadz4(4,i,150)
       if (dadz4(4,i,150) > zro) then
          atk   = dadz4(4,i,151)/dadz4(4,i,150)
          atks  = atks + dadz4(4,i,151)
          datk  = sqrt(abs(dadz4(4,i,149)/dadz4(4,i,150) - atk**2))
          datks = datks + dadz4(4,i,149)
       else
          atk = zro
          datk= zro
       endif
       printit = yy > zro
       if (printit) then
          npt = npt + 1
          write (31, 2200) iar(i), yy, dyy, atk, datk
       endif
    end do
!  End of i loop ^
    dyys = sqrt(abs(yys))*ss
    if (yys > zro) then
       atks = atks/yys
       datks= sqrt(abs(datks/yys - atks**2))
    endif
    yys = yys*ss
    write (31, 2300) npt, yys, dyys, atks, datks
    write (31, 3000)
    yys = zro
    atks = zro
    datks= zro
    npt = 0
    do i = 1,ipz+1
       yy = dadz4(4,351,i)*ss
       dyy = sqrt(abs(dadz4(4,351,i)))*ss
       yys = yys + dadz4(4,351,i)
       if (dadz4(4,351,i) > zro) then
          atk   = dadz4(4,350,i)/dadz4(4,351,i)
          atks  = atks + dadz4(4,350,i)
          datk  = sqrt(abs(dadz4(4,349,i)/dadz4(4,351,i) - atk**2))
          datks = datks + dadz4(4,349,i)
       else
          atk  = zro
          datk = zro
       endif
       printit = yy > zro
       if (printit) then
          npt = npt + 1
          write (31, 2500) izr(i), yy, dyy, atk, datk
       endif
    end do
!  End of i loop ^
    dyys = sqrt(abs(yys))*ss
    if (yys > zro) then
       atks = atks/yys
       datks= sqrt(abs(datks/yys - atks**2))
    endif
    yys = yys*ss
    write (31, 2600) npt, yys, dyys, atks, datks
!
!   Print out isotope average kinetic energies:
!
    do iz = 1,izmax + 1
       z(iz) = max (-one, zbeg - dble(iz-1))
       dasum(iz) = zro
       do ia = 1,iamax
          ar(ia) = max(abeg - dble(ia-1), zro)
          dasum(iz) = dasum(iz) + dadz(ia,iz)
       end do
       prtz(iz) = z(iz) >= zro .and. dasum(iz) > zro
    end do
    write (31, 3100)
    do lpr = 1,lprm
       iz1 = 3*lpr - 2
       iz2 = 3*lpr - 1
       iz3 = 3*lpr
       prtz1 = prtz(iz3)
       prtz2 = .not.prtz(iz3) .and. prtz(iz2)
       prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. prtz(iz1)
       if (prtz1) then
          write (31, 1100) z(iz1), z(iz2), z(iz3)
       elseif (prtz2) then
          write (31, 1200) z(iz1), z(iz2)
       elseif (prtz3) then
          write (31, 1300) z(iz1)
       endif
       ppt  = zro
       ys1  = zro
       dys1 = zro
       ys2  = zro
       dys2 = zro
       ys3  = zro
       dys3 = zro
       do k = 1,iamax
          if (dadz(k,iz1) > zro) then
             y1 = dadz4(1,k,iz1)/dadz(k,iz1)
             dy1 = sqrt(abs(dadz4(2,k,iz1)/dadz(k,iz1) - y1**2))
             ys1 = ys1 + dadz4(1,k,iz1)
             dys1= dys1+ dadz4(2,k,iz1)
          else
             y1 = zro
             dy1 = zro
          endif
          if (dadz(k,iz2) > zro) then
             y2   = dadz4(1,k,iz2)/dadz(k,iz2)
             dy2  = sqrt(abs(dadz4(2,k,iz2)/dadz(k,iz2) - y2**2))
             ys2  = ys2 + dadz4(1,k,iz2)
             dys2 = dys2 + dadz4(2,k,iz2)
          else
             y2 = zro
             dy2 = zro
          endif
          if (dadz(k,iz3) > zro) then
             y3 = dadz4(1,k,iz3)/dadz(k,iz3)
             dy3 = sqrt(abs(dadz4(2,k,iz3)/dadz(k,iz3) - y3**2))
             ys3 = ys3 + dadz4(1,k,iz3)
             dys3= dys3+ dadz4(2,k,iz3)
          else
             y3 = zro
             dy3 = zro
          endif
          printit = y1 > zro .or. y2 > zro .or. y3 > zro
          if (printit) then
             if (ar(k) > zro .and. ar(k) >= z(iz3)) then
                iaa = nint(ar(k))
                prtz1 = prtz(iz3)
                prtz2 = .not.prtz(iz3) .and. prtz(iz2)
                prtz3 = (.not.prtz(iz2) .and. .not.prtz(iz3)) .and. &
                     & prtz(iz1)
                if (prtz1) then
                   ppt = ppt + one
                   write (31, 1400) iaa, y1, dy1, y2, dy2, y3, dy3
                elseif (prtz2) then
                   ppt = ppt + one
                   write (31, 1500) iaa, y1, dy1, y2, dy2
                elseif (prtz3) then
                   ppt = ppt + one
                   write (31, 1600) iaa, y1, dy1
                endif
             endif
          endif
       end do
!  End of k loop ^
       if (ppt > zro) then
          ipt = nint(ppt)
          if (dasum(iz1) > zro) then
             ys1  = ys1/dasum(iz1)
             dys1 = sqrt(abs(dys1/dasum(iz1) - ys1**2))
          else
             ys1  = zro
             dys1 = zro
          endif
          if (dasum(iz2) > zro) then
             ys2  = ys2/dasum(iz2)
             dys2 = sqrt(abs(dys2/dasum(iz2) - ys2**2))
          else
             ys2  = zro
             dys2 = zro
          endif
          if (dasum(iz3) > zro) then
             ys3  = ys3/dasum(iz3)
             dys3 = sqrt(abs(dys3/dasum(iz3) - ys3**2))
          else
             ys3  = zro
             dys3 = zro
          endif
          if (ys3 > zro) then
             write (31, 1700) ipt, ys1, dys1, ys2, dys2, ys3, dys3
          elseif (ys2 > zro) then
             write (31, 1800) ipt, ys1, dys1, ys2, dys2
          elseif (ys1 > zro) then
             write (31, 1900) ipt, ys1, dys1
          endif
       endif
    end do
!  End of lpr loop ^
    write (31, 3200)
    ss = faca
    write (31, 3300)
    yys   = zro
    tets  = zro
    dtets = zro
    vzs   = zro
    dvzs  = zro
    yyfs  = zro
    yybs  = zro
    npt = 0
    do i = 1,ip
       yy   = dadz(i,150)*ss
       dyy  = sqrt(abs(dadz(i,150)))*ss
       yys  = yys  + dadz(i,150)
       yyfs = yyfs + dadz4(3,i,150)
       yybs = yybs + dadz4(4,i,150)
       if (dadz(i,150) > zro) then
          atet    = dadz4(1,i,151)/dadz(i,150)
          datet   = sqrt(abs(dadz4(2,i,151)/dadz(i,150) - atet**2))
          tets    = tets  + dadz4(1,i,151)
          dtets   = dtets + dadz4(2,i,151)
          avz(i)  = dadz4(1,i,150)/dadz(i,150)
          davz(i) = sqrt(abs(dadz4(2,i,150)/dadz(i,150) - avz(i)**2))
          vzs     = vzs  + dadz4(1,i,150)
          dvzs    = dvzs + dadz4(2,i,150)
          yyf     = dadz4(3,i,150)
          yyb     = dadz4(4,i,150)
          dyyf    = sqrt(abs(dadz4(3,i,150)))
          dyyb    = sqrt(abs(dadz4(4,i,150)))
!   Calculation of forward/backward ratio:
          if (yyb > zro) then
             fb(i)  = yyf/yyb
             dfb(i) = (dyyf*yyb + dyyb*yyf)/yyb**2
          else
             fb(i)  = one
             dfb(i) = zro
          endif
       else
          atet    = zro
          datet   = zro
          avz(i)  = zro
          davz(i) = zro
          fb(i)   = one
          dfb(i)  = zro
       endif
       printit = yy > zro
       if (printit) then
          npt = npt + 1
          write (31, 3400) iar(i), yy, dyy, atet, datet
       endif
    end do
!  End of i loop ^
    dyys = sqrt(abs(yys))*ss
    if (yys > zro) then
       tets  = tets/yys
       dtets = sqrt(abs(dtets/yys - tets**2))
       vzs   = vzs/yys
       dvzs  = sqrt(abs(dvzs/yys - vzs**2))
    endif
    if (yybs > zro) then
       fbs   = yyfs/yybs
       dyyfs = sqrt(abs(yyfs))
       dyybs = sqrt(abs(yybs))
       dfbs  = (dyyfs*yybs + dyybs*yyfs)/yybs**2
    else
       fbs  = one
       dfbs = zro
    endif
    yys = yys*ss
    write (31, 3500) npt, yys, dyys, tets, dtets
    write (31, 3350)
    do i = 1,ip
       yy   = dadz(i,150)*ss
       printit = yy > zro
       if (printit) write (31, 3400) iar(i), avz(i), davz(i), fb(i), &
            & dfb(i)
    end do
    write (31, 3500) npt, vzs, dvzs, fbs, dfbs
    write (31, 3600)
    yys   = zro
    tets  = zro
    dtets = zro
    vzs   = zro
    dvzs  = zro
    yyfs  = zro
    yybs  = zro
    npt = 0
    do i = 1,ipz+1
       yy   = dadz(351,i)*ss
       dyy  = sqrt(abs(dadz(351,i)))*ss
       yys  = yys  + dadz(351,i)
       yyfs = yyfs + dadz4(3,351,i)
       yybs = yybs + dadz4(4,351,i)
       if (dadz(351,i) > zro) then
          atet    = dadz4(1,351,i)/dadz(351,i)
          datet   = sqrt(abs(dadz4(2,351,i)/dadz(351,i) - atet**2))
          tets    = tets + dadz4(1,351,i)
          dtets   = dtets+ dadz4(2,351,i)
          avz(i)  = dadz4(1,350,i)/dadz(351,i)
          davz(i) = sqrt(abs(dadz4(2,350,i)/dadz(351,i) - avz(i)**2))
          vzs     = vzs + dadz4(1,350,i)
          dvzs    = dvzs+ dadz4(2,350,i)
          yyf     = dadz4(3,350,i)
          yyb     = dadz4(4,350,i)
          dyyf    = sqrt(abs(dadz4(3,350,i)))
          dyyb    = sqrt(abs(dadz4(4,350,i)))
!   Calculation of forward/backward ratio:
          if (yyb > zro) then
             fb(i)  = yyf/yyb
             dfb(i) = (dyyf*yyb + dyyb*yyf)/yyb**2
          else
             fb(i)  = one
             dfb(i) = zro
          endif
       else
          atet    = zro
          datet   = zro
          avz(i)  = zro
          davz(i) = zro
          fb(i)   = one
          dfb(i)  = zro
       endif
       printit = yy > zro
       if (printit) then
          npt = npt + 1
          write (31, 3700) izr(i), yy, dyy, atet, datet
       endif
    end do
!  End of i loop ^
    dyys = sqrt(abs(yys))*ss
    if (yys > zro) then
       tets  = tets/yys
       dtets = sqrt(abs(dtets/yys - tets**2))
       vzs   = vzs/yys
       dvzs  = sqrt(abs(dvzs/yys - vzs**2))
    endif
    if (yybs > zro) then
       fbs   = yyfs/yybs
       dyyfs = sqrt(abs(yyfs))
       dyybs = sqrt(abs(yybs))
       dfbs  = (dyyfs*yybs + dyybs*yyfs)/yybs**2
    else
       fbs  = one
       dfbs = zro
    endif
    yys = yys*ss
    write (31, 3800) npt, yys, dyys, tets, dtets
    write (31, 3350)
    do i = 1,ipz+1
       yy   = dadz(351,i)*ss
       printit = yy > zro
       if (printit) write (31, 3700) izr(i), avz(i), davz(i), fb(i), &
            & dfb(i)
    end do
    write (31, 3800) npt, vzs, dvzs, fbs, dfbs
    return

! ======================================================================
1000 format (/1x,'********** Nuclide yields [mb] in forward ', &
          &        'direction (theta_lab < 90) **********'/11x, &
          &        '(zero values suppressed)')
1100 format (/17x,'Z = ',f4.0,15x,'Z = ',f4.0,15x,'Z = ',f4.0)
1200 format (/17x,'Z = ',f4.0,15x,'Z = ',f4.0)
1300 format (/17x,'Z = ',f4.0)
1400 format (1x,'A =',i4,3(1x,1pe9.3,' +/- ',e8.2))
1500 format (1x,'A =',i4,2(1x,1pe9.3,' +/- ',e8.2))
1600 format (1x,'A =',i4,1(1x,1pe9.3,' +/- ',e8.2))
1700 format (1x,'S =',i4,3(1x,1pe9.3,' +/- ',e8.2))
1800 format (1x,'S =',i4,2(1x,1pe9.3,' +/- ',e8.2))
1900 format (1x,'S =',i4,1(1x,1pe9.3,' +/- ',e8.2))
2000 format (/1x,'End of nuclide yields (forward direction).')
2100 format (/' Mass yield [mb] and the mean and variance of the ', &
          &         'kinetic energy [MeV]'/1x,'of residual nuclei in the ', &
          &         'forward direction:')
2200 format (1x,'A =',i4,1x,1pe9.3,' +/- ',e8.2,2x,e9.3,' +/- ',e8.2)
2300 format (1x,'S =',i4,1x,1pe9.3,' +/- ',e8.2,2x,e9.3,' +/- ',e8.2)
2400 format (/'Charge yield [mb] and the mean and variance of the ', &
          &         'kinetic energy [MeV]'/1x,'of residual nuclei in the ', &
          &         'forward direction:')
2500 format (1x,'Z =',i3,1x,1pe9.3,' +/- ',e8.2,2x,e9.3,' +/- ',e8.2)
2600 format (1x,'S =',i3,1x,1pe9.3,' +/- ',e8.2,2x,e9.3,' +/- ',e8.2)
2700 format (/1x,'********** Nuclide yields [mb] in backward ', &
          &        'direction (theta_lab > 90) *********'/11x, &
          &        '(zero values suppressed)')
2800 format (/1x,'End of nuclide yields (backward direction).')
2900 format (/' Mass yield [mb] and the mean and variance of the ', &
          &         'kinetic energy [MeV]'/1x,'of residual nuclei in the ', &
          &         'backward direction:')
3000 format (/' Charge yield [mb] and the mean and variance of the ', &
          &         'kinetic energy [MeV]'/1x,'of residual nuclei in the ', &
          &         'backward direction:')
3100 format (/'******** Nuclide average kinetic energies ', &
          &    '[MeV] (zero yield suppressed) ********')
3200 format (/1x,'End of nuclide average kinetic energies.')
3300 format (/' Mass yield [mb] and the mean and variance of the ', &
          &         'emission angle [deg.]'/1x,'of residual nuclei:')
3350 format (/' The mean and variance of the z velocity [v/c] of ', &
          &         'residual nuclei,'/1x,'and the forward/backward ratio:')
3400 format (1x,'A =',i4,1x,1pe10.3,' +/- ',e8.2,2(2x,e10.3,' +/- ', &
          &        e8.2))
3500 format (1x,'S =',i4,1x,1pe10.3,' +/- ',e8.2,2(2x,e10.3,' +/- ', &
          &        e8.2))
3600 format (/' Charge yield [mb] and the mean and variance of the ', &
          &         'emission angle [deg.]'/1x,'of residual nuclei:')
3700 format (1x,'Z =',i3,1x,1pe10.3,' +/- ',e8.2,2(2x,e10.3,' +/- ', &
          &        e8.2))
3800 format (1x,'S =',i3,1x,1pe10.3,' +/- ',e8.2,2(2x,e10.3,' +/- ', &
          &        e8.2))
5000 format("Divide by zero error prevented in 'prtdadz.f90', line(s) ", A)
! ======================================================================
  end subroutine prtdadz
