
  subroutine fermid (fbuObj, nucleus, results)

! ======================================================================
!
!     Fermi breakup driver routine (executes fermi breakup simulation)
!
!     Produced fragments are in prod: ( no longer the case when switch to fermi_bnk() made)
!          i                          prod(i,k)
!     --------------------------------------------------------
!          1                         mass number
!          2                         charge
!          3                         Px (MeV/c)
!          4                         Py (MeV/c) Momentum
!          5                         Pz (MeV/c)
!          6                         Ek (MeV), kinetic energy
!          7                         M  (MeV), mass
!     --------------------------------------------------------
!     mv = nst - nst0 - number(multiplicity) of produced  fragments
!     --------------------------------------------------------
!
! ======================================================================
!
!    Called by: execute() [from class definition; execution in GSM in fermib]
!
!    Calls: RASTARFERMI
!
!    "Last" change: 13-Aug-2003 BY NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Modified 03 May 2006 by R E Prael to allow decay of multiple 
!       intermediate states. The only interediate states currently
!       allowed are He-5, Be-8, and B-9. Others may be added by defining
!       MS1 in subroutine GITAB.
!    Edited by A. J. Sierk, LANL T-2, February, 2009.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by CMJ, LANL XCP-3, July 2018 (creation of FermiBreakup class)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use fermiBreakupParams, only : zro, one, thousand

    implicit none
    class(FermiBreakup),       intent(inout) :: fbuObj
    type(fermiBreakUpNucleus), intent(in   ) :: nucleus
    type(fermiBreakUpResults), intent(inout) :: results

    real(real64), dimension(3) :: pf
    type(fermiBreakUpProgeny), dimension(results%maxProgeny) :: fermiBnk0

    integer(int32) ::  iaf, izf, ksf, m, m1, mf, mv, &
         & nst, nst0, nst1, ntry
    real(real64)   :: af, zf, uf

! ======================================================================

    ! ---------------------
    ! Start of simulation
    ! ---------------------
    nst = zro ! Assume no fragments at this time
    nst0 = nst
    nst1 = zro
    ! Use temporary array for rastar call:
    pf(1) = nucleus%linearXMom
    pf(2) = nucleus%linearYMom
    pf(3) = nucleus%linearZMom
    call fbuObj%rastar (nucleus%numBaryons, nucleus%numProtons, &
         & nucleus%kinEnergy, pf, nst1, 1, results)
    mv = nst1

    ! Initialize 'fermiBnk0'
    if ( mv > results%maxProgeny ) then
       ! Ensure arrays for progeny storage aren't exceeded
       write(fbuObj%io%message,1000)
       call fbuObj%io%print(2, 3, fbuObj%io%message)
       write(fbuObj%io%message,1010) (mv - results%maxProgeny)
       call fbuObj%io%print(2, 3, fbuObj%io%message)
       mv = results%maxProgeny
       results%simState = results%simState + exceededProgenyArray
    end if
    do m = 1,mv
       fermiBnk0(m) = results%progenyBnk(m)
    end do

10  continue
    ntry = 0
    m1 = mv
    do m = 1,mv
       af = fermiBnk0(m)%numBaryons
       zf = fermiBnk0(m)%numProtons
       iaf = nint(af)
       izf = nint(zf)
! All allowed intermediate states should be included below. More are
! included here than are presently allowed.
       if ((iaf == 5 .and. (izf == 2.or.izf == 3)) .or.    &
            (iaf == 6 .and. izf >= 4) .or.    &
            (iaf == 7 .and. izf >= 2) .or.    &
            (iaf == 8 .and. (izf == 4.or.izf == 6)) .or.    &
            (iaf == 9 .and. izf == 5)) then
          pf(1) = fermiBnk0(m)%linearXMom/thousand
          pf(2) = fermiBnk0(m)%linearYMom/thousand
          pf(3) = fermiBnk0(m)%linearZMom/thousand
          uf = 1.0d-6
          ksf = 0

          call fbuObj%rastar (af, zf, uf, pf, ksf, 1, results)
          if (ksf > 1) ntry = 1
          do mf = 1,ksf
             if (mf == 1) then
                ! Obtain first fragment
                fermiBnk0(m) = results%progenyBnk(mf)
             else
                ! Add fragment to list
                m1 = m1 + 1
                if ( m1 > results%maxProgeny ) then
                   ! If all fragments have been stored, print warning and continue in simulation
                   write(fbuObj%io%message,1000) (m1 - results%maxProgeny)
                   call fbuObj%io%print(2, 3, fbuObj%io%message)
                   results%simState = results%simState + exceededProgenyArray
                   m1 = m1 - 1
                   go to 20
                end if
                fermiBnk0(m1) = results%progenyBnk(mf)
             endif
          end do
20        continue
       endif
    end do
    mv = m1
    if (ntry > 0) go to 10
    nst = nst0 + mv
    ! Store produced fragments into the tallied array
    if ( mv > results%maxProgeny ) then
       write(fbuObj%io%message,1000) (mv - results%maxProgeny)
       call fbuObj%io%print(2, 3, fbuObj%io%message)
       mv = results%maxProgeny
    end if

    ! Transfer progeny over for final progeny bank:
    do m = 1,mv
       results%progenyBnk(m) = fermiBnk0(m)
    end do

    ! Store number of produced fragments
    results%numProgeny = nst

    return
! ====================================================================== 
1000 format("The fragment array during Fermi Break-Up was exceeded.")
1010 format("   The remaining ", i3, " fragment(s) will not be stored.")
! ====================================================================== 
  end subroutine fermid
