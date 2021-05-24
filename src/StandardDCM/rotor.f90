
  subroutine rotor (sDCM, ar, br, pstar, pr)

! ======================================================================
!
!    Coordinate rotation.
!
!    Called by: ABEL DIRECT8 ISOBAR PRECOF STAT
!
!    ar is the momentum of the emitting system in the CM frame:
!    it defines the z' axis.
!    br is the velocity vector of the CM frame in the lab frame.
!    ar and br define a plane, which contains the x' axis.
!    pstar is the momentum of one of the final state particles in the
!    CM' frame.
! 
!    This frame is defined with z' parallel to ar, x' perpendicular
!    to ar, in the plane defined by ar and br, and y' perpendicular
!    to the same plane.  The routine calculates the vector pr, which
!    is the vector of final state particle momentum in the CM frame 
!    which has its coordinate axes parallel to those of the lab frame,
!    so pr is easily Lorentz transformed into the lab frame.
!
!   CEM95 written by S. G. Mashnik
!
!   Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!   Modified by A. J. Sierk,  LANL  T-2  March, 1996.
!   Further editing by AJS, July, 1997.
!   Modified by AJS; August, 1997.
!   Comments corrected by AJS, December, 1999.
!   Modified for real*8 by A. J. Sierk, LANL T-16, October, 2003.
!   Edited by AJS, LANL T-2, December, 2011.
!   Edited by LMK, XCP-3, July 2013 (included error protection)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use standardDCMParams, only : zro, one

    implicit none
    class(StandardDCM), intent(inout) :: sDCM
    real(real64), dimension(3), intent(in)  :: ar
    real(real64), dimension(3), intent(in)  :: br
    real(real64), dimension(3), intent(in)  :: pstar
    real(real64), dimension(3), intent(out) :: pr

    integer(int32) :: ir
    real(real64)   :: alpha2, amod, bmod, cmod, cp, dmod, f2, fac
    real(real64), dimension(3) :: an=zro, bn=zro, cn=zro, dn=zro

! ======================================================================

    amod = sqrt(ar(1)**2 + ar(2)**2 + ar(3)**2)
    bmod = sqrt(br(1)**2 + br(2)**2 + br(3)**2)
    if (amod < div0Lim .and. amod > -div0Lim) then
       amod = div0Lim
       write(sDCM%io%message,1000) "62, 66"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    if (bmod < div0Lim .and. bmod > -div0Lim) then
       bmod = div0Lim
       write(sDCM%io%message,1000) "63, 66"
       call sDCM%io%print(4, 3, sDCM%io%message)
    end if
    cp = zro
    do ir = 1,3
       cp = cp + ar(ir)*br(ir)
       an(ir) = ar(ir)/amod
       bn(ir) = br(ir)/bmod
    end do
!   cp is the cosine of the angle between vectors ar and br.
    cp = cp/(amod*bmod)
    if (abs(cp) > one) cp = sign(one, cp)
!   cn is the cross product of two unit vectors; its magnitude is
!   sin (angle);  redefined cn is the unit vector perpendicular to 
!   the plane of ar and br.
    cn(1) = an(2)*bn(3) - an(3)*bn(2)
    cn(2) = an(3)*bn(1) - an(1)*bn(3)
    cn(3) = an(1)*bn(2) - an(2)*bn(1)
    cmod = sqrt(cn(1)**2 + cn(2)**2 + cn(3)**2)
    if (cmod > zro) then
!  The cross product is not exactly 0 to machine precision:
       if (cmod < div0Lim .and. cmod > -div0Lim) then
          cmod = div0Lim
          write(sDCM%io%message,1000) "81"
       call sDCM%io%print(4, 3, sDCM%io%message)
       end if
       fac = one/cmod
       do ir = 1,3
          cn(ir) = cn(ir)*fac
          dn(ir) = bn(ir) - cp*an(ir)
       end do
    else
!   bn is exactly (or anti-) parallel to an to within machine precision:
!   Define an arbitrary coordinate system with an as its z' axis.
       cn(2) = zro
       if (an(3).ne.zro) then
          fac = one/sqrt(an(1)**2 + an(3)**2)   !No check needed because of if statement
          cn(1) =  fac*an(3)
          cn(3) = -fac*an(1)
          if (an(2).ne.zro) then
             dn(1) = an(1)/an(3)   !No check needed because of if statement
             dn(2) = -(an(1)**2 + an(3)**2)/(an(2)*an(3))
             dn(3) = one
          else
             dn(1) = zro
             dn(2) = -one
             dn(3) = zro
          endif
       else
          cn(1) =  zro
          cn(3) = -one
          dn(1) =  an(2)
          dn(2) = -an(1)
          dn(3) =  zro
       endif
    endif
    alpha2 = sqrt(one - cp**2)
    if (alpha2 > 1.0d-6) then
       do ir = 1,3
          pr(ir) = pstar(1)*dn(ir)/alpha2 + pstar(3)*an(ir) &
               + pstar(2)*cn(ir)
       end do
    else
!  an and bn are nearly parallel or antiparallel; round-off may make
!  dmod = zro to computer word length.
       dmod = sqrt(dn(1)**2 + dn(2)**2 + dn(3)**2)
       if (dmod.ne.zro) then
          f2 = one/dmod
       else
          f2 = zro
       endif
       do ir = 1,3
          pr(ir) = pstar(3)*an(ir) + pstar(2)*cn(ir) + &
               pstar(1)*dn(ir)*f2
       end do
    endif

    return

! ======================================================================
1000 format ("Divide by zero error prevented in 'rotor.f90', line(s) ", A)
! ======================================================================
  
end subroutine rotor

