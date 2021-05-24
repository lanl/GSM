
  function rndm (rdummy)

! ======================================================================
!                                                                      *
!     This routine merely acts as an interface to F.James rm48  gen.   *
!                                                                      *
!     Created on 03  April  1992   by    Alfredo Ferrari & Paola Sala  *
!            INFN - Milan                                              *
!                                                                      *
!     Modified    on 16-Sep-93     by    Alfredo Ferrari               *
!                                                                      *
!     Revision: 13-FEB-1997    BY NVMokhov                             *
!     Edited by A. J. Sierk, LANL T-16, September-October, 2003.       *
!                                                                      *
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, int64, real64
    use randomNumberGenerator, only: rndmRANG => rang

!    implicit none
    implicit real(real64) (a-h, o-z)
    integer(int32) ::  cntr, idummy, ione, iseed1, iseed2

    real(real64) ::  rdummy
! ======================================================================

    parameter (pipipi = 3.1415926535897932270d0)

    dimension rndnum (2)
    common /countrn/ cntr

    abstract interface
       function RANDOM() result(rndm) BIND(C)
         use, intrinsic:: iso_C_binding, only: c_double
         implicit none
         real(c_double) :: rndm
       end function RANDOM
    end interface
    procedure(RANDOM), pointer, save :: rang => rndmRANG

    logical, parameter :: useGsmRng = .TRUE.

! ======================================================================

    ! Obtain Random Number from MCNP5 RNG and return control
    if ( useGsmRng ) then

       rndm = rang()
       return

    endif
! ----------------------------------------------------------------------

    call rm48 (rndnum, 1_int32)
    rndm = rndnum (1)
    cntr = cntr + 1
    return

    entry rd2in (iseed1, iseed2)
!  The following card just to avoid warning messages on the HP compiler
    rd2in = pipipi
    call rm48in (54217137, iseed1, iseed2)
    return

    entry rd2out (iseed1, iseed2)
!  The following card just to avoid warning messages on the HP compiler
    rd2out = pipipi
    call rm48ut (idummy, iseed1, iseed2)
    return

! ======================================================================
1000 format (3x, "comment: using the ", a, " random number generator ", &
         & "with the dcm.",/)
! ======================================================================
  end function rndm

  subroutine rm48 (rvec, lenv)

! ======================================================================
!
!     Double-precision version of
!    universal random number generator proposed by Marsaglia and Zaman
!    in report FSU-SCRI-87-50
!      Based on RANMAR, modified by F. James, to generate vectors
!      of pseudorandom numbers rvec of length lenv, where the numbers
!      in rvec are numbers with at least 48-bit mantissas.
!                                                                      *
!     Revision: 21-JAN-1997    BY NVMokhov                             *
!     Modified by A. J. Sierk, LANL T-16, September-October, 2003.
!                                                                      *
!   Input and output entry points: rm48in, rm48ut.
! !! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! !!  Calling sequences for rm48 :                                    ++
! !!      call rm48  (rvec, len)    returns a vector rvec of len      ++
! !!      64-bit random floating point numbers between                ++
! !!      zero and one.                                               ++
! !!      call rm48in (i1,n1,n2)  initializes the generator from one  ++
! !!      64-bit integer i1, and number counts n1,n2                  ++
! !!     (for initializing, set n1=n2=0, but to restart               ++
! !!       a previously generated sequence, use values                ++
! !!       output by rm48ut)                                          ++
! !!      call rm48ut (i1,n1,n2)  outputs the value of the original   ++
! !!     seed and the two number counts, to be used                   ++
! !!     for restarting by initializing to I1 and                     ++
! !!     skipping n2*100000000+n1 numbers.                            ++
! !! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit real(real64) (a-h, o-z), integer(int32) (i-n)

    character(len=15) :: chars1
    character(len=18) :: chars2
    character(len=22) :: chars3

! ======================================================================

    parameter (modcns=1000000000)

    dimension rvec(*)

    common /r48st1/ u(97), c, i97, j97

    save cd, cm, twom24,  zero, one, half, ntot, ntot2, ijkl

    data ntot, ntot2, ijkl /-1, 0, 0/
    data lunout /16/
    data one, half, zero /1.d0, 0.5d0, 0.d0/

! ======================================================================

    if (ntot >= 0) go to 20
!
!   Default initialization. User has called rm48  without rm48in.
    ijkl = 54217137
    ntot = 0
    ntot2 = 0
    lenv = max (lenv, 1)
    kalled = 0
    go to 10
!
    entry rm48in (ijklin, ntotin, ntot2n)
!    Initializing routine for rm48 , may be called before
!    generating pseudorandom numbers with rm48 .   The input
!    values should be in the ranges:  0<=ijklin<=900 OOO OOO
!           0<=ntotin<=999 999 999
!           0<=ntot2n<<999 999 999!
! To get the standard values in Marsaglia's paper, ijklin=54217137
!     ntotin,ntot2n=0
    ijkl = ijklin
    ntot = max(ntotin, 0)
    ntot2= max(ntot2n, 0)
    lenv = 1
    kalled = 1
!   Always come here to initialize:
10  continue
    ij = ijkl/30082
    kl = ijkl - 30082*ij
    i = mod(ij/177, 177) + 2
    j = mod(ij, 177)     + 2
    k = mod(kl/169, 178) + 1
    l = mod(kl, 169)
    chars2 = ' rm48 initialized:'
!     write (lunout, '(a,i10,2x,2i10)')
!    &       ' RM48  initialized:', ijkl, ntot, ntot2
    write (lunout, 1000) chars2, ijkl, ntot, ntot2
    do ii = 1, 97
       x = zero
       t = half
       do jj = 1, 48
          m = mod(mod(i*j,179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l + 1, 169)
          if (mod(l*m,64) >= 32)  x = x + t
          t = half*t
       end do
       u(ii) = x
    end do
    twom24 = one
    do i24 = 1, 24
       twom24 = half*twom24
    end do
    c  =   362436.d0*twom24
    cd =  7654321.d0*twom24
    cm = 16777213.d0*twom24
    i97 = 97
    j97 = 33
!   Complete initialization by skipping
!            (ntot2*modcns + ntot) random numbers
    do loop2 = 1, ntot2+1
       now = modcns
       if (loop2 == ntot2 + 1) now = ntot
       if (now > 0) then
          chars3 = ' rm48in skipping over '
!         write (lunout,'(a,i15)') ' RM48IN skipping over ', now
          write (lunout, 1400) chars3, now
          do idum = 1, ntot
             uni = u(i97) - u(j97)
             if (uni < zero) uni = uni + one
             u(i97) = uni
             i97 = i97 - 1
             if (i97 == 0) i97 = 97
             j97 = j97 - 1
             if (j97 == 0) j97 = 97
             c = c - cd
             if (c < zero) c = c + cm
          end do
       endif
    end do
    if (kalled == 1)  return
!
!    Normal entry to generate lenv random numbers:
20  continue
    do ivec = 1, lenv
       uni = u(i97) - u(j97)
       if (uni < zero) uni = uni + one
       u(i97) = uni
       i97 = i97 - 1
       if (i97 == 0) i97 = 97
       j97 = j97 - 1
       if (j97 == 0) j97 = 97
       c = c - cd
       if (c < zero) c = c + cm
       uni = uni - c
       if (uni < zero) uni = uni + one
       rvec(ivec) = uni
    end do
    ntot = ntot + lenv
    if (ntot >= modcns) then
       ntot2 = ntot2 + 1
       ntot = ntot - modcns
    endif
    return

!  Entry to output current status
    entry rm48ut (ijklut, ntotut, ntot2t)
    ijklut = ijkl
    ntotut = ntot
    ntot2t = ntot2
    chars1 = ' rm48ut output:'
!     write (lunout, '(//a,i10,2x,2i10)')
!    & ' RM48ut output:', ijkl, ntot, ntot2
!     write (lunout, 1000) chars1, ijkl, ntot, ntot2
    return

!   Output routine for RM48 , without skipping numbers:
    entry rm48wr (ioseed)
!     write (ioseed, '(2i10)') ntot, ntot2
!     write (ioseed, '(2i10,f18.16)') i97, j97, c
!     write (ioseed, '(24(4z16,/),z16)') u
    write (ioseed, 1100) ntot, ntot2
    write (ioseed, 1200) i97, j97, c
    write (ioseed, 1300) u
    return

!   Initializing routine for RM48 , without skipping numbers:
    entry  rm48rd (ioseed)
!     read (ioseed, '(2i10)') ntot, ntot2
!     read (ioseed, '(2i10,f18.16)') i97, j97, c
!     read (ioseed, '(24(4z16,/),z16)') u
    write (ioseed, 1100) ntot, ntot2
    write (ioseed, 1200) i97, j97, c
    write (ioseed, 1300) u
    close (unit = ioseed)
    ijkl = 54217137
    ij = ijkl/30082
    kl = ijkl - 30082*ij
    i = mod(ij/177, 177) + 2
    j = mod(ij, 177)     + 2
    k = mod(kl/169, 178) + 1
    l = mod(kl, 169)
    chars2 = ' rm48 initialized:'
!     write (lunout,'(a,i10,2x,2i10)')
!    &  ' RM48  initialized:', ijkl, ntot, ntot2
    write (lunout,1000) chars2, ijkl, ntot, ntot2
    twom24 = one
    do i24 = 1, 24
       twom24 = half*twom24
    end do
    cd =  7654321.d0*twom24
    cm = 16777213.d0*twom24
    return

! ======================================================================

1000 format (//a18,i10,2x,2i10)
1100 format (2i10)
1200 format (2i10,f18.16)
1300 format (24(4z16,/),z16)
1400 format (a22,i15)

! ======================================================================
  end subroutine rm48
! ***************************************************************

  subroutine rdmini
!   Last change: 17-DEC-2003 by NVM
!   Last Change: 07-SEP-2017 by CMJ (XCP-3 LANL) - added call to RDMIN

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit integer(int32) (i-n)

    iseed1=54217137
    iseed2=0
    iseed3=0
    call rdmin (iseed1,iseed2,iseed3)
!      CALL RM48IN (ISEED1,ISEED2,ISEED3)
    return
  end subroutine rdmini
! ----------------------------------------------------------------
! ******************************************
  subroutine rdmin (iseed1,iseed2,iseed3)

    use, intrinsic:: iso_fortran_env, only: int32, real64

!    implicit none
    implicit integer(int32) (i-n)
    if(iseed1 < 0.or.iseed1 > 900000000)  iseed1=54217137
    if(iseed2 < 0.or.iseed2 > 999999999)  iseed1=0
    if(iseed3 < 0.or.iseed3 > 900000000)  iseed1=0
    call rm48in (iseed1,iseed2,iseed3)
!     write(*,*) ' ISEED1,2,3=',ISEED1,ISEED2,ISEED3
    return
! ******************************************
    entry rdmout(iseed1,iseed2,iseed3)
    call rm48ut (iseed1,iseed2,iseed3)
    return
  end subroutine rdmin
