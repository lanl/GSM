
  subroutine setPhotoChannelSigma (sDCM, egamma, photoData)

! ======================================================================
!    Main routine to extract ds/do for channel 1-4 (single pion
!    production):
!
!    Written by K. K. Gudima, Fall 2003?
!    Modified by AJS, May, 2004.
!    Modified by KKG, Sep., 2004
!    Edited by AJS, January, 2005.
!    Modified by AJS, February, 2005.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use standardDCMParams, only: zro, one, two, radianToDegree, emnucg
    use standardDCMData,   only: domo2, domo2i, gppipn, gppi0p, elg, &
         & thetai, photonEG

    implicit none
    class(StandardDCM),            intent(inout) :: sDCM
    real(real64),                  intent(in   ) :: egamma
    type(sDCMPhotonCrossSections), intent(  out) :: photoData

    integer(int32) :: ie, ieg1, ieg2, j, jch
    real(real64)   :: deleg, den1, eg1, eg2, sig1, sig2, sint, temp

    real(real64), dimension(2, 20) :: sig = zro
    real(real64), dimension(2, 19) :: r   = zro

! ======================================================================

    ! Obtain cross section for gamma energy in each channel at each angle
    do jch = 1,2

       ! Obtain ieg1/2 values for nearest cross section data
       if (egamma <= elg(jch,2)) then
          ieg1 = 2
          ieg2 = 3
       elseif (egamma >= elg(jch,50)) then
          ieg1 = 49
          ieg2 = 50
       else
          do ie = 3,50
             if (egamma >= elg(jch,ie-1) .and. egamma <= elg(jch,ie)) then
                ieg1 = ie - 1
                ieg2 = ie
                go to 10
             endif
          enddo
       endif

10     continue

       ! Check ieg1/2 values for validity (i.e. within range)
       if (ieg1 < 2 .or. ieg1 > 49 .or. ieg2 < 3 .or. ieg2 > 50) then
          write(sDCM%io%message,2000) ieg1, ieg2
          call sDCM%io%print(3, 3, sDCM%io%message)
          write(sDCM%io%message,2010)
          call sDCM%io%print(3, 3, sDCM%io%message)

          ! Obtain nearest index
          if ( ieg1 <  2 ) ieg1 =  2
          if ( ieg1 > 49 ) ieg1 = 49
          if ( ieg2 <  3 ) ieg2 =  3
          if ( ieg2 > 50 ) ieg2 = 50
       endif

       ! Obtain interpolation energies for the channel
       eg1 = elg(jch,ieg1)
       eg2 = elg(jch,ieg2)
       temp = eg2 - eg1
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(sDCM%io%message,1000) "116"
          call sDCM%io%print(4, 3, sDCM%io%message)
       end if

       ! Obtain interpolation values (den1 and deleg)
       den1 = zro
       if (eg2.ne.eg1) den1 = one/(temp)
       deleg = egamma - eg1
       sint = zro

       ! Obtain cross section for all bins
       do j = 1,19

          ! Obtain points nearest data for the channel
          if (jch == 1)     then
             sig1 = gppipn(ieg1,j)
             sig2 = gppipn(ieg2,j)
          elseif (jch == 2) then
             sig1 = gppi0p(ieg1,j)
             sig2 = gppi0p(ieg2,j)
          endif

          ! Perform linear interpolation
          sig(jch,j) = sig1 + (sig2 - sig1)*(den1*deleg)

          ! Calculate the channel's angle-integrated cross section
          if (j >= 2) then
             sint = sint + domo2(j)*(sig(jch,j-1) + sig(jch,j))
          endif
          r(jch,j) = sint

       enddo
       sig(jch,20) = sint


       ! Normalize 'r" for all bins in the channel
       do j = 1,19
          if (sint > zro) r(jch,j) = r(jch,j)/sint
       enddo

    end do


    ! Obtain "fine" mesh cross section based on the "course" mesh cross sections
    do jch = 1,2
       sint = zro
       do j = 1,181
          ! Interpolate for angle "j" based on course mesh (theta, sig)
          photoData%si(jch,j) = photonEG%qintxs(thetai(j), sig, jch, 2, 19)

          ! Calculate angle-integrated cross section
          if (j >= 2) then
             sint = sint + domo2i(j)*(photoData%si(jch,j-1) + photoData%si(jch,j))
          endif
          photoData%ri(jch,j) = sint
       end do
       photoData%si(jch,182) = sint

       ! Normalize "ri" for each angle
       do j = 1,181
          if (sint > zro) photoData%ri(jch,j) = photoData%ri(jch,j)/sint
       enddo

    end do

    return

! ======================================================================
1000 format("Divide by zero error prevented in ", &
          & "'setPhotoChannelSigma.f90', line(s)", A)
2000 format("Photon cross section interpolation bounds (ieg1=", i3, &
          & ", ieg2=", i3, ") could not be obtained.")
2010 format("   Will round to nearest valid bound...")
! ======================================================================
  end subroutine setPhotoChannelSigma
