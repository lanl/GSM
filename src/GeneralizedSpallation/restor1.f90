
  subroutine restor1 (gsmObj, numBaryons, numProtons, recoilEnergy, restMass, v, results)

! ======================================================================
!
!    Storing the properties of the residual nucleus in the spt and parz
!    arrays for nuclei which do not evaporate (low E*)
!
!    Calls:  ATAN2, SIGN
!
!    Written by A. J. Sierk, LANL T-16, March, 2004.
!    Edited by AJS, LANL T-2, December, 2011.
!    Edited by LMK, XCP-3, uly 2013 (included error protection)
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use gsm_params, only : zro, one, two, thr, four, fiv, six, thousand, &
         & pi, twpi

    implicit none
    class(GSM),     intent(inout) :: gsmObj
    real(real64),   intent(in   ) :: numBaryons
    real(real64),   intent(in   ) :: numProtons
    real(real64),   intent(in   ) :: recoilEnergy
    real(real64),   intent(in   ) :: restMass
    real(real64),   intent(in   ) :: v(3)
    type(GSMResults), intent(inout) :: results

    integer(int32) :: inj, izj, particleID
    real(real64)   :: absv, cosPhi, cosTheta, phi, sinPhi, sinTheta, theta, &
         & temp, temp2

! ======================================================================

    ! Ensure bank isn't full:
    if ( results%numProgeny > results%maxProgenyM1 ) then
       write(gsmObj%io%message, 1000)
       call gsmObj%io%print(3, 3, gsmObj%io%message)
    else
       ! Obtain interim results (emission angles namely)
       temp = v(1)**2 + v(2)**2 + v(3)**2
       absv = sqrt (temp)
       temp2 = absv
       if ( abs(temp2) < div0Lim ) then
          temp2 = div0Lim
          write(gsmObj%io%message,2000) "456"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       cosTheta = v(3)/temp2
       if (abs(cosTheta) >= one) then
          cosTheta = sign(one,cosTheta)
          sinTheta = zro
          theta = zro
          if (cosTheta < zro) theta = pi 
       else 
          sinTheta = sqrt(abs(one - cosTheta**2))
          theta = atan2(sinTheta,cosTheta)
       endif
       if (abs(sinTheta) > 1.d-10) then
          cosPhi = v(1)/(absv*sinTheta)
          sinPhi = v(2)/(absv*sinTheta)
          if (abs(cosPhi) >= one) then
             cosPhi = sign(one,cosPhi)
             if (cosPhi <= zro) phi = pi
             if (cosPhi > zro) phi = zro
             sinPhi = zro
          else
             phi = atan2(sinPhi,cosPhi)
             if (phi < zro) phi = twpi + phi
          endif
       else
          phi = zro
          cosPhi = one
          sinPhi = zro
       endif

       ! Obtain particle ID:
       izj = nint(numProtons)
       inj = nint(numBaryons) - izj
       particleID = 0
       if (izj == 0 .and. inj == 1) particleID = one 
       if (izj == 1 .and. inj == 0) particleID = two
       if (izj == 1 .and. inj == 1) particleID = thr
       if (izj == 1 .and. inj == 2) particleID = four
       if (izj == 2 .and. inj == 1) particleID = fiv
       if (izj == 2 .and. inj == 2) particleID = six
       if ( particleID == 0 ) particleID = 1000*izj + inj


       ! Tally particle:
       results%numProgeny = results%numProgeny + 1
       results%progenyBnk(results%numProgeny)%numBaryons = numBaryons
       results%progenyBnk(results%numProgeny)%numProtons = numProtons
       results%progenyBnk(results%numProgeny)%kinEnergy  = recoilEnergy / thousand
       results%progenyBnk(results%numProgeny)%restMass   = restMass
       results%progenyBnk(results%numProgeny)%phi        = phi
       results%progenyBnk(results%numProgeny)%theta      = theta
       results%progenyBnk(results%numProgeny)%sinTheta   = sinTheta
       results%progenyBnk(results%numProgeny)%cosTheta   = cosTheta
       results%progenyBnk(results%numProgeny)%typeID     = particleID
       results%progenyBnk(results%numProgeny)%prodMech   = 1000
    end if

    return
! ======================================================================
1000 format("The GSM progeny bank is full. The stable nucleus cannot be tallied.")
2000 format("Divide by zero error prevented in 'restor1.f90', line ", A)
! ======================================================================
  end subroutine restor1
