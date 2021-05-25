
  subroutine coalescenceInterface ( gsmObj, rxnPhysics, results, outData )

! ======================================================================
!
! Interface routine for GSM to the coalescence library.
!
! Interface stores main array values into coalescence library data
! structures, then simulates and moves resulting fragments (from data
! structures) to the main GSM fragment arrays (parz/spt)
!
!
! Written (prior to class-based model) by KKG (setting and getting of partBnk members [previously in arrays])
! Modifed by CMJ, XCP-3, July 2018 (Coalscence Class creation).
! Modifed by CMJ, XCP-3, Jan. 2019 (Modified class initialization/structure).
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use gsm_params, only : zro, hlf, one, two, thr, four, fiv, &
         & six, seven, eight, nine, twpi
    use coalescenceClass, only: &
         & coalescenceParticle, &                    ! Particle type for coalescence sim.
         & coalescenceResults, newCoalescenceResults ! Results class and its constructor

    implicit none
    class(GSM),              intent(inout) :: gsmObj
    type(rxnSpecificModels), intent(inout) :: rxnPhysics
    type(GSMResults),        intent(inout) :: results
    class(OutputData),       intent(inout), optional :: outData

    integer(int32) :: k, lastProgenyIndx, particleID
    real(real64)   :: restMass, cosPhi, cosTheta, phi, totLinMome, &
         & sinPhi, sinTheta, kinEnergy, temp

! ======================================================================

    ! Class and its associated data
    type(coalescenceParticle), dimension(results%numProgeny) :: partBnk
    type(coalescenceResults) :: coalResults

! ======================================================================

    ! Copy fragment information to the coalescence data structures
    setupProgeny: do  k = 1,results%numProgeny

       ! Obtain interim results:
       kinEnergy  = results%progenyBnk(k)%kinEnergy     ! Kinetic energy
       restMass = results%progenyBnk(k)%restMass     ! Rest mass
       temp = kinEnergy*(kinEnergy + two*restMass)
       if (temp < 0.0d0) then
          temp = 0.01d0
          write(gsmObj%io%message,2100) "1149"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       totLinMome  = sqrt(temp)   ! total momentum

       ! A/Z of particle (and quantum decay number [strangeness])
       partBnk(k)%numBaryons= results%progenyBnk(k)%numBaryons
       partBnk(k)%charge = nint(results%progenyBnk(k)%numProtons)
       partBnk(k)%strangeness = zro
       ! Energy information
       partBnk(k)%kinEnergy  = kinEnergy    ! Kinetic energy [GeV]
       partBnk(k)%restMass   = restMass   ! Rest mass      [GeV/c**2]
       ! Momentum of particle
       partBnk(k)%linearMomX = totLinMome * &
           & results%progenyBnk(k)%sinTheta * cos(results%progenyBnk(k)%phi)
       partBnk(k)%linearMomY = totLinMome * &
           & results%progenyBnk(k)%sinTheta * sin(results%progenyBnk(k)%phi)
       partBnk(k)%linearMomZ = totLinMome * &
           & results%progenyBnk(k)%cosTheta
    end do setupProgeny


    ! Construct results object
    coalResults = newCoalescenceResults( partBnk, results%numProgeny )


    ! Simulate coalescence
    call rxnPhysics%coales%simulate( coalResults )


    if ( coalResults%simState > 0 ) then
       ! Check why simulation had an error:
       if ( coalResults%simState == 1 ) then
          ! Attempted simulation using more particles than the array actually holds
          ! NOTE: This error should ONLY occur in the event the client's particle bank
          !       has been exceeded, or if the client isn't properly tracking the number
          !       of produced particles.
       else if ( coalResults%simState == 2 ) then
          ! Attempted simulation using a negative amount of particles.
          ! NOTE: This should only occur if the client fails to properly track its own
          !       progeny information
       else if ( coalResults%simState >= 10 ) then
          ! =10 means that the Coalescence object wasn't constructed (requires data object)
          ! =11 means that the Results object wasn't constructed     (requires particle bank)
          stop
       end if
    end if


    ! Re-tally fragments (exit if nothing coalesced)
    if ( results%numProgeny == coalResults%numParticles  ) then
       return ! Nothing coalesced
    end if


    ! Update fragments:
    updateProgeny: do k = 1, coalResults%numParticles

       ! Obtain interim results:
       totLinMome  = sqrt(partBnk(k)%linearMomX**2 + partBnk(k)%linearMomY**2 + &
           & partBnk(k)%linearMomZ**2)
       restMass = partBnk(k)%restMass
       kinEnergy  = sqrt(totLinMome**2 + restMass**2) - restMass
       temp = totLinMome
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(gsmObj%io%message,2000) "1010"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       cosTheta  = partBnk(k)%linearMomZ/temp
       temp = one - cosTheta*cosTheta
       if (temp < 0.0d0) then
          temp = 0.01d0
          write(gsmObj%io%message,2100) "1016"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       sinTheta  = sqrt(temp)
       temp = totLinMome*sinTheta
       if (temp < div0Lim .and. temp > -div0Lim) then
          temp = div0Lim
          write(gsmObj%io%message,2000) "1022, 1023"
          call gsmObj%io%print(4, 3, gsmObj%io%message)
       end if
       cosPhi  = partBnk(k)%linearMomX/(temp)
       sinPhi  = partBnk(k)%linearMomY/(temp)
       phi  = atan2(sinPhi,cosPhi)
       if (phi < zro) phi = twpi + phi

       particleID = 0
       if (partBnk(k)%charge == 0 .and. partBnk(k)%numBaryons == 1) then
          particleID = one     
       elseif (partBnk(k)%charge == 1 .and. partBnk(k)%numBaryons == 1) then
          particleID = two    
       elseif (partBnk(k)%charge == 1 .and. partBnk(k)%numBaryons == 2) then
          particleID = thr   
          if (present(outData)) outData%ncoal(1) = outData%ncoal(1) + 1
        elseif (partBnk(k)%charge == 1 .and. partBnk(k)%numBaryons == 3)  then
          particleID = four
          if (present(outData)) outData%ncoal(2) = outData%ncoal(2) + 1
        elseif (partBnk(k)%charge == 2 .and. partBnk(k)%numBaryons == 3)  then
          particleID = fiv
          if (present(outData)) outData%ncoal(3) = outData%ncoal(3) + 1
        elseif (partBnk(k)%charge == 2 .and. partBnk(k)%numBaryons == 4)  then
          particleID = six
          if (present(outData)) outData%ncoal(4) = outData%ncoal(4) + 1
        elseif (partBnk(k)%charge == -1 .and. partBnk(k)%numBaryons == 0)  then
          particleID = seven
       elseif (partBnk(k)%charge ==  0 .and. partBnk(k)%numBaryons == 0)  then
          particleID = eight
       elseif (partBnk(k)%charge ==  1 .and. partBnk(k)%numBaryons == 0)  then
          particleID = nine
       else if (partBnk(k)%charge == 2 .and. partBnk(k)%numBaryons == 6) then
          particleID = 2004.d0   ! 6He
          if (present(outData)) outData%ncoal(5) = outData%ncoal(5) + 1
       else if (partBnk(k)%charge == 3 .and. partBnk(k)%numBaryons == 6) then
          particleID = 3003.d0   ! 6Li
          if (present(outData)) outData%ncoal(6) = outData%ncoal(6) + 1
       else if (partBnk(k)%charge == 3 .and. partBnk(k)%numBaryons == 7) then
          particleID = 3004.d0   ! 7Li
          if (present(outData)) outData%ncoal(7) = outData%ncoal(7) + 1
       else if (partBnk(k)%charge == 4 .and. partBnk(k)%numBaryons == 7) then
          particleID = 4003.d0   ! 7Be
          if (present(outData)) outData%ncoal(8) = outData%ncoal(8) + 1
       else
          ! LMK 12/2014
          write(gsmObj%io%message,1200)
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          write(gsmObj%io%message,1300) partBnk(k)%numBaryons, partBnk(k)%charge, kinEnergy
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          particleID = 1000 * partBnk(k)%charge + (partBnk(k)%numBaryons-partBnk(k)%charge)
       end if


       ! Store progeny information:
       results%progenyBnk(k)%numBaryons = partBnk(k)%numBaryons
       results%progenyBnk(k)%numProtons = partBnk(k)%charge
       results%progenyBnk(k)%kinEnergy  = kinEnergy
       results%progenyBnk(k)%restMass   = restMass
       results%progenyBnk(k)%phi        = phi
       results%progenyBnk(k)%theta      = atan2(sinTheta, cosTheta)
       results%progenyBnk(k)%sinTheta   = sinTheta
       results%progenyBnk(k)%cosTheta   = cosTheta
       results%progenyBnk(k)%typeID     = particleID
       if (partBnk(k)%numBaryons > 1)  then
          results%progenyBnk(k)%prodMech = 200.d0
       else
          results%progenyBnk(k)%prodMech = 0.0_real64
       endif

    end do updateProgeny

    ! Remove all coalesced INC secondaries (from coalesResults%numParticle up to results%numProgeny)
    lastProgenyIndx = results%numProgeny
    results%numProgeny = coalResults%numParticles
    removeCoalescedProgeny: do k = results%numProgeny+1,lastProgenyIndx
       results%progenyBnk(k) = GSMProgeny()
    end do removeCoalescedProgeny

    return
! ======================================================================
1200 format("Unidentified fragment after coalescence and INC physics.")
1300 format("   Unphysical product has properties (A = ", &
          & i3, ",  Z = ", i3, ", Ek = ", f8.3, ").")
2000 format("Divide by zero error prevented in ", &
          & "'coaleslInterface.f90', line ", A)
2100 format("Square root error prevented in ", &
          & "'coaleslInterface.f90', line ", A)
! ======================================================================
  end subroutine coalescenceInterface
