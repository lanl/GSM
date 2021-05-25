
  subroutine gemdec (evapObj, clientCompound, results )

! ======================================================================
!
!   ex is excitation energy PER NUCLEON initially.  AJS 10/29/03
!
!   Called by: EQDECY
!   Calls: STDCAY
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, October, 2003.
!    Changed by K.K. Gudima, using new fit by M.I. Baznat, Nov., 2004
!    bf12 added to common "/ifissc/" by KKG, December, 2004.
!    Edited by A. J. Sierk, LANL T-16, January, 2005.
!    Changed by K.K.G.  06/26/06
!    Editted by SGM 07/09/06 to include chages by KKG of 06/26/06
!       considering Fermi break-up in Preco and Evap when A < 13.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, two, thsn, amu
    use evaporationDataClass, only: newEvaporationData

    implicit none
    class(Evaporation),       intent(inout) :: evapObj
    type(evapCompound),       intent(in   ) :: clientCompound
    type(evaporationResults), intent(inout) :: results

    integer(int32) :: iaf, iar, iflag, ik, &
         & izf, izr, k, k2, kf, kk
    real(real64)   :: am, pmod, pnx, pny, pnz, temp, wtot
    logical        :: fisinh

    type(evapCompound) :: compound

    real(real64)               :: ekin = zro   ! Kinetic energy in the lab system
    real(real64), dimension(3) :: bet1 = zro   ! Unit vector of momentum of emitted particle
    type(evaporationCalculation) :: calcVars
                                                                                   
! ======================================================================

    ! Ensure evaporation object was constructed properly
    if ( .not. evapObj%constructed ) then
       write(evapObj%io%message, 3000)
       call evapObj%io%print(1, 1, evapObj%io%message)
       write(evapObj%io%message,3999)
       call evapObj%io%print(1, 1, evapObj%io%message)
       results%simState = 10
       return
    end if


    ! Ensure results were constructed:
    if ( .not. results%constructed ) then
       write(evapObj%io%message,3100)
       call evapObj%io%print(1, 1, evapObj%io%message)
       write(evapObj%io%message,3999)
       call evapObj%io%print(1, 1, evapObj%io%message)
       results%simState = 11
       return
    end if


    ! Create local copy of the compound (to be changed):
    compound = clientCompound


    ! Reset some important variables
    iflag = 0
    fisinh = .false.
    wtot = one


!   Start evaporation calculation:
!   iflag = 1 : No more emission
!   iflag = 2 : Fission
!   iflag = 3 : Emission occured
!     npref = 0
    do ik = 1, results%progenyBnkSize
       results%fissionProbability = zro

       ! Check for Fermi break-up calculation for small nuclei :
       if ( evapObj%fbuObj%recommendFermiBreakUp( compound%numBaryons ) )  then
          izf = nint(compound%numProtons)
          iaf = nint(compound%numBaryons)
          am = compound%numBaryons*amu + evapObj%energy (izf, iaf)
          temp = compound%recEnergy**2 + two*am*compound%recEnergy
          if (temp < 0.0d0) then
             temp = 0.01d0
             write(evapObj%io%message,1000) '100'
             call evapObj%io%print(4, 3, evapObj%io%message)
          end if
          pmod = sqrt(temp)/thsn  ! GeV/c
          pnx = pmod*compound%linearMomFrac(1)
          pny = pmod*compound%linearMomFrac(2)
          pnz = pmod*compound%linearMomFrac(3)

          ! Flag residual as having undergone Fermi Breakup (set values before Fermi Break-up occurred)
          results%compound%numBaryons   = compound%numBaryons
          results%compound%numProtons   = compound%numProtons
          results%compound%kinEnergy    = compound%kinEnergy
          results%compound%recEnergy    = 0.0_real64
          results%compound%linearMomTot = pmod
          results%compound%physicsFlag  = fermiBreakUpFlag

          ! Simulate Fermi Break-up physics:
          call evapObj%evapFermiInterface (compound%numBaryons, compound%numProtons, &
               & compound%kinEnergy, pnx, pny, pnz, results)

          return
       endif


       ! If bank nearly full, store residual instead
       if ( results%numProgeny == results%progenyBnkSizeM1 ) then
          ! Store residual as-is and exit
          write(evapObj%io%message,2100) compound%numBaryons, compound%numProtons, &
               & compound%kinEnergy, results%progenyBnkSize
          call evapObj%io%print(4, 3, evapObj%io%message)
          iflag = 1
       else
          ! Evaporate a random fragment
          call evapObj%stdcay (compound, iaf, izf, iflag, fisinh, ekin, bet1, &
               & results, calcVars)
          wtot = wtot*(one - results%fissionProbability)
       end if



       ! Store fragment into progeny array
       if (.not.fisinh) then
          ! No Fission (evaporated or done evaporating)
          if (iflag == 3) then

             ! Fragment was evaporated from compound nucleus; store fragment information
             results%numProgeny = results%numProgeny + 1
             results%progenyBnk(results%numProgeny)%numProtons = dble(izf)
             results%progenyBnk(results%numProgeny)%numBaryons = dble(iaf)
             results%progenyBnk(results%numProgeny)%kinEnergy  = ekin
             results%progenyBnk(results%numProgeny)%linearMomX = bet1(1)
             results%progenyBnk(results%numProgeny)%linearMomY = bet1(2)
             results%progenyBnk(results%numProgeny)%linearMomZ = bet1(3)
             results%progenyBnk(results%numProgeny)%restMass   = &
                  & evapObj%progenyMass( &
                  & results%progenyBnk(results%numProgeny)%numBaryons, &
                  & results%progenyBnk(results%numProgeny)%numProtons  &
                  & )

             results%progenyBnk(results%numProgeny)%physicsFlag = evaporationFlag


             ! Emit photons if client specified:
             if ( evapObj%usePhotoEmission ) then
                call evapObj%photonEmission( results%progenyBnk(results%numProgeny),  &
                     & compound%numBaryons, compound%numProtons, compound%kinEnergy )
             end if


          elseif (iflag == 1) then

             ! No more evaporation; store evaporated compound nucleus information
             results%compound%numBaryons  = compound%numBaryons
             results%compound%numProtons  = compound%numProtons
             results%compound%kinEnergy   = compound%kinEnergy
             results%compound%recEnergy   = compound%recEnergy
             !  Recoil momentum (nonrel):
             am = compound%numBaryons*amu + &
                  & evapObj%energy (nint(clientCompound%numProtons), nint(clientCompound%numBaryons))
             results%compound%linearMomTot = sqrt( &
                  & abs(compound%recEnergy**2 + two*am*compound%recEnergy) )
             results%compound%physicsFlag = evaporationFlag

             ! Probability of fission for evaporated compound nucleus
             results%fissionProbability = one - wtot


! AJS: Store evaporated compound nucleus in fragment bank
!   (03/04/04)
             results%numProgeny = results%numProgeny + 1
             results%progenyBnk(results%numProgeny)%numProtons = compound%numProtons
             results%progenyBnk(results%numProgeny)%numBaryons = compound%numBaryons
             results%progenyBnk(results%numProgeny)%kinEnergy  = compound%recEnergy
             results%progenyBnk(results%numProgeny)%linearMomX = compound%linearMomFrac(1)
             results%progenyBnk(results%numProgeny)%linearMomY = compound%linearMomFrac(2)
             results%progenyBnk(results%numProgeny)%linearMomZ = compound%linearMomFrac(3)
             results%progenyBnk(results%numProgeny)%restMass = &
                  & evapObj%progenyMass( &
                  & results%progenyBnk(results%numProgeny)%numBaryons, &
                  & results%progenyBnk(results%numProgeny)%numProtons  &
                  & )

             ! Evaporated nucleus; not from fission
             results%progenyBnk(results%numProgeny)%physicsFlag = evaporationFlag


             ! Emit photons if client specified:
             if ( evapObj%usePhotoEmission ) then
                call evapObj%photonEmission( results%progenyBnk(results%numProgeny),  &
                     & compound%numBaryons, compound%numProtons, compound%kinEnergy )
             end if

             return
          endif

       else
          ! Fission occurred, check for fission fragment evaporation/Fermi Breakup
          ! Evaporate fission fragments:
          do k = 1, results%numFissionFragments
             ! Obtain fission fragment (A,Z), energy, recoil energy, momentum and evaporate as if its a compound nucleus:
             compound%numBaryons   = results%fissionBnk(k)%numBaryons
             compound%numProtons   = results%fissionBnk(k)%numProtons
             compound%kinEnergy   = results%fissionBnk(k)%kinEnergy
             compound%recEnergy = results%fissionBnk(k)%rotEnergy
             do kk = 1,3
                compound%linearMomFrac(kk) = results%fissionBnk(k)%linearMomFrac(kk)
             end do

             ! Evaporate fission fragment
             do k2 = results%numProgeny, results%progenyBnkSize

                ! Check if small enough to undergo Fermi Breakup process
                if ( evapObj%fbuObj%recommendFermiBreakUp( compound%numBaryons ) )  then
                   izf = nint(compound%numProtons)
                   iaf = nint(compound%numBaryons)
                   am = compound%numBaryons*amu + evapObj%energy (izf, iaf)
                   temp = compound%recEnergy**2 + two*am*compound%recEnergy
                   if (temp < 0.0d0) then
                      temp = 0.01d0
                      write(evapObj%io%message,1000) '100'
                      call evapObj%io%print(4, 3, evapObj%io%message)
                   end if
                   pmod = sqrt(temp)/thsn  ! GeV/c
                   pnx = pmod*compound%linearMomFrac(1)
                   pny = pmod*compound%linearMomFrac(2)
                   pnz = pmod*compound%linearMomFrac(3)
                   call evapObj%evapFermiInterface (compound%numBaryons, compound%numProtons, &
                        & compound%kinEnergy, pnx, pny, pnz, results)

                   ! Flag the fission fragment as having undergone Fermi Breakup
                   results%fissionBnk(k)%physicsFlag = fermiBreakUpFlag

                   ! Move on to next fission fragment (or exit)
                   go to 100

                endif


                ! Ensure bank limit is not exceeded
                if ( results%numProgeny >= results%progenyBnkSizeM1) then
                   ! Bank limit at maximum; store last residual and return
                   write(evapObj%io%message,2100) clientCompound%numBaryons, clientCompound%numProtons, &
                        & clientCompound%kinEnergy, results%progenyBnkSize
                   call evapObj%io%print(3, 3, evapObj%io%message)
                   write(evapObj%io%message, 2110)
                   call evapObj%io%print(3, 3, evapObj%io%message)
                   write(evapObj%io%message, 2999)
                   call evapObj%io%print(3, 3, evapObj%io%message)
                   iflag = 1
                else
                   ! Check for evaporation or fission progeny
                   call evapObj%stdcay (compound, iaf, izf, iflag, fisinh, &
                        & ekin, bet1, results, calcVars)
                end if


                ! If the fragment evaporated progeny, or if it's done evaporating, store [else fission (ignore)]
                ! BUG: The fission fragment information IS updated if one of the original fission fragments fissions (see fisgem routine);
                !       thus providing the user incorrect information regarding fission fragments (if only allowing one).
                ! SUGGESTION: Move this information to a subroutine and RECURSIVELY call it (i.e. if fission fragments fission, then call the routine)
                !             Additionally, create a LIST of fission fragments (user to provide a pointer to the class for the array)
                if (iflag == 3) then

                   ! Evaporation progeny from fission fragment
                   results%numProgeny = results%numProgeny + 1
                   results%progenyBnk(results%numProgeny)%numProtons = dble(izf)
                   results%progenyBnk(results%numProgeny)%numBaryons = dble(iaf)
                   results%progenyBnk(results%numProgeny)%kinEnergy  = ekin
                   results%progenyBnk(results%numProgeny)%linearMomX = bet1(1)
                   results%progenyBnk(results%numProgeny)%linearMomY = bet1(2)
                   results%progenyBnk(results%numProgeny)%linearMomZ = bet1(3)
                   results%progenyBnk(results%numProgeny)%restMass = &
                        & evapObj%progenyMass( &
                        & results%progenyBnk(results%numProgeny)%numBaryons, &
                        & results%progenyBnk(results%numProgeny)%numProtons  &
                        & )
                   ! Fission fragment residue:
                   results%progenyBnk(results%numProgeny)%physicsFlag = fissionFlag


                   ! Emit photons if client specified:
                   if ( evapObj%usePhotoEmission ) then
                      call evapObj%photonEmission( results%progenyBnk(results%numProgeny),  &
                           & compound%numBaryons, compound%numProtons, compound%kinEnergy )
                   end if


                else if ( iflag == 2 ) then
                   ! Flag fission fragment as supposed to undergo fission again; go on to next fragment
                   results%fissionBnk(k)%physicsFlag = fissionFlag
                   go to 100
                elseif (iflag == 1) then

                   ! No more evaporation; store fission fragment into progeny array
                   izr = nint(compound%numProtons)
                   iar = nint(compound%numBaryons)
                   results%numProgeny = results%numProgeny + 1
                   results%progenyBnk(results%numProgeny)%numProtons = dble(izr)
                   results%progenyBnk(results%numProgeny)%numBaryons = dble(iar)
                   results%progenyBnk(results%numProgeny)%kinEnergy  = compound%recEnergy
                   results%progenyBnk(results%numProgeny)%linearMomX = compound%linearMomFrac(1)
                   results%progenyBnk(results%numProgeny)%linearMomY = compound%linearMomFrac(2)
                   results%progenyBnk(results%numProgeny)%linearMomZ = compound%linearMomFrac(3)
                   results%progenyBnk(results%numProgeny)%restMass   = &
                        & evapObj%progenyMass( &
                        & results%progenyBnk(results%numProgeny)%numBaryons, &
                        & results%progenyBnk(results%numProgeny)%numProtons  )
                   ! Fission fragment residue:
                   results%progenyBnk(results%numProgeny)%physicsFlag = fissionFlag


                   ! Update the state of the fission fragment (A,Z)
                   results%fissionBnk(k)%numBaryons = compound%numBaryons
                   results%fissionBnk(k)%numProtons = compound%numProtons
                   results%fissionBnk(k)%kinEnergy  = compound%kinEnergy
                   results%fissionBnk(k)%rotEnergy  = compound%recEnergy
                   do kf = 1, 3
                      results%fissionBnk(k)%linearMomFrac(kf) = compound%linearMomFrac(kf)
                   end do


                   ! Emit photons if client specified:
                   if ( evapObj%usePhotoEmission ) then
                      call evapObj%photonEmission( results%progenyBnk(results%numProgeny),  &
                           & compound%numBaryons, compound%numProtons, compound%kinEnergy )
                   end if

                   ! Evaporate next fission fragment (or exit)
                   go to 100

                endif
             end do
             write(evapObj%io%message, 2000) results%progenyBnkSize
             call evapObj%io%print(3, 3, evapObj%io%message)
             write(evapObj%io%message, 2999) results%progenyBnkSize
             call evapObj%io%print(3, 3, evapObj%io%message)
100          continue
          end do

          ! Update fission probability
          results%fissionProbability = one - wtot

          return
       endif
    end do
    write(evapObj%io%message, 2000) results%progenyBnkSize
    call evapObj%io%print(3, 3, evapObj%io%message)
    write(evapObj%io%message, 2999) results%progenyBnkSize
    call evapObj%io%print(3, 3, evapObj%io%message)

    return
! ======================================================================
1000 format("Square root error prevented in 'gemdec.f90', line ", A)
2000 format("Evaporation progeny exceed the received bank size.", &
          & " Please increase the bank size (currently ", i5, ").")
2100 format("Evaporation of residual (A=", f5.1, ", Z=", f5.1, &
          & ", U=", f5.1, ") has created more progeny (", i5, ") than allowed ", &
          & "by the progeny bank size.")
2110 format("   Residual will be stored and evaporation/fission will ", &
          & "be assumed to cease.")
2999 format("   Results may be suspect.")
3000 format("The evaporation object was not constructed ", &
          & "prior to simulation.")
3100 format("The results object was not constructed for the ", &
          & "evaporation and fission processes.")
3999 format("   The evaporation and fission processes ", &
          & "will not be simulated.")
! ======================================================================
  end subroutine gemdec


! ======================================================================

  function progenyMass(evapObj, af, zf) result(rstMass)

! ======================================================================
!
! Function returns rest mass of fragment using a 'simple' formula
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use evaporationParams, only: emnuct, emelec, thsn

    implicit none
    class(Evaporation), intent(in   ) :: evapObj
    real(real64),       intent(in   ) :: af
    real(real64),       intent(in   ) :: zf
    real(real64)                      :: rstMass

    integer(int32) :: izf, inf
    real(real64)   :: emxp

! ======================================================================

    izf = nint(zf)
    inf = nint(af) - izf

    if (izf.ne.0) then
       emxp = evapObj%evapMolnix%defineEnergy(izf, inf, 2)
    else
       emxp = 8.071d0
    endif

    rstMass = emnuct*af + emxp/thsn - zf*emelec

    return 
! ======================================================================
  end function progenyMass
