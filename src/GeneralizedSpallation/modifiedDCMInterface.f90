
  subroutine modifiedDCMInterface( gsmObj, gsmRxn )

! ======================================================================
!
! Interface for Modified DCM (INC) physics - performs simulation of mDCM
! and obtains results for GSM to use.
!    NOTE: This model is NOT parallelized at this time!
!
!
! Written by CMJ, XCP-3 (04/2019)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, int64, real64
    use gsm_params, only: hlf
    use modifiedDCMData, only: jamin, jamax
    use modifiedDCMClass, only: &
         & mDCMProgeny, &
         & results, newMDCMResults
    
    implicit none
    class(GSM),         intent(inout) :: gsmObj
    class(GSMReaction), intent(inout) :: gsmRxn

    integer(int32) :: iret, ja1, ja2, jz1, jz2, numElastic, mDCMProgIndx, &
         & xyz, numGSMProg, particleID

! ======================================================================

    ! Construct the global results object within the mDCM class module:
    type(mDCMProgeny), dimension( gsmRxn%results%maxProgeny ) :: progenyBnk

! ======================================================================

    ! NOTE: EVERYTHING ENCLOSED BY THIS 'IF' STATEMENT SHOULD BE MOVED
    !       INTO THE MDCM (a mDCM interface for clients to use - check
    !       proper object construction/setup, simulate, check results)!!!
    logical, parameter :: checkErrors = .TRUE.

    ! ----------------------------------------
    ! Common blocks to interface to the MDCM:
    ! ----------------------------------------
    ! FOR RESULTS FROM MDCM:
    real(real64)   :: an1,an2,zn1,zn2,enext1,enext2,pnucl1, pnucl2, amnuc1, amnuc2
    common/resultlaq/an1,an2,zn1,zn2,enext1,enext2,pnucl1(3),&
         & pnucl2(3), amnuc1(3),amnuc2(3)

    ! OPTIONS FOR THE SIMULATION:
    real(real64)   :: rn, delta
    common/rndelt/ rn,delta

! ======================================================================

    ! Allow only one thread at a time here
    !$OMP critical

    ! INC for LAQGSM type calculations
    ! NOTE: This has access to the 'results' object through the module!
35  continue
    results = newMDCMResults( progenyBnk )
    call cascaw(numElastic,rn,delta,iret)

    if ( checkErrors ) then
       ! Error occurred in cascade, re-running.
       if(iret == 1 ) then
          write(gsmObj%io%message, 3200) "An unknown Modified DCM error occurred", gsmRxn%eventNum
          call gsmObj%io%print(5, 3, gsmObj%io%message)
          go  to  35
       endif

       ! sgm, 12/04/05: A, Z values for "new" particles
       jz1 = nint(zn1)   ! z for projectile
       jz2 = nint(zn2)   ! z for target
       ja1 = nint(an1)   ! a for projectile
       ja2 = nint(an2)   ! a for target

       ! Checking validity of residual nuclei
       ! Neutron star in either target or projectile
       if( (ja1 >= 12 .and. jz1 == 0) .or. (ja2 >= 12 .and. jz2 == 0) ) then
          write(gsmObj%io%message,3200) "A large neutron cluster (A>11) remains", gsmRxn%eventNum
          call gsmObj%io%print(5, 3, gsmObj%io%message)
          go to 35
       endif

       ! NO Charge of projectile or target depleted (this is allowed with LAQGSM, NOT an error really)
       if(jz1 == 0 .or. jz2 == 0) then
          ! Prints when a residual does NOT exist (handled later)
          ! write(gsmObj%io%message,3200) "There is no charge in the residual nucleus", gsmRxn%eventNum
          ! call gsmObj%io%print(6, 4, gsmObj%io%message)
          go to 37
       endif

       ! Projectile out of specified j- A value bounds
       if( (ja1 >= 12 .and. ja1 < jamin(jz1) ) .or. &
            & (ja1 >= 12 .and. ja1 > jamax(jz1) ) ) then
          write(gsmObj%io%message,3200) "Projectile residual has A>11; out of bounds", gsmRxn%eventNum
          call gsmObj%io%print(5, 3, gsmObj%io%message)
          go to 35
       endif

       ! Target out of specified j- A value bounds
       if( (ja2 >= 12 .and. ja2 < jamin(jz2) ) .or. &
            & (ja2 >= 12 .and. ja2 > jamax(jz2) ) ) then
          write(gsmObj%io%message,3200) "Target residual has A>11; out of bounds", gsmRxn%eventNum
          call gsmObj%io%print(5, 3, gsmObj%io%message)
          go to 35
       endif
    end if

    ! Valid Residual Nuclei
37  continue


    ! Obtain progeny:
    progenyLoop: do mDCMProgIndx = 1, results%numProgeny

       ! Ensure arrays aren't exceeded:
       if ( gsmRxn%results%numProgeny > gsmRxn%results%maxProgenyM1 ) then
          write(gsmObj%io%message, 1000)
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          write(gsmObj%io%message, 1010) results%numProgeny - mDCMProgIndx + 1
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          exit progenyLoop
       end if

       ! Obtain interim results:
       particleID = 0
       if ( progenyBnk(mDCMProgIndx)%restMass < hlf ) then
          if ( progenyBnk(mDCMProgIndx)%numProtons < 0 ) then
             particleID = 7   ! Pi-
          else if ( progenyBnk(mDCMProgIndx)%numProtons == 0 ) then
             particleID = 8   ! Pi0
          else
             particleID = 9   ! Pi+
          end if
       else
          if ( progenyBnk(mDCMProgIndx)%numProtons <= 0 ) then
             particleID = 1   ! Neutron
          else
             particleID = 2   ! Proton
          end if
       end if
       if ( particleID == 0 ) then
          write(gsmObj%io%message, 1100)
          call gsmObj%io%print(3, 3, gsmObj%io%message)
          particleID = 1   ! Assume a neutron
       end if

       ! Tally progeny:
       gsmRxn%results%numProgeny = gsmRxn%results%numProgeny + 1
       numGSMProg = gsmRxn%results%numProgeny
       gsmRxn%results%progenyBnk( numGSMProg )%numBaryons = progenyBnk(mDCMProgIndx)%numBaryons
       gsmRxn%results%progenyBnk( numGSMProg )%numProtons = progenyBnk(mDCMProgIndx)%numProtons
       gsmRxn%results%progenyBnk( numGSMProg )%kinEnergy  = progenyBnk(mDCMProgIndx)%kinEnergy
       gsmRxn%results%progenyBnk( numGSMProg )%restMass   = progenyBnk(mDCMProgIndx)%restMass
       gsmRxn%results%progenyBnk( numGSMProg )%phi        = progenyBnk(mDCMProgIndx)%phi
       gsmRxn%results%progenyBnk( numGSMProg )%theta      = progenyBnk(mDCMProgIndx)%theta
       gsmRxn%results%progenyBnk( numGSMProg )%sinTheta   = sin( progenyBnk(mDCMProgIndx)%theta )
       gsmRxn%results%progenyBnk( numGSMProg )%cosTheta   = cos( progenyBnk(mDCMProgIndx)%theta )
       gsmRxn%results%progenyBnk( numGSMProg )%typeID     = particleID
       gsmRxn%results%progenyBnk( numGSMProg )%prodMech   = 0
    end do progenyLoop


    ! Setup target residual:
    gsmRxn%results%targRes%kinEnergy  = enext2   ! in [GeV] currently!
    gsmRxn%results%targRes%numBaryons = nint(an2)
    gsmRxn%results%targRes%numProtons = nint(zn2)
    do xyz = 1, 3
       gsmRxn%results%targRes%linearMom(xyz)  = pnucl2(xyz)
       gsmRxn%results%targRes%angularMom(xyz) = amnuc2(xyz)
    enddo
    gsmRxn%results%targExc%numTotal    = results%targExc%numTotal
    gsmRxn%results%targExc%numNeutrons = results%targExc%numTotal - &
         & results%targExc%numHoles - results%targExc%numProtons
    gsmRxn%results%targExc%numProtons  = results%targExc%numProtons
    gsmRxn%results%targExc%numHoles    = results%targExc%numHoles

    ! Setup projectile residual:
    gsmRxn%results%projRes%kinEnergy  = enext1   ! in [GeV] currently!
    gsmRxn%results%projRes%numBaryons = nint(an1)
    gsmRxn%results%projRes%numProtons = nint(zn1)
    do xyz = 1, 3
       gsmRxn%results%projRes%linearMom(xyz)  = pnucl1(xyz)
       gsmRxn%results%projRes%angularMom(xyz) = amnuc1(xyz)
    end do
    gsmRxn%results%projExc%numTotal    = results%projExc%numTotal
    gsmRxn%results%projExc%numNeutrons = results%projExc%numTotal - &
         & results%projExc%numHoles - results%projExc%numProtons
    gsmRxn%results%projExc%numProtons  = results%projExc%numProtons
    gsmRxn%results%projExc%numHoles    = results%projExc%numHoles

    ! Setting total number of elastic events and resetting number
    gsmRxn%results%numElasticEvents = gsmRxn%results%numElasticEvents + numElastic
    gsmRxn%outData%intel = gsmRxn%outData%intel + numElastic

    !$OMP end critical

    return
! ======================================================================
1000 format("The GSM progeny bank was filled during the mDCM simulation.")
1010 format("   The last ", i3, " progeny cannot be tallied.")
1100 format("An unknown particle was found in the mDCM-GSM interface.")
3200 format (A, " (event ", i8, ").  Event will be restarted...")
! ======================================================================
  end subroutine modifiedDCMInterface
