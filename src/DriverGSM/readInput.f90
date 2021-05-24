
  subroutine readInput( inputFile )

! ============================================================================
!
! Reads the input file
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64, int64
    use gsm_params, only: parName, numParticleNames
    use generalizedSpallationData, only: initializeGSMData
    use gsmClass, only: &
         & GSM, newGSM, &
         & GSMProjectile, GSMTarget,  &
         & gsmOptions, &
         & gsmOutput, &
         & stopName
    use defaults_mod, only : verbose

    implicit none
    character(LEN=*), intent(in   ) :: inputFile

    integer(int32) :: iter = 0_int32

    integer(int32) :: ityp, &
         & nevtype, &
         & npreqtyp, aproj0, zproj0, ms0, ih, indx, j
    real(real64)   :: atarg0, dang, dt0, t0max, tempr, ztarg0, &
         & t0mev

      ! For reading ZAID and particle verification
    integer(int32) :: iostat, zaid

    character(len=6) :: pnametemp
    character(len=4) :: pname

    logical, save :: first = .TRUE.

    ! Deprecated
    integer(int64) :: limc
    integer(int32) :: mang, mdubl, misy, mchy, mpyld, mspec
    integer(int32) :: ipisa, ihist, ipar1, ipar2
    real(real64) :: dteta
    real(real64), dimension(10) :: tet1, tet2, ang
    real(real64), dimension( 4) :: tmin, tmax, dt


! ======================================================================

    ! GSM object structure:
    type(GSM), allocatable :: gsmObj
    type(gsmOptions)       :: options

    ! For a single simulation
    type(GSMProjectile) :: proj
    type(GSMTarget)     :: targ

    ! For the output file specification
    type(GSMOutput) :: output

! ======================================================================

    logical, parameter :: &
         & read_input = .TRUE.
    character(len=128) :: dummy

! ======================================================================

    ! Allocate data for GSM model physics:
    if (.not.allocated(gsmObj)) allocate(gsmObj)

5   continue

    if(read_input) then
! Initial files and tables
       if(first) then
          first = .false.

          output%inputFile = inputFile
          open (15, file=inputFile, status='old', form="formatted", action="read")
          if ( verbose >= 2 ) write(*,1225) 'reading input file "', trim(adjustl(inputFile)), '".'
          read (15, 1200) output%auxiliaryFile
          read (15, 1200) output%outputFile

          ! Remove all trailing text when 2 or more contiguous spaces are found
          do indx = 1, len_trim(output%outputFile)-2
             if (len_trim(output%outputFile(indx:indx+2)) == 0) &
                 & output%outputFile = trim(adjustl(output%outputFile(:indx)))
          end do
          do indx = 1, len_trim(output%auxiliaryFile)-2
             if (len_trim(output%auxiliaryFile(indx:indx+2)) == 0) &
                 & output%auxiliaryFile = trim(adjustl(output%auxiliaryFile(:indx)))
          end do
          open (31, file=output%outputFile, status='unknown')
          open (16, file=output%auxiliaryFile, status='unknown')
          write (16, 1200) trim(output%outputFile)
          write (16, 1200) trim(output%auxiliaryFile)

       endif

         ! Primary read input
       read (15, *) pnametemp
       ms0 = 0

         ! Checking For Standard Particle Name, skips ZAID determination if a match is found
       do indx = 1, numParticleNames
          if ( trim(pnametemp) == parname(indx) ) then
             aproj0 = 0
             zproj0 = 0
             pname = parname(indx)
             go to 10
          endif
       enddo

         ! Converting String to Real Number if common parname not found
       read(pnametemp, *, iostat = iostat) tempr
       if ( iostat.ne.0 ) then
          write(*,1225) "error reading zaid from ", trim(pnametemp), "."
          return
       endif

         ! Getting A, Z Values from ZAID
       zaid = nint(tempr)

         ! Obtaining A, Z, ms0 based on ZAID
       if ( zaid <= 1000 .and. zaid >= -1 ) then

          ! Go back one line
          backspace(15)

            ! Using "A, Z, ms0" card (A = -1 for pions)
          read (15, *) aproj0, zproj0, ms0
          pname = parname(8)

       elseif ( zaid > 1000 ) then

            ! Using "ZAID" card (ms0 = 0 assumed)
          aproj0 = mod(zaid, 1000)
          zproj0 = floor(zaid/1000.0)
          ms0 = 0
          pname = parname(8)

       else

          write (*, 1211) 'unable to reading incident particle information.  please correct'
          write (*, 1211) 'to the standard "zaid" or "a, z, (stin)" format.'
          write (*, 1211) 'program exiting...'
          return

       endif

10     continue

       ! Check for program end
       if (pname == parname(numParticleNames)) then
          ! return
          proj%particleName = stopName
          call gsmObj%generateOutput(proj, targ, output) ! calls to main software library, ends program there.
          return
       endif

       ! Reset the ouput object
       output = GSMOutput()

       read (15, 1200) dummy   ! continuing list of pname options
       read (15, 1200) dummy   ! continuing list of pname options
       read (15, *) t0mev
       read (15, *) atarg0
       read (15, *) ztarg0
       read (15, *) output%numInelasticEvents
       read (15, *) dt0
       read (15, *) t0max
       read (15, *) output%deltaTheta
       read (15, *) output%energySpectra
       read (15, *) output%multiplicities
       read (15, *) output%channelCrossSection
       read (15, *) output%nuclideCrossSection
       read (15, *) output%doubleDiffSpectra
       read (15, *) output%angularSpectra
       read (15, *) output%minEjectileRange, output%maxEjectileRange
       output%spectraEjectiles(:) = .FALSE.
       do iter = output%minEjectileRange, output%maxejectileRange
          output%spectraEjectiles(iter) = .TRUE.
       end do
       read (15, *) (output%angleBins(j)%lowerBound, &
            & output%angleBins(j)%upperBound, j = 1,10)
       read (15, *) (output%energyBins(j)%lowerBound, output%energyBins(j)%upperBound, &
            & output%energyBinSubStep(j), j = 1,4)
       read (15, *) nevtype
       read (15, *) npreqtyp
       read (15, *) ipisa
       if(ipisa > 0) output%printPISA = .TRUE.
       if (output%printPISA) read(15,*) output%pisaAngles, output%pisaDTheta
       read (15, *) ihist
       if ( ihist > 0 ) output%printHist = .TRUE.
       read (15, *) ityp
!  Read in up to ten lines of text to write at the top of
!  the output/auxiliary files and on the terminal as comments on the run.
       read (15, *) output%numComments
       if (output%numComments > 0) then
          read (15, 1200) (output%comments(ih), ih=1,output%numComments)
          if ( verbose >= 0 ) then
             write ( *, 1900)
             write ( *, 2100)
             do ih = 1,output%numComments
                write ( *, 1205) output%comments(ih)
             end do
             write ( *, 2100)
             write ( *, 1900)
          endif
       endif
    endif
    if ( verbose >= 2 ) write (*, 1250) trim(inputFile)


    ! Ensure energies are all positive:
    if ( t0mev < 0 ) then
       ! Has units of GeV/A:
       if ( aproj0 > 0 ) then
          t0mev = abs(t0mev) * dble(aproj0) * 1000.0_real64
       else
          t0mev = abs(t0mev) * 1000.0_real64
       end if

       ! Check if dt0 > 0 [i.e. not used]
       if ( dt0 > 0 ) then
          dt0 = 0
       end if
    else
       ! Positive energy [units=MeV]
       if ( dt0 < 0 ) then
          dt0 = 0
       end if
    end if


    ! Construct the GSM object and its data:
    call initializeGSMData(ityp)

    ! Update physics options:
    options = gsmOptions()   ! Reset the object
    options%preeqOpts%numPreeqType = npreqtyp   ! Apply input specified values:
    options%evapOpts%numEvapType   = nevtype    ! Apply input specified values:
    gsmObj = newGSM( clientOptions = options )


    ! Setup projectile and target objects:
    ! (projectile)
    proj%particleName = pname
    proj%numBaryons = aproj0
    proj%numProtons = zproj0
    proj%decayNumber = ms0
    proj%kinEnergy  = t0mev / 1000.0_real64
    if ( dt0 > 0 ) then
       proj%dKinEnergy = dt0
    end if
    proj%kinEnergyMax = t0max / 1000.0_real64
    ! (target)
    write(targ%particleName, 3000) nint(1000.0*ztarg0) + nint(atarg0)
    targ%numBaryons = atarg0
    targ%numProtons = ztarg0

    ! Simulate physics now:
    call gsmObj%generateOutput(proj, targ, output)


    ! go to top and re-read input
    go to 5

! ======================================================================
1200 format (a)
1205 format (1x, a)
1211 format (5x, 'Error: ', a)
1225 format (3x, "Comment: ", a, a, a)
1250 format (3x, 'Comment: successfully read input file "', a, '".')
1900 format (1x)
2100 format (1x, "-----------------------------------")
3000 format (i6)
! ======================================================================
  end subroutine readInput
