
  function initializeEvaporationData( &
       & massFile, &
       & levelFile, &
       & shellFile) &
       & result(errorFlag)

! ======================================================================
!
!  SETUP
!      SETUP Nuclear Data and some parameters
! =====================================================================
!
! <Files> 
!  mass.tbl       :  Mass Excess Table    
!  level.tbl      :  Nuclear levels of the ejectiles
!  shell.tbl      :  Shell effect data for Myers-Swiatecki fission
!                    barrier calculation
!
!    "Last" change: 13-AUG-2003 by NVMokhov
!    Edited by A. J. Sierk, LANL T-16, September, 2003.
!    Edited by AJS, LANL T-2, December, 2011.
!
! ======================================================================

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use evaporationParams, only: zro, one, ato3rd

    implicit none
    character(LEN=*), intent(in   ), optional :: massFile
    character(LEN=*), intent(in   ), optional :: levelFile
    character(LEN=*), intent(in   ), optional :: shellFile
    integer(int32)                            :: errorFlag

    ! Variables used for the files (names, unit numbers)
    integer(int32) :: massFUnit, levelFUnit, shellFUnit
    character(LEN=128) :: massFileName = defaultMassFile
    character(LEN=128) :: levelFileName = defaultLevelFile
    character(LEN=128) :: shellFileName = defaultShellFile
    logical            :: fileExists

    integer(int32) :: i, ia, iaa0, iz, j, jmax, k
    real(real64)   :: sh

! ======================================================================

    integer(int32) :: &
         & massFailed  = noErrorsFound, &
         & shellFailed = noErrorsFound, &
         & levelFailed = noErrorsFound

! ======================================================================

    ! Assume all data will be initialized without errors
    errorFlag = noErrorsFound


    ! ----------------------------------------------------------------
    ! MASS DATA FILE
    ! ----------------------------------------------------------------
    if ( present(massFile) ) massFileName = massFile
5   continue

    ! Determine if file exists
    inquire(file = massFileName, exist = fileExists)
    if ( fileExists ) then
       open (newunit = massFUnit, &
            & file = massFileName, &
            & status="old", action="read")
    else
       ! File not in directory
       if ( massFileName /= defaultMassFile ) then
          ! Look for default file instead
          write(*,2100) trim(massFileName), trim(defaultMassFile)
          massFileName = defaultMassFile
          go to 5
       else
          write(*,2150) trim(massFileName)
          massFailed = massDataFailed
       end if
    end if

    ! ----------------------------------------------------------------
    ! Read mass excess table:
    ! ----------------------------------------------------------------
    if (.not.allocated(wapsm)) allocate(wapsm(0:150, 0:250))
    wapsm(:,:) = zro
    if ( massFailed /= massDataFailed ) then
       do j = 1,10000
          read (massFUnit, *, end = 10) iz, ia, wapsm(iz,ia-iz)
          if (ia-iz > 250 .or. iz > inn) then
             write(*,1000) iz, (ia-iz)
             massFailed = massDataFailed
          endif
       end do
       write (*,1100) iz, ia, wapsm(ia,ia-iz)
10     close (massFUnit)
    end if
    if ( massFailed /= noErrorsFound ) errorFlag = errorFlag + massDataFailed



    ! ----------------------------------------------------------------
    ! LEVEL DATA FILE
    ! ----------------------------------------------------------------
    if ( present(levelFile) ) levelFileName = levelFile
6   continue

    ! Determine if file exists
    inquire(file = levelFileName, exist = fileExists)
    if ( fileExists ) then
       open (newunit = levelFUnit, &
            & file = levelFileName, &
            & status="old", action="read")
    else
       ! File not in directory
       if ( levelFileName /= defaultLevelFile ) then
          ! Look for default file instead
          write(*,2100) trim(levelFileName), trim(defaultLevelFile)
          levelFileName = defaultLevelFile
          go to 6
       else
          write(*,2150) trim(levelFileName)
          levelFailed = levelDataFailed
       end if
    end if

    ! ----------------------------------------------------------------
    ! Read data for the excited states of the ejectiles:
    ! ----------------------------------------------------------------
    if (.not.allocated(exm))   allocate(exm(maxSize, 200))
    if (.not.allocated(spin))  allocate(spin(maxSize, 200))
    if (.not.allocated(width)) allocate(width(maxSize, 200))
    exm(:, :) = zro
    spin(:, :) = zro
    width(:, :) = zro

    if ( levelFailed /= levelDataFailed ) then
       do i = 1,maxSize
          read (levelFUnit, *, end=20, err=20) jmax
          if (jmax > 0) then
             do j = 1,jmax
                read (levelFUnit, *) width(i,j), spin(i,j), exm(i,j)
             end do
          endif
       end do
20     close (levelFUnit)
    end if
    if ( levelFailed /= noErrorsFound ) errorFlag = errorFlag + levelDataFailed



    ! ----------------------------------------------------------------
    ! SHELL DATA FILE
    ! ----------------------------------------------------------------
    if ( present(shellFile) ) shellFileName = shellFile
7   continue

    ! Determine if file exists
    inquire(file = shellFileName, exist = fileExists)
    if ( fileExists ) then
       open (newunit = shellFUnit, &
            & file = shellFileName, &
            & status="old", action="read" )
    else
       ! File not in directory
       if ( shellFileName /= defaultShellFile ) then
          ! Look for default file instead
          write(*,2100) trim(shellFileName), trim(defaultShellFile)
          shellFileName = defaultShellFile
          go to 7
       else
          write(*,2150) trim(shellFileName)
          shellFailed = shellDataFailed
       end if
    end if

    ! ----------------------------------------------------------------
    ! Read shell effect data for fission-barrier calculations:
    ! ----------------------------------------------------------------
    if (.not.allocated(shellc)) allocate(shellc(inn, 250))
    shellc(:, :) = zro

    if ( shellFailed /= shellDataFailed ) then
       do i = 1,10000
          read (shellFUnit, *, end=30) j, k, sh
          if (k-j < 1) then
             write (*, 1200) j, k
             shellFailed = shellDataFailed
          else
             shellc (j,k-j) = sh
          endif
       end do
30     close (shellFUnit)
    end if
    if ( shellFailed /= noErrorsFound ) errorFlag = errorFlag + shellDataFailed


    ! ----------------------------------------------------------------
    ! Establish parameterized data now
    ! ----------------------------------------------------------------

    !  Set Dostrovsky's parameter set, the footnote of
    !     PR116(1959)683 (p.699)
    !  . kp
    t(1,1) = 0.51d0
    t(1,2) = 0.60d0
    t(1,3) = 0.66d0
    t(1,4) = 0.68d0
    !  . k_alpha
    t(2,1) = 0.81d0
    t(2,2) = 0.85d0
    t(2,3) = 0.89d0
    t(2,4) = 0.93d0
    !  . cp
    t(3,1) =  0.0d0 
    t(3,2) = -0.06d0
    t(3,3) = -0.10d0
    t(3,4) = -0.10d0


    ! Setup 'paire0' and 'st0' arrays
    if (.not.allocated(paire0)) allocate(paire0(iiz, inn))
    if (.not.allocated(st0))    allocate(st0(iiz, inn))
    do i = 1, iiz
       do j = 1, inn
          paire0(i, j) = pz(i) + pn(j)
          st0(i, j)    = sz(i) + sn(j)
       end do
    end do


    !   alpha and beta for neutrons; the precise parameter set:
    do iaa0 = 1,aMax
       ! For alpp
       alpp(iaa0) = 0.76d0 + 1.93d0/ato3rd(iaa0)

       ! For betap
       betap(iaa0) = (1.66d0/ato3rd(iaa0)**2 - 5.d-2)/alpp(iaa0)
       betap(iaa0) = max(betap(iaa0), zro)
    end do


    if ( errorFlag == noErrorsFound ) then
       evaporationDataEstablished = .TRUE.
       errorFlag = evaporationDataInitialized
    end if

    return

! ======================================================================
1000 format(5x, "Fatal error: Mass table reading error in subroutine", &
          & "'inigem'.", /, 5x, "Fatal error:", 3x, &
          & "Z (", i4, ") or A-Z (", i4, ") in ", &
          & "mass. tbl is larger than 150 or 250")
1100 format(3x, "Warning: The last record read in subroutine ", &
          & "'inigem' is for A = ", i4, ", Z = ", i4, ", with a ", &
          & "value of ", f8.4, ".", /, 3x, &
          "Warning:   mass.tbl should have less than 10000 records.")
1200 format(5x, "Fatal error: Shell table reading error;  iz = ", i4, &
          & ", ia = ", i4)
2100 format(3x, "Warning: The file '", A, "' does not exist in the ", &
          & "program's directory.", /, 3x, "Warning:    Will look for ", &
          & "the default data file, '", A, "', to use in the Evaporation ", &
          & "simulation.")
2150 format(5x, "Error: The file '", A, "' does not exist in the ", &
          & "program's directory.", /, 5x, "Error:    Cannot obtain ", &
          & "data for the Evaporation simulation; simulation will stop.")
! ======================================================================
  end function initializeEvaporationData
