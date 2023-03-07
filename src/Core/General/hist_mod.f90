
module hist_mod

! ==============================================================================
! Gathers histogram data on p, M, E*, A, Z, exciton numbers n, p, pz, and h, and writes tables to file_name
!
! See hist_Init to change bin number, bin min, bin max, step size, and variables that are tallied, etc.
!
! See hist_Print to change file_name output goes to.
!
! For all arrays, the ith variable types are:
!    1 - Momentum (total)
!    2 - Angular momentum (total)
!    3 - E*
!    4 - A number
!    5 - Z number
!    6 - total number of excitons (n)
!    7 - total number of particle excitons (p)
!    8 - total number of charged particle excitons (pz)
!    9 - total number of hole excitons (h)
!
! Leslie M. Kerby, 07/2012, LANL
!
! ==============================================================================

  use, intrinsic:: iso_fortran_env, only: int32, real64

  implicit none
  private

! SUBROUTINE declarations
  public :: hist_tally        !tallies particles in appropriate property bins
  public :: hist_init        !initializes hist arrays to zero, and sets init variables
  public :: hist_print        !prints hist data in file_name


  ! For divide by zero limits:
  real(real64), private, parameter :: div0Lim = 1.0d-15



! CHARACTER variables

! CHARACTER arrays
  character(len=20), private, allocatable, save :: variable_name(:)    !name of the ith variable we are plotting
  character(len=20), private, allocatable, save :: unit_name(:)        !units of the ith variable in histogram

! INTEGER variables
  integer(int32), private, save :: num_var         !number of different variables to plot

! INTEGER arrays
  integer(int32), private, allocatable, save :: num_bins(:)      !number of bins for the ith variable histogram

! REAL variables


! REAL arrays
  real(real64), private, allocatable, save :: momentum_hist1(:)    !momentum histogram tallies
  real(real64), private, allocatable, save :: ang_mom_hist2(:)    !angular momentum histogram tallies
  real(real64), private, allocatable, save :: e_hist3(:)    !e* histogram tallies
  real(real64), private, allocatable, save :: a_hist4(:)    !a number histogram tallies
  real(real64), private, allocatable, save :: z_hist5(:)    !z number histogram tallies
  real(real64), private, allocatable, save :: tot_exciton_hist6(:)    !total excitons histogram tallies
  real(real64), private, allocatable, save :: part_ex_hist7(:)    !total particle-excitons histogram tallies
  real(real64), private, allocatable, save :: proton_ex_hist8(:)    !total charged-particle-excitons histogram tallies
  real(real64), private, allocatable, save :: holes_hist9(:)    !total holes histogram tallies

  real(real64), private, allocatable, save :: particles_out_of_bounds(:)    !array of particles not tallied
  real(real64), private, allocatable, save :: mean(:)            !mean for ith variable type
  real(real64), private, allocatable, save :: step(:)            !step size, or bin width for ith variable type
  real(real64), private, allocatable, save :: minimum(:), maximum(:)    !min and max values for bins for ith variable type



contains



  subroutine hist_tally(momentum_vector, ang_mom_vector, e_residual, a, z, neutron_excitons, &
       & proton_excitons, holes)

! ==============================================================================
!
! Tallies particle data
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), intent(in) :: neutron_excitons,proton_excitons,holes
    integer(int32)    :: excitons_tot, excitons_part

    real(real64), intent(in) :: momentum_vector(3), ang_mom_vector(3), e_residual, a, z

    integer(int32) :: i, bin_number

    real(real64) :: mom_total, ang_mom_tot, w

! ==============================================================================

    excitons_part = neutron_excitons + proton_excitons
    excitons_tot = excitons_part + holes
    mom_total = sqrt((momentum_vector(1)*1000.d0)**2 + (momentum_vector(2)*1000.d0)**2 + &
         & (momentum_vector(3)*1000.d0)**2)    !lmk convert gev to mev
    ang_mom_tot = sqrt(ang_mom_vector(1)**2 + ang_mom_vector(2)**2 + ang_mom_vector(3)**2)
    w = 1.d0        !for simple tallies w=1

 ! Tallying routine
    !Momentum
    i = 1
    if (step(i) < div0Lim .and. step(i) > -div0Lim) then
       step(i) = div0Lim
       write(*,1000) "95"
    end if
    bin_number = ceiling((mom_total-minimum(i))/step(i))
    if (mom_total==minimum(i)) then
       bin_number = 1
    end if
    if (0 < bin_number .and. bin_number <= num_bins(i)) then
       momentum_hist1(bin_number) = momentum_hist1(bin_number) + w
    else
       particles_out_of_bounds(i) = particles_out_of_bounds(i) + w
    end if
    mean(i) = mean(i) + mom_total
    !Angular momentum
    i = 2
    if (step(i) < div0Lim .and. step(i) > -div0Lim) then
       step(i) = div0Lim
       write(*,1000) "111"
    end if
    bin_number = ceiling((ang_mom_tot-minimum(i))/step(i))
    if (ang_mom_tot==minimum(i)) then
       bin_number = 1
    end if
    if (0 < bin_number .and. bin_number <= num_bins(i)) then
       ang_mom_hist2(bin_number) = ang_mom_hist2(bin_number) + w
    else
       particles_out_of_bounds(i) = particles_out_of_bounds(i) + w
    end if
    mean(i) = mean(i) + ang_mom_tot
    !E*
    i = 3
    if (step(i) < div0Lim .and. step(i) > -div0Lim) then
       step(i) = div0Lim
       write(*,1000) "127"
    end if
    bin_number = ceiling((e_residual-minimum(i))/step(i))
    if (e_residual==minimum(i)) then
       bin_number = 1
    end if
    if (0 < bin_number .and. bin_number <= num_bins(i)) then
       e_hist3(bin_number) = e_hist3(bin_number) + w
    else
       particles_out_of_bounds(i) = particles_out_of_bounds(i) + w
    end if
    mean(i) = mean(i) + e_residual
    !A number
    i = 4
    if (step(i) < div0Lim .and. step(i) > -div0Lim) then
       step(i) = div0Lim
       write(*,1000) "143"
    end if
    bin_number = ceiling((a-minimum(i))/step(i))
    if (a==minimum(i)) then
       bin_number = 1
    end if
    if (0 < bin_number .and. bin_number <= num_bins(i)) then
       a_hist4(bin_number) = a_hist4(bin_number) + w
    else
       particles_out_of_bounds(i) = particles_out_of_bounds(i) + w
    end if
    mean(i) = mean(i) + a
    !Z number
    i = 5
    if (step(i) < div0Lim .and. step(i) > -div0Lim) then
       step(i) = div0Lim
       write(*,1000) "159"
    end if
    bin_number = ceiling((z-minimum(i))/step(i))
    if (z==minimum(i)) then
       bin_number = 1
    end if
    if (0 < bin_number .and. bin_number <= num_bins(i)) then
       z_hist5(bin_number) = z_hist5(bin_number) + w
    else
       particles_out_of_bounds(i) = particles_out_of_bounds(i) + w
    end if
    mean(i) = mean(i) + z
    !Total number of excitons
    i = 6
    if (step(i) < div0Lim .and. step(i) > -div0Lim) then
       step(i) = div0Lim
       write(*,1000) "175"
    end if
    bin_number = ceiling((excitons_tot-minimum(i))/step(i))
    if (excitons_tot<minimum(i)+0.01d0) then        !don't want nuclei with no excitons being counted in the #excitons=1 bin
       bin_number = 0
    end if
    if (0 < bin_number .and. bin_number <= num_bins(i)) then
       tot_exciton_hist6(bin_number) = tot_exciton_hist6(bin_number) + w
    else
       particles_out_of_bounds(i) = particles_out_of_bounds(i) + w
    end if
    mean(i) = mean(i) + excitons_tot
    !Total number of particle excitons
    i = 7
    if (step(i) < div0Lim .and. step(i) > -div0Lim) then
       step(i) = div0Lim
       write(*,1000) "191"
    end if
    bin_number = ceiling((excitons_part-minimum(i))/step(i))
    if (excitons_tot<minimum(i)+0.01d0) then        !don't want nuclei with no excitons being counted in the #excitons=1 bin
       bin_number = 0
    end if
    if (0 < bin_number .and. bin_number <= num_bins(i)) then
       part_ex_hist7(bin_number) = part_ex_hist7(bin_number) + w
    else
       particles_out_of_bounds(i) = particles_out_of_bounds(i) + w
    end if
    mean(i) = mean(i) + excitons_part
    !Total number of charged particle excitons
    i = 8
    if (step(i) < div0Lim .and. step(i) > -div0Lim) then
       step(i) = div0Lim
       write(*,1000) "207"
    end if
    bin_number = ceiling((proton_excitons-minimum(i))/step(i))
    if (excitons_tot<minimum(i)+0.01d0) then        !don't want nuclei with no excitons being counted in the #excitons=1 bin
       bin_number = 0
    end if
    if (0 < bin_number .and. bin_number <= num_bins(i)) then
       proton_ex_hist8(bin_number) = proton_ex_hist8(bin_number) + w
    else
       particles_out_of_bounds(i) = particles_out_of_bounds(i) + w
    end if
    mean(i) = mean(i) + proton_excitons
    !Total number of hole excitons
    i = 9
    if (step(i) < div0Lim .and. step(i) > -div0Lim) then
       step(i) = div0Lim
       write(*,1000) "223"
    end if
    bin_number = ceiling((holes-minimum(i))/step(i))
    if (excitons_tot<minimum(i)+0.01d0) then        !don't want nuclei with no excitons being counted in the #excitons=1 bin
       bin_number = 0
    end if
    if (0 < bin_number .and. bin_number <= num_bins(i)) then
       holes_hist9(bin_number) = holes_hist9(bin_number) + w
    else
       particles_out_of_bounds(i) = particles_out_of_bounds(i) + w
    end if
    mean(i) = mean(i) + holes

    return
! ==============================================================================
1000 format("Divide by zero error prevented in 'hist_mod.f90', line(s) ", A)
! ==============================================================================
  end subroutine hist_tally


  subroutine hist_init()

! ==============================================================================
!
! Allocate and initialize hist variables
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    implicit none

    integer(int32) ::  is, i

! ==============================================================================

    num_var = 9

! Allocate arrays of length num_va
    if ( .not.allocated(variable_name) ) then
       allocate (variable_name(num_var), unit_name(num_var), num_bins(num_var), &
            & minimum(num_var), maximum(num_var), step(num_var), mean(num_var), &
            & particles_out_of_bounds(num_var), stat=is)
       if (is/=0) then
          print *, 'error allocating arrays.'
          return
       end if
    end if

    ! Set all values to 0
    mean(:) = 0.0_real64
    particles_out_of_bounds(:) = 0.0_real64

! Momentum (total)
    i = 1
    variable_name(i) = 'Momentum (total)'
    unit_name(i) = 'MeV/c'
    num_bins(i) = 100
    minimum(i) = 600.d0
    maximum(i) = 700.d0
    step(i) = (maximum(i)-minimum(i))/real(num_bins(i))
  !Allocate momentum array
    if ( .not. allocated(momentum_hist1) ) then
       allocate (momentum_hist1(num_bins(i)), stat=is)
       if (is/=0) then
          print *, 'Error allocating momentum array.'
          return
       end if
    end if
  !Set histogram tallies to zero
    momentum_hist1 = 0.d0

! Angular momentum (total)
    i = 2
    variable_name(i) = 'Angular mom (total)'
    unit_name(i) = 'h-bar/(2*pi)'
    num_bins(i) = 20
    minimum(i) = 0.d0
    maximum(i) = 20.d0
    step(i) = (maximum(i)-minimum(i))/real(num_bins(i))
  !Allocate angular momentum array
    if ( .not. allocated(ang_mom_hist2) ) then
       allocate (ang_mom_hist2(num_bins(i)), stat=is)
       if (is/=0) then
          print *, 'Error allocating angular momentum array.'
          return
       end if
    end if
  !Set histogram tallies to zero
    ang_mom_hist2 = 0.d0

! E*
    i = 3
    variable_name(i) = 'E*'
    unit_name(i) = 'MeV'
    num_bins(i) = 100
    minimum(i) = 150.d0
    maximum(i) = 250.d0
    step(i) = (maximum(i)-minimum(i))/real(num_bins(i))
  !Allocate E* array
    if ( .not.allocated(e_hist3) ) then
       allocate (e_hist3(num_bins(i)), stat=is)
       if (is/=0) then
          print *, 'Error allocating E* array.'
          return
       end if
    end if
  !Set histogram tallies to zero
    e_hist3 = 0.d0

! A number
    i = 4
    variable_name(i) = 'A number'
    unit_name(i) = ''
    num_bins(i) = 250
    minimum(i) = 0.d0
    maximum(i) = 250.d0
    step(i) = (maximum(i)-minimum(i))/real(num_bins(i))
  !Allocate A number array
    if ( .not.allocated(a_hist4) ) then
       allocate (a_hist4(num_bins(i)), stat=is)
       if (is/=0) then
          print *, 'Error allocating A number array.'
          return
       end if
    end if
  !Set histogram tallies to zero
    a_hist4 = 0.d0

! Z number
    i = 5
    variable_name(i) = 'Z number'
    unit_name(i) = ''
    num_bins(i) = 100
    minimum(i) = 0.d0
    maximum(i) = 100.d0
    step(i) = (maximum(i)-minimum(i))/real(num_bins(i))
  !Allocate Z number array
    if ( .not.allocated(z_hist5) ) then
       allocate (z_hist5(num_bins(i)), stat=is)
       if (is/=0) then
          print *, 'Error allocating Z number array.'
          return
       end if
    end if
  !Set histogram tallies to zero
    z_hist5 = 0.d0

! Number of excitons (n)
    i = 6
    variable_name(i) = 'Number of excitons (n)'
    unit_name(i) = ''
    num_bins(i) = 200
    minimum(i) = 0.d0
    maximum(i) = 200.d0
    step(i) = (maximum(i)-minimum(i))/real(num_bins(i))
  !Allocate total excitons array
    if ( .not.allocated(tot_exciton_hist6) ) then
       allocate (tot_exciton_hist6(num_bins(i)), stat=is)
       if (is/=0) then
          print *, 'Error allocating total excitons array.'
          return
       end if
    end if
  !Set histogram tallies to zero
    tot_exciton_hist6 = 0.d0

! Number of particle excitons (p)
    i = 7
    variable_name(i) = 'Number of particle excitons (p)'
    unit_name(i) = ''
    num_bins(i) = 100
    minimum(i) = 0.d0
    maximum(i) = 100.d0
    step(i) = (maximum(i)-minimum(i))/real(num_bins(i))
  !Allocate particle excitons array
    if ( .not.allocated(part_ex_hist7) ) then
       allocate (part_ex_hist7(num_bins(i)), stat=is)
       if (is/=0) then
          print *, 'Error allocating particle excitons array.'
          return
       end if
    end if
  !Set histogram tallies to zero
    part_ex_hist7 = 0.d0

! Number of charged particle excitons (pz)
    i = 8
    variable_name(i) = 'Number of charged particle excitons (pz)'
    unit_name(i) = ''
    num_bins(i) = 100
    minimum(i) = 0.d0
    maximum(i) = 100.d0
    step(i) = (maximum(i)-minimum(i))/real(num_bins(i))
  !Allocate charged particle excitons array
    if ( .not.allocated(proton_ex_hist8) ) then
       allocate (proton_ex_hist8(num_bins(i)), stat=is)
       if (is/=0) then
          print *, 'Error allocating charged particle excitons array.'
          return
       end if
    end if
  !Set histogram tallies to zero
    proton_ex_hist8 = 0.d0

! Number of hole excitons (h)
    i = 9
    variable_name(i) = 'Number of hole excitons (h)'
    unit_name(i) = ''
    num_bins(i) = 150
    minimum(i) = 0.d0
    maximum(i) = 150.d0    !
    step(i) = (maximum(i)-minimum(i))/real(num_bins(i))
  !Allocate hole excitons array
    if ( .not.allocated(holes_hist9) ) then
       allocate (holes_hist9(num_bins(i)), stat=is)
       if (is/=0) then
          print *, 'Error allocating hole excitons array.'
          return
       end if
    end if
  !Set histogram tallies to zero
    holes_hist9 = 0.d0

    return
! ==============================================================================
  end subroutine hist_init


  subroutine hist_print(num_inelastic_events)

! ==============================================================================
!
! Prints hist data tables to file
!
! ==============================================================================

    use, intrinsic:: iso_fortran_env, only: int32, int64
    implicit none

    integer(int64), intent(in   ) :: num_inelastic_events

    integer(int32) :: ierror, i, k
!    character(20) :: file_name        !file name for output, not used presently

! ==============================================================================

! File 31 = fres is already opened from cem03.f90
!  OPEN (unit=3, file=file_name, action='write', iostat = ierror)
!  openif: if (ierror==0) then
    !Open was fine. Write data.

    mean = mean/num_inelastic_events

    write (31, '(/a)') '***********************************************************'

    if (num_inelastic_events < div0Lim .and. num_inelastic_events > -div0Lim) then
       write(*,1000) "450, 460, and 470"
    end if

    !Momentum (total)
    i=1
    momentum_hist1 = momentum_hist1/real(num_inelastic_events)
    write (31,1010) variable_name(i), unit_name(i)
    write (31,1020) unit_name(i)
    do k=1,num_bins(i)
       write (31,1100) step(i)*k + minimum(i),momentum_hist1(k)
    end do
    write (31,1030) mean(i)
    write (31,1040) particles_out_of_bounds(i)
    !Angular momentum (total)
    i=2
    ang_mom_hist2 = ang_mom_hist2/real(num_inelastic_events)
    write (31,1010) variable_name(i), unit_name(i)
    write (31,1020) unit_name(i)
    do k=1,num_bins(i)
       write (31,1100) step(i)*k + minimum(i),ang_mom_hist2(k)
    end do
    write (31,1030) mean(i)
    write (31,1040) particles_out_of_bounds(i)
    !E*
    i=3
    e_hist3 = e_hist3/real(num_inelastic_events)
    write (31,1010) variable_name(i), unit_name(i)
    write (31,1020) unit_name(i)
    do k=1,num_bins(i)
       write (31,1100) step(i)*k + minimum(i),e_hist3(k)
    end do
    write (31,1030) mean(i)
    write (31,1040) particles_out_of_bounds(i)
    !A number
    i=4
    a_hist4 = a_hist4/real(num_inelastic_events)
    write (31,1010) variable_name(i), unit_name(i)
    write (31,1020) unit_name(i)
    do k=1,num_bins(i)
       write (31,1100) step(i)*k + minimum(i),a_hist4(k)
    end do
    write (31,1030) mean(i)
    write (31,1040) particles_out_of_bounds(i)
    !Z number
    i=5
    z_hist5 = z_hist5/real(num_inelastic_events)
    write (31,1010) variable_name(i), unit_name(i)
    write (31,1020) unit_name(i)
    do k=1,num_bins(i)
       write (31,1100) step(i)*k + minimum(i),z_hist5(k)
    end do
    write (31,1030) mean(i)
    write (31,1040) particles_out_of_bounds(i)
    !Total number of excitons
    i=6
    tot_exciton_hist6 = tot_exciton_hist6/real(num_inelastic_events)
    write (31,1010) variable_name(i), unit_name(i)
    write (31,1020) unit_name(i)
    do k=1,num_bins(i)
       if (tot_exciton_hist6(k)/=0.d0) then
          write (31,1100) step(i)*k + minimum(i),tot_exciton_hist6(k)
       end if
    end do
    write (31,1030) mean(i)
    write (31,1040) particles_out_of_bounds(i)
    !Total number of particle excitons
    i=7
    part_ex_hist7 = part_ex_hist7/real(num_inelastic_events)
    write (31,1010) variable_name(i), unit_name(i)
    write (31,1020) unit_name(i)
    do k=1,num_bins(i)
       if (part_ex_hist7(k)/=0.d0) then
          write (31,1100) step(i)*k + minimum(i),part_ex_hist7(k)
       end if
    end do
    write (31,1030) mean(i)
    write (31,1040) particles_out_of_bounds(i)
    !Total number of charged particle excitons
    i=8
    proton_ex_hist8 = proton_ex_hist8/real(num_inelastic_events)
    write (31,1010) variable_name(i), unit_name(i)
    write (31,1020) unit_name(i)
    do k=1,num_bins(i)
       if (proton_ex_hist8(k)/=0.d0) then
          write (31,1100) step(i)*k + minimum(i),proton_ex_hist8(k)
       end if
    end do
    write (31,1030) mean(i)
    write (31,1040) particles_out_of_bounds(i)
    !Total number of charged particle excitons
    i=9
    holes_hist9 = holes_hist9/real(num_inelastic_events)
    write (31,1010) variable_name(i), unit_name(i)
    write (31,1020) unit_name(i)
    do k=1,num_bins(i)
       if (holes_hist9(k)/=0.d0) then
          write (31,1100) step(i)*k + minimum(i),holes_hist9(k)
       end if
    end do
    write (31,1030) mean(i)
    write (31,1040) particles_out_of_bounds(i)

1010 format (/,'Normalized histogram data for ', a, ' with units of ', a)
1020 format (a, 'histogram height')
1030 format ('Mean: ', e15.8)
1040 format ('Number of particles not counted in bins: ', e15.8)
1100 format (2e15.8)

!  else openif
!    write (*,*) 'Error opening file.'
!  end if openif
!  CLOSE (unit=3)

!Deallocate arrays
    deallocate (variable_name, unit_name, num_bins, minimum, maximum, step, &
         & particles_out_of_bounds, momentum_hist1, ang_mom_hist2, e_hist3, a_hist4, &
         & z_hist5, tot_exciton_hist6, part_ex_hist7, proton_ex_hist8, holes_hist9, &
         & stat=ierror)
    if (ierror/=0) then
       write (*,*) 'Error deallocating arrays.'
    end if

    return
! ==============================================================================
1000 format("Divide by zero error prevented in 'hist_mod.f90', line(s) ", A)
! ==============================================================================
  end subroutine hist_print

end module hist_mod
