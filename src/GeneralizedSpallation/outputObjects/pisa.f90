! =======================================================================
! | LMK: PISA calculates double differential cross sections and angular |
! | distributions for the fragments listed below. It also prints them.  |
! |                                                                     |
! | The 30 fragment calculations are made for:                          |
! |                     Min &   Max Kinetic energy      (MeV)           |
! |     1 = n           0       160                                     |
! |     2 = p           2       160                                     |
! |     3 = d           2.6     215                                     |
! |     4 = t           3       250                                     |
! |     5 = 3He         2       580                                     |
! |     6 = 4He         2.5     650                                     |
! |     7 = 6He         2.5     2500                                    |
! |     8 = 6Li         4       2500                                    |
! |     9 = 7Li         4       2500                                    |
! |     10 = 8Li        4       2500                                    |
! |     11 = 9Li        4       2500                                    |
! |     12 = 7Be        4.5     2500                                    |
! |     13 = 9Be        4.5     2500                                    |
! |     14 = 10Be       4.5     2500                                    |
! |     15 = 9B         9       2500                                    |
! |     16 = 10B        9       2500                                    |
! |     17 = 11B        9       2500                                    |
! |     18 = 12B        9       2500                                    |
! |     19 = 11C        11      2500                                    |
! |     20 = 12C        11      2500                                    |
! |     21 = 13C        11      2500                                    |
! |     22 = 14C        11      2500                                    |
! |     23 = N          14      2500                                    |
! |     24 = O          16      2500                                    |
! |     25 = F          16      2500                                    |
! |     26 = Ne         16      2500                                    |
! |     27 = Na         16      2500                                    |
! |     28 = Mg         16      2500                                    |
! |     29 = Al         16      2500                                    |
! |     30 = Si         16      2500                                    |
! |                                                                     |
! | Particle type code (for ipar or parz(1,m))                          |
! |     1 = n                                                           |
! |     2 = p                                                           |
! |     3 = d                                                           |
! |     4 = t                                                           |
! |     5 = 3He                                                         |
! |     6 = 4He                                                         |
! |     7 = Pi-                                                         |
! |     8 = Pi0                                                         |
! |     9 = Pi+                                                         |
! |     1000Z + N = A + 999Z, for products heavier than 4He             |
! |                                                                     |
! |     Limitations of the code:                                        |
! |             **Total number of Tk bins must be <=150.                |
! |             Bin 151 is used for particles with energy above the max |
! |             Tk. Bin 152 is used as the total angle-and-energy-      |
! |             integrated cross section (integ only over Tmin to Tmax).|
! |             **dtet=5.0. If this changes, the code will need to be   |
! |             changed (replace 36===>180/dtet, etc.)                  |
! |                                                                     |
! |     KKG, 2006                                                       |
! |     Modified, L. M. Kerby, June 2012.                               |
! =======================================================================

  subroutine PisaSpectra (gsmObj, output, results)
! LMK, Calculates unnormalized particle tallies for d2fr, dofr, yfr


    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro, one, two, thr, four, thousand, radianToDegree

    implicit none
    class(GSM),       intent(inout) :: gsmObj
    class(GSMOutput), intent(inout) :: output
    type(GSMResults), intent(inout) :: results

! LMK, Variable descriptions
    integer(int32) ::  ia, iz        ! a and z of particle m, in int form
    integer(int32) ::  ipar          ! particle m type = parz(1,m), as int, see notes above
    integer(int32) ::  itet          ! = (theta of particle m)/dtet + 1, gives absolute theta bin number
                                        !For d2fr we do not tally all theta bins, only
                                        !those of the ten angles in ang()
    integer(int32) ::  kfr           ! particle fragment number (1-30) of particle m
    integer(int32) ::  kt            ! tk bin number of particle m
    integer(int32) ::  ka            ! counting variable going through the 10 theta bins we are tallying
    integer(int32) ::  m             ! counting variable going through the 150 particles in spt()
    integer(int32) ::  k             ! counting variable going through the 30 possible fragment types
    real(real64) ::  d2fr             ! bin count of all events in specific theta and tk bins, ie
                                        ! unnormalized Double Differential cross section
    real(real64) ::  dofr             ! bin count of all events in a given theta bin, irrespective of tk
                                        ! ie unnormalized Angular Distribution
    real(real64) ::  yfr              ! number of particles tallied:
                                        ! yfr(1,kfr) = total # of particles tallied,
                                        ! yfr(2,kfr) = # of particles with Tmin<=Tk<=Tmax
!    real(real64) ::  ang              ! ten user input angles to calculate d2fr at
!    real(real64) ::  dang             ! delta angle. for d2fr: theta bin width is 2*dang (ang +/- dang)
!                                        ! For dofr: width of theta bin is dang
!                                      ! Migrated to GSMOutput [see pisaAngles, pisaDTheta]
    real(real64) ::  t                ! kinetic energy (mev) of particle m = parz(2,m)
    real(real64) ::  tet              ! theta (degrees) of particle m = parz(3,m)
    real(real64) ::  am, zm           ! a and z of particle m, zm = spt(4,m)

    common /spefra/  d2fr(30,11,152), dofr(30,37), yfr(2,30)

! ======================================================================

    ! a and z of the 30 different particles to be considered
    integer(int32), parameter, dimension(30) :: izfr = [ &
         & 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, &
         & 5, 5, 5, 6, 6, 6, 6, 7, 8, 9,10,11,12,13,14]
    integer(int32), parameter, dimension(30) :: iafr = [ &
         & 1, 1, 2, 3, 3, 4, 6, 6, 7, 8, 9, 7, 9,10, 9, &
         & 10,11,12,11,12,13,14,14,16,18,20,22,24,26,28]
    ! min and max kinetic energy to consider, for each of the
    ! 30 different particles (neutrons through C14)
    real(real64), parameter, dimension(30) :: tmin =  [ & 
         &  0.0d0, 2.0d0, 2.6d0, 3.0d0, 2.0d0, 2.5d0, &
         & 2.5d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.5d0,  &
         & 4.5d0, 4.5d0, 9.0d0, 9.0d0, 9.0d0, 9.0d0,  &
         & 11.0d0,11.0d0,11.0d0,11.0d0,14.0d0,16.0d0, &
         & 16.0d0,16.0d0,16.0d0,16.0d0,16.0d0,16.0d0  ]
    real(real64), parameter, dimension(30) :: tmax = [ &
         & 160.d0,160.d0,215.d0,250.d0,580.d0,650.d0, &
         & 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,  &
         & 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,  &
         & 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,  &
         & 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3   ]
    ! delta theta for bins, set as a constant
    real(real64), parameter :: dtet = 5.0_real64
    real(real64), parameter, dimension(4) :: &
         &  tgr = [50.0d0, 350.0d0,  950.0d0, 2600.0d0], & ! defines break points in kinetic energy, for bins
         & dtgr = [ 2.0d0,   5.0d0,   20.0d0,   50.0d0]    ! defines size of tk bins between tgr(n-1) and tgr(n)
    integer(int32), parameter, dimension(4) :: &
         & jgr  = [ 1,      26,       86,      116    ] ! gives the index of the first bin in the nth interval

! ======================================================================

    do m = 1, results%numProgeny
! LMK, First find A and Z for particle m
       am = zro
       zm = zro
       if (results%progenyBnk(m)%restMass < 0.001d0) go to 40
                ! LMK, results%progenyBnk(m)%restMass is the rest mass of particle m. I think maybe this is designed to
                ! be an error catcher, which ends the calculation if mass is unphysical? <----****
       ipar = nint(results%progenyBnk(m)%typeID)
       zm = results%progenyBnk(m)%numProtons

! LMK, Section assign A number using algorithm from particle type code
       if (ipar < 3)     then
          am = one
       elseif (ipar < 4) then
          am = two
       elseif (ipar < 6) then
          am = thr
       elseif (ipar == 6) then
          am = four
       elseif(ipar >= 7.and.ipar <= 9)  then
          go  to  30
       elseif (ipar > 9) then
!   Modified calculation of mass number of residual nucleus:
          am = results%progenyBnk(m)%typeID - 999.d0*zm
       end if

       tet = results%progenyBnk(m)%theta*radianToDegree   ! deg. lmk, copy particle m theta angle
       t   = results%progenyBnk(m)%kinEnergy * thousand   ! in MeV
       iz = nint(zm)
       ia = nint(am)
! LMK, Section identifies what fragment number particle m is, and skips the loop
! and goes to the next m=m+1 particle if the particle is heavier than 28Si
! Assigns particle fragment number to kfr
       kfr = 0
       do k=1,30
          if (k <= 23)  then            !lmk, if particle type<23 then a and z have to match
              if (iz == izfr(k) .and. ia == iafr(k)) kfr = k
           else                  !lmk, if particle type>=23 then only z has to match (no isotopes)
              if(iz == izfr(k)) kfr = k
          end if
       end do
       if (kfr == 0)  go  to  30

! LMK, Finds the Tk bin of particle m and assigns it to kt
       if (t < tgr(1))      then
          kt = int(t/dtgr(1)) + jgr(1)
       elseif (t < tgr(2))  then
          kt = int((t - tgr(1))/dtgr(2)) + jgr(2)
       elseif (t < tgr(3))  then
          kt = int((t - tgr(2))/dtgr(3)) + jgr(3)
       elseif (t < tgr(4))  then
          kt = int((t - tgr(3))/dtgr(4)) + jgr(4)
       else
          kt = 151      ! lmk, bad form as this number will change depending on tk bins and
                        ! break points. This should be calculated!      <------******
       end if
! LMK, Determine if theta of particle m is in one of the theta bins we are tallying, and if so
        ! add particle m to associated Tk and theta bins in d2fr
       do ka=1,10
          if(tet >= (output%pisaAngles(ka)-output%pisaDTheta) .and. &
               & tet <= (output%pisaAngles(ka)+output%pisaDTheta)) then
             d2fr(kfr,ka,kt) = d2fr(kfr,ka,kt) + one
             if(t >= tmin(kfr).and.t <= tmax(kfr)) &
                  & d2fr(kfr,ka,152) = d2fr(kfr,ka,152) + one
          end if
       end do
! LMK, d2fr(:,11,:) is the angle-integrated cross section (mb/MeV) for each Tk bin kt
       d2fr(kfr,11,kt) = d2fr(kfr,11,kt) + one         !lmk, add particle m to the correct tk bin
! LMK, d2fr(:,11,152) is the angle-and-energy-integrated cross section (mb) for each particle
        !fragment number, that has min<=Tk<=max. Should match dofr(kfr,37).
       if (t >= tmin(kfr).and.t <= tmax(kfr)) &
            & d2fr(kfr,11,152) = d2fr(kfr,11,152) + one
       itet = int(tet/dtet) + 1        !lmk, find absolute theta bin number itet
       if (itet > 36) itet = 36       !lmk, assuming dtet=5, this corresponds to tet>180*. this
                                        !is probably an error control statement. Again, this "36"
                        !is data specific and it would be better to use 180/dtet. <-----*****
! LMK, Add event in appropriate absolute theta bin itet, in dofr(kfr, itet), used for angular dist.
       dofr(kfr,itet) = dofr(kfr,itet) + one
! LMK, dofr(kfr, 37) is the angle-and-energy-integrated cross section (mb) for each fragment number
        ! that has min<=Tk<=max. Should match d2fr(:,11,152).
       if (t >= tmin(kfr).and.t <= tmax(kfr)) &
            & dofr(kfr,37) = dofr(kfr,37) + one
       yfr(1,kfr) = yfr(1,kfr) + one
! LMK, Add event to yfr(2,kfr) if particle m has energy Tmin<=Tk<=Tmax
       if (t >= tmin(kfr).and.t <= tmax(kfr)) &
            & yfr(2,kfr) = yfr(2,kfr) + one
30     continue
    end do
40  continue
    return
! ======================================================================
  end subroutine PisaSpectra



  subroutine pisaInit (gsmObj)

! ======================================================================
! LMK, Reads ang(10) and dang from the input file.
! LMK, Initializes values of d2fr, dofr, yfr to zero.
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64
    use gsm_params, only: zro

    implicit none
    class(GSM),       intent(inout) :: gsmObj


! LMK, Variable descriptions
    integer(int32) ::  kt            ! counting variable for the 152 tk bins
    integer(int32) ::  ka            ! counting variable for the 37 (=180/(dtet=5) + 1) theta bins
    integer(int32) ::  kf            ! counting variable for the 30 particle fragment types
    real(real64) ::  d2fr             ! bin count of all events in specific theta and tk bins, ie
                                        ! unnormalized Double Differential cross section
    real(real64) ::  dofr             ! bin count of all events in a given theta bin, irrespective of tk
                                        ! ie unnormalized Angular Distribution
    real(real64) ::  yfr              ! number of particles tallied:
                                        ! yfr(1,kfr) = total # of particles tallied,
                                        ! yfr(2,kfr) = # of particles with Tmin<=Tk<=Tmax
!    real(real64) ::  ang              ! ten user input angles to calculate d2fr at
!    real(real64) ::  dang             ! delta angle, for d2fr: theta bin width is 2*dang (ang +/- dang)
!                                        ! for dofr: width of theta bin is dang
!                                      ! Migrated to GSMOutput [see pisaAngles, pisaDTheta]

    common /spefra/  d2fr(30,11,152), dofr(30,37), yfr(2,30)

! ======================================================================

    do kf=1,30
       yfr(1,kf) = zro
       yfr(2,kf) = zro
       do ka=1,37
          dofr(kf,ka) = zro
          if(ka <= 11) then
             do kt=1,152
                d2fr(kf,ka,kt) = zro
             end do
          end if
       end do
    end do
    return
! ======================================================================
  end subroutine pisainit


  subroutine pisaprint (gsmObj, output, sigin, ncas)

! ======================================================================
! LMK, Normalizes d2fr, dofr and finishes double differential cross section and angular
! distribution calculations, and then prints them.
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, real64, int64
    use gsm_params, only: zro, radianToDegree, pi, twpi

    implicit none
    class(GSM),       intent(inout) :: gsmObj
    class(GSMOutput), intent(in   ) :: output
    real(real64),     intent(in   ) :: sigin   ! inelastic cross section (mb) (total, across all events?)
    integer(int64),   intent(in   ) :: ncas    ! number of inelastic events simulated


! LMK, Variable descriptions
    integer(int32) ::  ia            ! counting variable for the 10 theta bins in d2fr
    integer(int32) ::  ka            ! counting variable for the 37 (=180/(dtet=5) + 1) theta bins
    integer(int32) ::  kfr           ! particle fragment number (1-30) of particle m
    integer(int32) ::  it            ! counting variable for the 150 tk bins
    real(real64) ::  fn               ! # of inelastic events in double
    real(real64) ::  fac              ! = inelastic cross section / number of inelastic events
    real(real64) ::  ang1, ang2       ! lower and upper angles used to find dom and domp
    real(real64) ::  dom(36)          ! delta omega (sr) for the 36 theta bins in dofr calculations
    real(real64) ::  domp(10)            ! delta omega (sr) for the 10 theta bins in d2fr calculations
    real(real64) ::  c(11)            ! double differential cross section of the ia-th theta bin and
                                        ! it-th Tk bin (mb/MeV/sr)
    real(real64) ::  sfr1             ! total cross section for all energies (mb) (= yfr(1,kfr)*fac)
                        ! Cross section of fragment k = (N_k)*(Sigma_inelastic)/(N_inelastic)
    real(real64) ::  sfr2             ! total cross section for tmin<=tk<=tmax (mb) (= yfr(2,kfr)*fac)
    real(real64) ::  dtk              ! bin width of the it-th tk bin
    real(real64) ::  tk1, tk2         ! lower and upper energies of the it-th tk bin
    real(real64) ::  summ             ! used to see if there were any particle tallies in the tk bin, and
                                        ! therefore whether or not to print x-section info
    real(real64) ::  d2fr             ! bin count of all events in specific theta and tk bins, ie
                                        ! unnormalized Double Differential cross section
    real(real64) ::  dofr             ! bin count of all events in a given theta bin, irrespective of tk
                                        ! ie unnormalized Angular Distribution
    real(real64) ::  yfr              ! number of particles tallied:
                                        ! yfr(1,kfr) = total # of particles tallied,
                                        ! yfr(2,kfr) = # of particles with Tmin<=Tk<=Tmax
!    real(real64) ::  ang              ! ten user input angles to calculate d2fr at
!    real(real64) ::  dang             ! delta angle. for d2fr: theta bin width is 2*dang (ang +/- dang)
!                                        ! For dofr: width of theta bin is dang
!                                      ! Migrated to GSMOutput [see pisaAngles, pisaDTheta]

! ======================================================================
 
   common /spefra/  d2fr(30,11,152), dofr(30,37), yfr(2,30)

   character(LEN=4), parameter, dimension(30) :: fname =  [ &
        & "   n", "   p", "   d", "   t", " He3", " He4", &
        & " He6", " Li6", " Li7", " Li8", " Li9", " Be7", &
        & " Be9", "Be10", "  B9", " B10", " B11", " B12", &
        & " C11", " C12", " C13", " C14", " Z=7", " Z=8", &
        & " Z=9", "Z=10", "Z=11", "Z=12", "Z=13", "Z=14"  ]
   ! min and max kinetic energy to consider, for each of the
   ! 30 different particles (neutrons through C14)
   real(real64), parameter, dimension(30) :: tmin =  [ & 
         &  0.0d0, 2.0d0, 2.6d0, 3.0d0, 2.0d0, 2.5d0, &
         & 2.5d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.5d0,  &
         & 4.5d0, 4.5d0, 9.0d0, 9.0d0, 9.0d0, 9.0d0,  &
         & 11.0d0,11.0d0,11.0d0,11.0d0,14.0d0,16.0d0, &
         & 16.0d0,16.0d0,16.0d0,16.0d0,16.0d0,16.0d0  ]
    real(real64), parameter, dimension(30) :: tmax = [ &
         & 160.d0,160.d0,215.d0,250.d0,580.d0,650.d0, &
         & 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,  &
         & 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,  &
         & 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,  &
         & 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3   ]
    ! delta theta for bins, set as a constant
    real(real64), parameter :: dtet = 5.0_real64
    real(real64), parameter, dimension(4) :: &
         &  tgr = [50.0d0, 350.0d0,  950.0d0, 2600.0d0], & ! defines break points in kinetic energy, for bins
         & dtgr = [ 2.0d0,   5.0d0,   20.0d0,   50.0d0]    ! defines size of tk bins between tgr(n-1) and tgr(n)
    integer(int32), parameter, dimension(4) :: &
         & jgr  = [ 1,      26,       86,      116    ] ! gives the index of the first bin in the nth interval


! ======================================================================

    fn = dble(ncas)
    fac = sigin/fn
! LMK, Calculates delta omega: domp for double differential and dom for angular distributions
        ! delta omega = 2*Pi*[cos(theta1) - cos(theta2)]
    do ka=1,36
! LMK, calculates domp for the 10 theta bins in double differential cross section calculations
       if (ka <= 10) then
          ang1 = (output%pisaAngles(ka) - output%pisaDTheta)/radianToDegree
          if(ang1 < zro) ang1 = zro
          ang2 = (output%pisaAngles(ka) + output%pisaDTheta)/radianToDegree
          if(ang2 > pi) ang2 = pi
          domp(ka) = twpi*(cos(ang1) - cos(ang2))
       end if
! LMK, calculates dom for the 36 theta bins in angular distribution calculations
       ang1 = dble(ka-1)*dtet/radianToDegree
       ang2 = dble(ka)*dtet/radianToDegree
       dom(ka) = twpi*(cos(ang1) - cos(ang2))
    end do

! LMK, Main double differential cross section loop.....finishes calculations and prints
! LMK, Loop over 30 different fragment types
    do kfr=1,30
       if(yfr(1,kfr) <= zro) go to 50  ! lmk, if there are no particles tallied for fragment type
                                        ! kfr, skip to next fragment type.
! LMK, First calculate and print total cross section
       sfr1 = yfr(1,kfr)*fac           ! lmk, total cross section for fragment kfr
       sfr2 = yfr(2,kfr)*fac           ! lmk, total x section for fragment kfr in tmin<=tk<=tmax
       write (31, 101) fname(kfr), sfr1, tmin(kfr), tmax(kfr), &
            & sfr2, output%pisaAngles
101    format(//25x,'Double differential cross-section d2S/dTdO ', &
            & '(mb/MeV/sr) of ',a4/12x, &
            & 'prod. xsec for all energies=',1pe11.4,' mb, for ', &
            & 0pf5.0,'<T<',f5.0,' xsec=',1pe11.4,'mb'/2x, &
            & 'T(MeV)/angle:',9(0pf5.0,6x),f5.0,4x,'dS/dT(mb/MeV)')
! LMK, Loop over 150 Tk bins
       do it=1,150
          summ = zro            !lmk, used to see if there are any particles to print x section for
! LMK, First determine the bin width and lower/upper energies for each Tk bin
! LMK, Check to see which of the four intervals the it-th Tk bin is in, and assign dtk and tk1
          if (it >= jgr(1) .and. it < jgr(2))            then
             dtk = dtgr(1)               ! lmk, bin width
             tk1 = dble(it-jgr(1))*dtk   ! lmk, lower energy
          elseif (it >= jgr(2) .and. it < jgr(3))        then
             dtk = dtgr(2)
             tk1 = dble(it-jgr(2))*dtk + tgr(1)
          elseif (it >= jgr(3) .and. it < jgr(4))        then
             dtk = dtgr(3)
             tk1 = dble(it-jgr(3))*dtk + tgr(2)
          elseif (it >= jgr(4))                           then
             dtk = dtgr(4)
             tk1 = dble(it-jgr(4))*dtk + tgr(3)
          else
          end if
          tk2 = tk1 + dtk               ! lmk, upper energy
! LMK, Loop over 10 theta bins
! LMK, Calculate double diff cross section for each of the ten theta bins, using d2fr tallies
          do ia=1,10
             c(ia) = d2fr(kfr,ia,it)*fac/(domp(ia)*dtk)  !lmk, (mb/MeV/sr)
        ! LMK, c(ia) = Double diff x section at ia-th theta bin, it-th Tk bin, for kfr fragment
             summ = summ + c(ia)
          end do
          c(11) = d2fr(kfr,11,it)*fac/dtk       !lmk, double diff x section integrated over angle
                                                ! (mb/MeV)
          summ = summ + c(11)
          if(summ > zro) write (31, 102) tk1, tk2, c
102       format(f5.0,'-',f5.0,11(1pe11.4))
       end do
! LMK, Calculate cross sections for all particles Tmin<=Tk<=Tmax in the 10 theta bins
       do ia=1,10
          c(ia) = d2fr(kfr,ia,152)*fac/domp(ia) !lmk, dd x section integrated over energy (mb/sr)
       end do
       c(11) = d2fr(kfr,11,152)*fac            !lmk, total cross section (mb)
       write (31, 103) c
103    format(' Energ. int.',11(1pe11.4))
50     continue
    end do
! LMK, Done with double differential cross sections

! LMK, Now calculate and print angular distributions for the first 10 fragment types (through 8Li)
    write (31, 104) (tmin(kfr), tmax(kfr),kfr=1,10), &
         & (fname(kfr),kfr=1,10)
104 format(//10x,'Angular distribution of produced fragments', &
         & ' dS/dOm [mb/sr] for energy range(MeV)'/'Tmin-Tmax', &
         & 10(f5.0,'-',f5.0)/ 'ang1-ang2',10(3x,a4,4x))
! LMK, Loop over 36 theta bins
    do ka=1,36
       ang1 = dble(ka-1)*dtet         !lmk, calculate lower and upper angles
       ang2 = dble(ka)*dtet
       do kfr=1,10
          c(kfr) = dofr(kfr,ka)*fac/dom(ka)     !lmk, angular distribution of fragment kfr,
                                                ! ka-th theta bin (mb/sr)
       end do
       write (31, 105) ang1, ang2, (c(kfr),kfr=1,10)
105    format(f4.0,'-',f4.0,10(1pe11.4))
    end do
! LMK, Print total cross section for fragment kfr
    do kfr=1,10
       c(kfr) = dofr(kfr,37)*fac
    end do
    write (31, 106) (c(kfr),kfr=1,10)
106 format(' Int. x sec',10(1pe11.4))
! LMK, Now calculate and print angular distributions for the 10-20 fragment types (9Li-12C)
    write (31, 107) (tmin(kfr), tmax(kfr),kfr=11,20), &
         & (fname(kfr),kfr=11,20)
107 format(//10x,'Angular distribution of produced fragments', &
         & ' dS/dOm [mb/sr] for energy range(MeV)'/'Tmin-Tmax', &
         & 10(f5.0,'-',f5.0)/'ang1-ang2',10(3x,a4,4x))
! LMK, Loop over 36 theta bins
    do ka=1,36
       ang1 = dble(ka-1)*dtet         !lmk, calculate lower and upper angles
       ang2 = dble(ka)*dtet
       do kfr=11,20
          c(kfr-10) = dofr(kfr,ka)*fac/dom(ka)  !lmk, angular distribution of fragment kfr,
                                                ! ka-th theta bin (mb/sr)
       end do
       write(31,108) ang1,ang2,(c(kfr),kfr=1,10)
108    format(f4.0,'-',f4.0,10(1pe11.4))
    end do
! LMK, Print total cross section for fragment kfr
    do kfr=11,20
       c(kfr-10) = dofr(kfr,37)*fac
    end do
    write (31, 109) (c(kfr),kfr=1,10)
109 format(' Int. xsec',10(1pe11.4))
! LMK, Now calculate and print angular distributions for the last 10 fragment types (13C-Z14/28Si)
    write (31, 110) (tmin(kfr), tmax(kfr),kfr=21,30), &
         & (fname(kfr),kfr=21,30)
110 format(//10x,'Angular distribution of produced fragments', &
         & ' dS/dOm [mb/sr] for energy range(MeV)'/'Tmin-Tmax', &
         & 10(f5.0,'-',f5.0)/'ang1-ang2',10(3x,a4,4x))
! LMK, Loop over 36 theta bins
    do ka=1,36
       ang1 = dble(ka-1)*dtet         !lmk, calculate lower and upper angles
       ang2 = dble(ka)*dtet
       do kfr=21,30
          c(kfr-20) = dofr(kfr,ka)*fac/dom(ka)  !lmk, angular distribution of fragment kfr,
                                                ! ka-th theta bin (mb/sr)
       end do
       write (31, 111) ang1, ang2, (c(kfr),kfr=1,10)
111    format(f4.0,'-',f4.0,10(1pe11.4))
    end do
! LMK, Print total cross section for fragment kfr
    do kfr=21,30
       c(kfr-20) = dofr(kfr,37)*fac
    end do
    write (31, 112) (c(kfr),kfr=1,10)
112 format(' Int. xsec',10(1pe11.4))
    return

! ======================================================================
  end subroutine pisaprint
