c =======================================================================
c | LMK: PISA calculates double differential cross sections and angular |
c | distributions for the fragments listed below. It also prints them.	|
c |									|
c | The 30 fragment calculations are made for:				|
c |			Min &	Max Kinetic energy	(MeV)		|
c |	1 = n		0	160					|
c |	2 = p		2	160					|
c |	3 = d		2.6	215					|
c |	4 = t		3	250					|
c |	5 = 3He		2	580					|
c |	6 = 4He		2.5	650					|
c |	7 = 6He		2.5	2500					|
c |	8 = 6Li		4	2500					|
c |	9 = 7Li		4	2500					|
c |	10 = 8Li	4	2500					|
c |	11 = 9Li	4	2500					|
c |	12 = 7Be	4.5	2500					|
c |	13 = 9Be	4.5	2500					|
c |	14 = 10Be	4.5	2500					|
c |	15 = 9B		9	2500					|
c |	16 = 10B	9	2500					|
c |	17 = 11B	9	2500					|
c |	18 = 12B	9	2500					|
c |	19 = 11C	11	2500 					|
c |	20 = 12C	11	2500					|
c |	21 = 13C	11	2500					|
c |	22 = 14C	11	2500					|
c |	23 = N		14	2500					|
c |	24 = O		16	2500					|
c |	25 = F		16	2500					|
c |	26 = Ne		16	2500					|
c |	27 = Na		16	2500					|
c |	28 = Mg		16	2500					|
c |	29 = Al		16	2500					|
c |	30 = Si		16	2500					|
c |									|
c | Particle type code (for ipar or parz(1,m))				|
c |	1 = n								|
c |	2 = p								|
c |	3 = d								|
c |	4 = t								|
c |	5 = 3He								|
c |	6 = 4He								|
c |	7 = Pi-								|
c |	8 = Pi0								|
c |	9 = Pi+								|
c |	1000Z + N = A + 999Z, for products heavier than 4He		|
c |									|
c |	Limitations of the code: 					|
c |		**Total number of Tk bins must be <=150.		|
c |		Bin 151 is used for particles with energy above the max	|
c |		Tk. Bin 152 is used as the total angle-and-energy-	|
c |		integrated cross section (integ only over Tmin to Tmax).|
c |		**dtet=5.0. If this changes, the code will need to be	|
c |		changed (replace 36===>180/dtet, etc.)			|
c |									|
c |	KKG, 2006							|
c |	Modified, L. M. Kerby, June 2012.				|
c =======================================================================

      subroutine PisaSpectra ()
c LMK, Calculates unnormalized particle tallies for d2fr, dofr, yfr

c LMK, Variable descriptions
     	integer*4 ia, iz	! A and Z of particle m, in INT form
	integer*4 iafr, izfr	! A and Z of the 30 different particles to be considered
	integer*4 ipar		! Particle m type = parz(1,m), as INT, see notes above
	integer*4 itet		! = (theta of particle m)/dtet + 1, gives absolute theta bin number
					!For d2fr we do not tally all theta bins, only
					!those of the ten angles in ang() 
	integer*4 kfr		! Particle fragment number (1-30) of particle m
	integer*4 kt		! Tk bin number of particle m
	integer*4 ka		! Counting variable going through the 10 theta bins we are tallying
	integer*4 jgr		! Gives the index of the first bin in the nth interval
	integer*4 m		! Counting variable going through the 150 particles in spt()
	integer*4 k		! Counting variable going through the 30 possible fragment types
	real*8 d2fr		! Bin count of all events in specific theta and Tk bins, ie 
					! unnormalized Double Differential cross section
	real*8 dofr		! Bin count of all events in a given theta bin, irrespective of Tk
					! ie unnormalized Angular Distribution
	real*8 yfr		! Number of particles tallied: 
					! yfr(1,kfr) = total # of particles tallied, 
					! yfr(2,kfr) = # of particles with Tmin<=Tk<=Tmax
	real*8 degrad		! Converts from radians to degrees (I think)
	real*8 pi, twpi		! Constants Pi and 2*Pi
	real*8 ang		! Ten user input angles to calculate d2fr at
	real*8 dang		! Delta angle. For d2fr: theta bin width is 2*dang (ang +/- dang)
					! For dofr: width of theta bin is dang
	real*8 spt, parz	! Particle information tallies, used all through CEM03 program
	real*8 tgr		! Defines break points in kinetic energy, for bins
	real*8 dtgr		! Defines size of Tk bins between tgr(n-1) and tgr(n)
	real*8 tmin, tmax	! Min and max kinetic energy to consider, for each of the 
					! 30 different particles (neutrons through C14)
	real*8 t		! Kinetic Energy (MeV) of particle m = parz(2,m)
	real*8 tet		! Theta (degrees) of particle m = parz(3,m)
	real*8 dtet		! Delta theta for bins, set as a constant
	real*8 am, zm		! A and Z of particle m, zm = spt(4,m)

      common /blok77/  spt(5,150)
      common /zapp/    parz(6,150)
      common /spefra/  d2fr(30,11,152), dofr(30,37), yfr(2,30),
     &		ang(10), dang		! LMK, added ang and dang to spefra	06/2012
      common /degrad/  degrad		
      common /pi/      pi, twpi		

      dimension tgr(4), dtgr(4), jgr(4)
      dimension izfr(30), iafr(30), tmin(30), tmax(30)
      data izfr/ 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5,
     &           5, 5, 5, 6, 6, 6, 6, 7, 8, 9,10,11,12,13,14/
      data iafr/ 1, 1, 2, 3, 3, 4, 6, 6, 7, 8, 9, 7, 9,10, 9,
     &          10,11,12,11,12,13,14,14,16,18,20,22,24,26,28/
      data tmin/  0.0d0, 2.0d0, 2.6d0, 3.0d0, 2.0d0, 2.5d0, 
     &            2.5d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.5d0,
     &            4.5d0, 4.5d0, 9.0d0, 9.0d0, 9.0d0, 9.0d0,
     &           11.0d0,11.0d0,11.0d0,11.0d0,14.0d0,16.0d0,
     &           16.0d0,16.0d0,16.0d0,16.0d0,16.0d0,16.0d0/ 
      data tmax/ 160.d0,160.d0,215.d0,250.d0,580.d0,650.d0, 
     &            2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,
     &            2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,
     &            2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,
     &            2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3/
      data dtet /5.0d0/  
      data  tgr/50.0d0, 350.0d0, 950.0d0, 2600.0d0/,
     &     dtgr/ 2.0d0,   5.0d0,   20.0d0,   50.0d0/,
     &      jgr/     1,      26,       86,      116/     
      data zro, one, two, thr, for, mil /0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 
     &                             1.d3/

c ======================================================================
      
        do m = 1,150 
c LMK, First find A and Z for particle m
        am = zro
        zm = zro
        if (spt(5,m) < 0.001d0) go to 40
		! LMK, spt(5,m) is the rest mass of particle m. I think maybe this is designed to
		! be an error catcher, which ends the calculation if mass is unphysical? <----****
        ipar = nint(parz(1,m))		!LMK, Copy particle type code
        zm = spt(4,m)			!LMK, Copy Z number

c LMK, Section assign A number using algorithm from particle type code
        if (ipar < 3)     then
          am = one
        elseif (ipar < 4) then
          am = two
        elseif (ipar < 6) then
          am = thr
        elseif (ipar == 6) then
          am = for
        elseif(ipar >= 7.and.ipar <= 9)  then
          go  to  30
        elseif (ipar > 9) then
c   Modified calculation of mass number of residual nucleus:
          am = parz(1,m) - 999.d0*zm 
        endif

        tet = parz(3,m)*degrad      ! deg. LMK, Copy particle m theta angle
        t   = parz(2,m)*mil         ! MeV LMK, Copy particle m kinetic energy
        iz = nint(zm)
        ia = nint(am)
c LMK, Section identifies what fragment number particle m is, and skips the loop 
c and goes to the next m=m+1 particle if the particle is heavier than 28Si
c Assigns particle fragment number to kfr
        kfr = 0 
          do k=1,30
          if (k <= 23)  then		!LMK, If particle type<23 then A and Z have to match
            if (iz == izfr(k) .and. ia == iafr(k)) kfr = k
          else			!LMK, if particle type>=23 then only Z has to match (no isotopes)
            if(iz == izfr(k)) kfr = k
          endif
          enddo
        if (kfr == 0)  go  to  30

c LMK, Finds the Tk bin of particle m and assigns it to kt
        if (t < tgr(1))      then
          kt = int(t/dtgr(1)) + jgr(1)
        elseif (t < tgr(2))  then
          kt = int((t - tgr(1))/dtgr(2)) + jgr(2)
        elseif (t < tgr(3))  then
          kt = int((t - tgr(2))/dtgr(3)) + jgr(3)
        elseif (t < tgr(4))  then
          kt = int((t - tgr(3))/dtgr(4)) + jgr(4)
        else
          kt = 151	! LMK, Bad form as this number will change depending on Tk bins and
			! break points. This should be calculated!	<------******
        endif
c LMK, Determine if theta of particle m is in one of the theta bins we are tallying, and if so
	! add particle m to associated Tk and theta bins in d2fr
          do ka=1,10
          if(tet >= (ang(ka)-dang).and.tet <= (ang(ka)+dang)) then
            d2fr(kfr,ka,kt) = d2fr(kfr,ka,kt) + one
            if(t >= tmin(kfr).and.t <= tmax(kfr))
     &      d2fr(kfr,ka,152) = d2fr(kfr,ka,152) + one
          endif
          enddo
c LMK, d2fr(:,11,:) is the angle-integrated cross section (mb/MeV) for each Tk bin kt
        d2fr(kfr,11,kt) = d2fr(kfr,11,kt) + one		!LMK, Add particle m to the correct Tk bin
c LMK, d2fr(:,11,152) is the angle-and-energy-integrated cross section (mb) for each particle 
	!fragment number, that has min<=Tk<=max. Should match dofr(kfr,37).
        if (t >= tmin(kfr).and.t <= tmax(kfr))
     &    d2fr(kfr,11,152) = d2fr(kfr,11,152) + one
        itet = int(tet/dtet) + 1	!LMK, Find absolute theta bin number itet
        if (itet > 36) itet = 36	!LMK, Assuming dtet=5, this corresponds to tet>180*. This
					!is probably an error control statement. Again, this "36"
			!is data specific and it would be better to use 180/dtet. <-----*****
c LMK, Add event in appropriate absolute theta bin itet, in dofr(kfr, itet), used for angular dist.
        dofr(kfr,itet) = dofr(kfr,itet) + one	
c LMK, dofr(kfr, 37) is the angle-and-energy-integrated cross section (mb) for each fragment number 
	! that has min<=Tk<=max. Should match d2fr(:,11,152).
        if (t >= tmin(kfr).and.t <= tmax(kfr))
     &    dofr(kfr,37) = dofr(kfr,37) + one
        yfr(1,kfr) = yfr(1,kfr) + one
c LMK, Add event to yfr(2,kfr) if particle m has energy Tmin<=Tk<=Tmax
        if (t >= tmin(kfr).and.t <= tmax(kfr))
     &    yfr(2,kfr) = yfr(2,kfr) + one
   30   continue
        enddo
   40 continue
      return
	end
c

c 
      subroutine PisaInit ()
c LMK, Reads ang(10) and dang from the input file.	
c LMK, Initializes values of d2fr, dofr, yfr to zero. 	

c LMK, Variable descriptions
	integer*4 kt		! Counting variable for the 152 Tk bins 
	integer*4 ka		! Counting variable for the 37 (=180/(dtet=5) + 1) theta bins
	integer*4 kf		! Counting variable for the 30 particle fragment types
	real*8 d2fr		! Bin count of all events in specific theta and Tk bins, ie 
					! unnormalized Double Differential cross section
	real*8 dofr		! Bin count of all events in a given theta bin, irrespective of Tk
					! ie unnormalized Angular Distribution
	real*8 yfr		! Number of particles tallied: 
					! yfr(1,kfr) = total # of particles tallied, 
					! yfr(2,kfr) = # of particles with Tmin<=Tk<=Tmax
	real*8 ang		! Ten user input angles to calculate d2fr at
	real*8 dang		! Delta angle, for d2fr: theta bin width is 2*dang (ang +/- dang)
					! for dofr: width of theta bin is dang

	common /spefra/  d2fr(30,11,152), dofr(30,37), yfr(2,30),
     &		ang(10), dang		!LMK, added ang and dang to spefra	06/2012

	data zro /0.d0/

c ========================

      read (15, *) ang, dang
c 
        do kf=1,30
        yfr(1,kf) = zro
        yfr(2,kf) = zro
          do ka=1,37
          dofr(kf,ka) = zro
          if(ka <= 11) then
            do kt=1,152
            d2fr(kf,ka,kt) = zro
            enddo
          endif
          enddo
        enddo
      return
	end
c 


      subroutine PisaPrint (sigin, ncas)
c LMK, Normalizes d2fr, dofr and finishes double differential cross section and angular
c distribution calculations, and then prints them.

c LMK, Variable descriptions
	integer*8 ncas		! Number of inelastic events simulated
	integer*4 ia		! Counting variable for the 10 theta bins in d2fr
	integer*4 ka		! Counting variable for the 37 (=180/(dtet=5) + 1) theta bins
	integer*4 kfr		! Particle fragment number (1-30) of particle m
	integer*4 it		! Counting variable for the 150 Tk bins
	integer*4 jgr		! Gives the index of the first bin in the nth interval
	real*8 sigin		! Inelastic cross section (mb) (Total, across all events?)
	real*8 fn		! # of inelastic events in DOUBLE
	real*8 fac		! = Inelastic cross section / Number of inelastic events
	real*8 ang1, ang2	! Lower and upper angles used to find dom and domp
	real*8 dom		! Delta omega (sr) for the 36 theta bins in dofr calculations
	real*8 domp		! Delta omega (sr) for the 10 theta bins in d2fr calculations
	real*8 c		! Double differential cross section of the ia-th theta bin and
					! it-th Tk bin (mb/MeV/sr)
	real*8 sfr1		! Total cross section for all energies (mb) (= yfr(1,kfr)*fac)
			! Cross section of fragment k = (N_k)*(Sigma_inelastic)/(N_inelastic)
	real*8 sfr2		! Total cross section for Tmin<=Tk<=Tmax (mb) (= yfr(2,kfr)*fac)
	real*8 dtk		! Bin width of the it-th Tk bin
	real*8 tk1, tk2		! Lower and upper energies of the it-th Tk bin
	real*8 summ		! Used to see if there were any particle tallies in the Tk bin, and
					! therefore whether or not to print x-section info
	real*8 d2fr		! Bin count of all events in specific theta and Tk bins, ie 
					! unnormalized Double Differential cross section
	real*8 dofr		! Bin count of all events in a given theta bin, irrespective of Tk
					! ie unnormalized Angular Distribution
	real*8 yfr		! Number of particles tallied: 
					! yfr(1,kfr) = total # of particles tallied, 
					! yfr(2,kfr) = # of particles with Tmin<=Tk<=Tmax
	real*8 tmin, tmax	! Min and max kinetic energy to consider, for each of the 
					! 30 different particles (neutrons through C14)
	real*8 tgr		! Defines break points in kinetic energy, for bins
	real*8 dtgr		! Defines size of Tk bins between tgr(n-1) and tgr(n)
	real*8 ang		! Ten user input angles to calculate d2fr at
	real*8 dang		! Delta angle. For d2fr: theta bin width is 2*dang (ang +/- dang)
					! For dofr: width of theta bin is dang
	real*8 dtet		! Delta theta for bins, set as a constant
	real*8 degrad		! Converts from radians to degrees (I think)
	real*8 pi, twpi		! Constants Pi and 2*Pi
	character*4 fname(30)

      common /spefra/  d2fr(30,11,152), dofr(30,37), yfr(2,30),
     &		ang(10), dang		!LMK, added ang and dang to spefra	06/2012
      common /degrad/  degrad
      common /pi/      pi, twpi

      dimension tgr(4), dtgr(4), jgr(4)
      dimension dom(36), domp(10),
     &          tmin(30), tmax(30), c(11)
      data fname/'   n','   p','   d','   t',' He3',' He4',
     &           ' He6',' Li6',' Li7',' Li8',' Li9',' Be7',
     &           ' Be9','Be10','  B9',' B10',' B11',' B12',
     &           ' C11',' C12',' C13',' C14',' Z=7',' Z=8',
     &           ' Z=9','Z=10','Z=11','Z=12','Z=13','Z=14'/
      data tmin/  0.0d0, 2.0d0, 2.6d0, 3.0d0, 2.0d0, 2.5d0, 
     &            2.5d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0, 4.5d0,
     &            4.5d0, 4.5d0, 9.0d0, 9.0d0, 9.0d0, 9.0d0,
     &           11.0d0,11.0d0,11.0d0,11.0d0,14.0d0,16.0d0,
     &           16.0d0,16.0d0,16.0d0,16.0d0,16.0d0,16.0d0/ 
      data tmax/ 160.d0,160.d0,215.d0,250.d0,580.d0,650.d0, 
     &            2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,
     &            2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,
     &            2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3,
     &            2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3, 2.5d3/
      data dtet /5.0d0/  
      data  tgr/50.0d0, 350.0d0, 950.0d0, 2600.0d0/,
     &     dtgr/ 2.0d0,   5.0d0,   20.0d0,   50.0d0/,
     &      jgr/     1,      26,       86,      116/     
      data zro, one, two, thr, for, mil /0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 
     &                             1.d3/

c ======================================================================

      fn = dble(ncas)
      fac = sigin/fn 
c LMK, Calculates delta omega: domp for double differential and dom for angular distributions
	! delta omega = 2*Pi*[cos(theta1) - cos(theta2)]
        do ka=1,36
c LMK, calculates domp for the 10 theta bins in double differential cross section calculations
        if (ka <= 10) then
          ang1 = (ang(ka) - dang)/degrad
          if(ang1 < zro) ang1 = zro
          ang2 = (ang(ka) + dang)/degrad
          if(ang2 > pi) ang2 = pi
          domp(ka) = twpi*(cos(ang1) - cos(ang2))
        endif
c LMK, calculates dom for the 36 theta bins in angular distribution calculations
        ang1 = dble(ka-1)*dtet/degrad         
        ang2 = dble(ka)*dtet/degrad         
        dom(ka) = twpi*(cos(ang1) - cos(ang2))
        enddo

c LMK, Main double differential cross section loop.....finishes calculations and prints
c LMK, Loop over 30 different fragment types
        do kfr=1,30
        if(yfr(1,kfr) <= zro) go to 50	! LMK, If there are no particles tallied for fragment type 
					! kfr, skip to next fragment type.
c LMK, First calculate and print total cross section
        sfr1 = yfr(1,kfr)*fac		! LMK, Total cross section for fragment kfr
        sfr2 = yfr(2,kfr)*fac		! LMK, Total x section for fragment kfr in Tmin<=Tk<=Tmax
        write (31, 101) fname(kfr), sfr1, tmin(kfr), tmax(kfr),
     &                  sfr2, ang
  101 format(//25x,'Double differential cross-section d2S/dTdO ',
     &       '(mb/MeV/sr) of ',a4/12x,
     &       'prod. xsec for all energies=',1pe11.4,' mb, for ',
     &        0pf5.0,'<T<',f5.0,' xsec=',1pe11.4,'mb'/2x,
     &       'T(MeV)/angle:',9(0pf5.0,6x),f5.0,4x,'dS/dT(mb/MeV)')
c LMK, Loop over 150 Tk bins
          do it=1,150
          summ = zro 		!LMK, used to see if there are any particles to print x section for
c LMK, First determine the bin width and lower/upper energies for each Tk bin
c LMK, Check to see which of the four intervals the it-th Tk bin is in, and assign dtk and tk1
          if (it >= jgr(1) .and. it < jgr(2))            then
            dtk = dtgr(1)		! LMK, Bin width
            tk1 = dble(it-jgr(1))*dtk	! LMK, Lower energy
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
          endif
          tk2 = tk1 + dtk		! LMK, Upper energy
c LMK, Loop over 10 theta bins
c LMK, Calculate double diff cross section for each of the ten theta bins, using d2fr tallies
            do ia=1,10
            c(ia) = d2fr(kfr,ia,it)*fac/(domp(ia)*dtk)	!LMK, (mb/MeV/sr)
	! LMK, c(ia) = Double diff x section at ia-th theta bin, it-th Tk bin, for kfr fragment
            summ = summ + c(ia)
            enddo
          c(11) = d2fr(kfr,11,it)*fac/dtk	!LMK, Double diff x section integrated over angle
						! (mb/MeV)
          summ = summ + c(11)
          if(summ > zro) write (31, 102) tk1, tk2, c
  102 format(f5.0,'-',f5.0,11(1pe11.4))
          enddo
c LMK, Calculate cross sections for all particles Tmin<=Tk<=Tmax in the 10 theta bins
          do ia=1,10
          c(ia) = d2fr(kfr,ia,152)*fac/domp(ia)	!LMK, DD x section integrated over energy (mb/sr)
          enddo
        c(11) = d2fr(kfr,11,152)*fac		!LMK, Total cross section (mb)
        write (31, 103) c
  103 format(' energ. int.',11(1pe11.4))
   50   continue
        enddo 
c LMK, Done with double differential cross sections

c LMK, Now calculate and print angular distributions for the first 10 fragment types (through 8Li)
      write (31, 104) (tmin(kfr), tmax(kfr),kfr=1,10),
     &                (fname(kfr),kfr=1,10)
  104 format(//10x,'Angular distribution of produced fragments',
     &  ' dS/dOm [mb/sr] for energy range(MeV)'/'Tmin-Tmax',
     & 10(f5.0,'-',f5.0)/ 'ang1-ang2',10(3x,a4,4x)) 
c LMK, Loop over 36 theta bins
        do ka=1,36
        ang1 = dble(ka-1)*dtet         !LMK, Calculate lower and upper angles
        ang2 = dble(ka)*dtet         
          do kfr=1,10
          c(kfr) = dofr(kfr,ka)*fac/dom(ka) 	!LMK, Angular distribution of fragment kfr, 
						! ka-th theta bin (mb/sr)
          enddo
        write (31, 105) ang1, ang2, (c(kfr),kfr=1,10)
  105 format(f4.0,'-',f4.0,10(1pe11.4)) 
        enddo
c LMK, Print total cross section for fragment kfr
        do kfr=1,10
        c(kfr) = dofr(kfr,37)*fac
        enddo 
      write (31, 106) (c(kfr),kfr=1,10)
  106 format(' Int. x sec',10(1pe11.4))    
c LMK, Now calculate and print angular distributions for the 10-20 fragment types (9Li-12C)      
      write (31, 107) (tmin(kfr), tmax(kfr),kfr=11,20),
     &                (fname(kfr),kfr=11,20)
  107 format(//10x,'Angular distribution of produced fragments',
     &  ' dS/dOm [mb/sr] for energy range(MeV)'/'Tmin-Tmax',
     & 10(f5.0,'-',f5.0)/'ang1-ang2',10(3x,a4,4x)) 
c LMK, Loop over 36 theta bins
        do ka=1,36
        ang1 = dble(ka-1)*dtet         !LMK, Calculate lower and upper angles
        ang2 = dble(ka)*dtet         
          do kfr=11,20
          c(kfr-10) = dofr(kfr,ka)*fac/dom(ka) 	!LMK, Angular distribution of fragment kfr, 
						! ka-th theta bin (mb/sr)
          enddo
        write(31,108) ang1,ang2,(c(kfr),kfr=1,10)
  108 format(f4.0,'-',f4.0,10(1pe11.4)) 
        enddo
c LMK, Print total cross section for fragment kfr
        do kfr=11,20
        c(kfr-10) = dofr(kfr,37)*fac
        enddo 
      write (31, 109) (c(kfr),kfr=1,10)
  109 format(' Int. xsec',10(1pe11.4))          
c LMK, Now calculate and print angular distributions for the last 10 fragment types (13C-Z14/28Si)
      write (31, 110) (tmin(kfr), tmax(kfr),kfr=21,30),
     &                (fname(kfr),kfr=21,30)
  110 format(//10x,'Angular distribution of produced fragments',
     &  ' dS/dOm [mb/sr] for energy range(MeV)'/'Tmin-Tmax',
     & 10(f5.0,'-',f5.0)/'ang1-ang2',10(3x,a4,4x)) 
c LMK, Loop over 36 theta bins
        do ka=1,36
        ang1 = dble(ka-1)*dtet         !LMK, Calculate lower and upper angles
        ang2 = dble(ka)*dtet        
          do kfr=21,30
          c(kfr-20) = dofr(kfr,ka)*fac/dom(ka) 	!LMK, Angular distribution of fragment kfr, 
						! ka-th theta bin (mb/sr)
          enddo
        write (31, 111) ang1, ang2, (c(kfr),kfr=1,10)
  111 format(f4.0,'-',f4.0,10(1pe11.4)) 
        enddo
c LMK, Print total cross section for fragment kfr
        do kfr=21,30
        c(kfr-20) = dofr(kfr,37)*fac
        enddo 
      write (31, 112) (c(kfr),kfr=1,10)
  112 format(' Int. xsec',10(1pe11.4))
      return

c ======================================================================
      end

