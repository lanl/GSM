
  subroutine gsmMain( gsmObj, projNucleus, targNucleus, output, maxProgeny )

! ======================================================================
!
!                              NOTICE
!
!    This software and ancillary information (herein called "SOFTWARE")
!    named CEM03.03 is made available under the terms described here.
!    The SOFTWARE has been approved for release with associated LA-CC
!    number LA-CC-04-085.
!
!    Copyright (2012). Los Alamos National Security, LLC.
!
!    This material was produced under U.S. Government contract
!    DE-AC52-06NA25396 for Los Alamos National Laboratory, which is
!    operated by Los Alamos National Security, LLC, for the U.S.
!    Department of Energy. The Government is granted for itself and
!    others acting on its behalf a paid-up, nonexclusive, irrevocable
!    worldwide license in this material to reproduce, prepare derivative
!    works, and perform publicly and display publicly.
!
!    NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF
!    ENERGY, NOR LOS ALAMOS NATIONAL SECURITY LLC, NOR ANY OF THEIR
!    EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
!    LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
!    OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS
!    DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY
!    OWNED RIGHTS.
!
!    Additionally, this program is free software; you can redistribute
!    it and/or modify it under the terms of the GNU General Publi!
!    License as published by the Free Software Foundation; either
!    version 2 of the License, or (at your option) any later version.
!    Accordingly, this program is distributed in the hope that it will
!    be useful, but WITHOUT ANY WARRANTY; without even the implied
!    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!    See the GNU General Public License for more details.
!
! ======================================================================
!
!    The primary authors of CEM03.03 are:
!    S. G. Mashnik (LANL), K. K. Gudima (IAP), and A. J. Sierk (LANL);
!    with important contributions from R. E. Prael (LANL), M. I. Baznat (IAP),
!    and N. Mokhov (Fermilab).  (IAP = IAP, Academy of Science of Moldova.)
!
! ======================================================================
!
!    Final version for CEM03.03 release, AJS, LANL T-2, February, 2012.
!
!    INCOMPLETE LIST OF CHANGES:
!    1) Restricted asymmetric fission fragments to A >= 13 to prevent
!       a very rare problem.
!    2) Added Fermi breakup after each stage of the reaction to decay
!       any excited nuclei with A < 13 produced at the stage considered.
!    3) Replaced the RNDM generator with the RANG generator, which is
!       about 5 times faster because it can be compiled inline, and
!       has a period of about 2^63 ~ 10^19. See comments in the module
!       ranmc.f.
!    4) An early error derived wrong values for the a_f/a_n
!       multipliers found in the routine FITAFPA for the nucleus
!       181 Ta, which would affect fission cross sections for Z = 72 and
!       73. New values which now actually do reproduce the systematics
!       of Prokofiev are included in the routine.
!    5) Some fixes to the Fermi breakup were made by K. K. Gudima,
!       R. Prael of LANL, and S. G. Mashnik.
!    6) Integers are now explicitly typed as Real*4 or real*8 as needed.
!    7) Several other bug fixes.
!    8) Some additional clarifying comments were added to various SRs.
!
! ======================================================================
!    Various minor mods. were made; 2005-2009
! ======================================================================
!
!   October 7, 2005.
!
!   Incorporated several fixes to problems cropping up when Mokhov
!   runs a CEM03 event generator inside MARS in an MPP implementation of
!   MARS.  Also, some revised values of af/an for fission of Tungsten
!   isotopes were added to reflect recent data. AJS
!
! ======================================================================
!
!   February 14, 2005:
!
!   There have been several major changes introduced by Gudima and
!   Baznat, including redetermining af/an values for fission.
!     **** See 4) above ****                                  Mashnik
!   introduced Kalbach systematics for preequilibrium decay, some
!   improved angular distribution treatments for N-N and gamma-N
!   reactions and incorporation of systematics for gamma-A total cross
!   sections were implemented by Gudima and Mashnik. Sierk has
!   rewritten much of this later coding to vastly improve the execution
!   speed of the code, and also made some speedups in GEM2. Many
!   additional output options have been implemented. The important
!   improvements are discussed in the CEM03.01 Manual.
!
! ======================================================================
!
!   March 11, 2004.
!
!   THE CODE SHOULD BE COMPILED TO USE 8-BYTE INTEGERS!!!!!
!    ****  See 6) above ****
!   (Large simulations exhaust the range of 32-bit integers.)
!
!   Improvements made by A. J. Sierk, LANL T-16:
!
!   Edited various modules to have a consistent style (lower case,
!   indentation, etc.).
!   Corrected a number of logic errors and outright mistakes.
!   Fixed the photonuclear part of the cascade, which was incorrectly
!   implemented previously.
!   Removed calculation of several arrays from inside SIMULATEDECAY and put
!   the quantities in larger arrays to save time (ALJ and GB).
!   Broke up SIMULATEDECAY into a number of smaller and more manageable
!   pieces; (new) SIMULATEDECAY, PREQAUX, AUXL, PEQEMT, EQDECY.
!   Removed a large part of TYPEOUT into 2 separate subroutines PRTDIST
!   and PRTMULT to remove redundant code fragments.
!   Put in GAUSS2, a much more efficient gaussian-distributed random
!   number generator in place of GAUSSN, used by GEM2.
!   Enhanced the breakdown of emitted particles into those from
!   spallation, prefission, fission fragments, coalescence, etc.
!   Fixed the GSMOBJ%PRINTPARTPROP routine so it works properly for debugging the
!   preequilibrium, evaporation, and fission fragment production.
!   Converted many (all, I hope, except large tables in MOLLNIX, etc.)
!   constants to double precision.
!   Used the array a3rd to remove many calculations of A^1/3.
!   Put the reduced mass calculation into INITIAL to save time, also
!   ignoring the mass excess to remove the slight Z-dependence of
!   the reduced masses. **** THIS IS NOT IMPLEMENTED IN CEM03.03 ****
!   Fixed the calculation of statistics of the nuclei at various
!   stages of the collision. (See TYPEOUT, PRTMULT, PRTDIST)
!   Completely removed the iopt = 1 option, which gives the CEM95
!   version of the code.  This is obsolete and unuseable with the
!   modified GEM2 code modules.
!   Removed and/or commented some portions of unused code.
!   Replaced many INT statements with NINT where appropriate.
!   Removed all reference to file 14, the old cemxx.out file.
!   Replaced the Wapstra + Cameron calculations in the ENERGY function
!   by Wapstra and Moller-Nix, where MN are defined; use MNMACRO for
!   heavier nuclei outside the MN table, and Cameron (based on M(O16)
!   = 16.000!!) only for EXTREMELY neutron-rich isotopes of Z = 1-7.
!   Removed the request for the name of the input file; this is always
!   cem03.inp in the current version.
!   Switched to REAL*8 typing for more precision in variable typing.
!   Renamed the auxiliary ....tbl files to lower case.
!   Updated the constants in md in GITAB to be current with a newer
!   version of the Wapstra-Audi mass table.
!   Removed all EQUIVALENCE statements as a token move toward f90/f95.
!   Fixed the coalescence module; the trajectory separation was not
!   defined, so was always zero and the comparison to a radius
!   parameter was always satisfied. This leads to MUCH less coalescence.
!
! ======================================================================

!*******************************************************************
! Cascade-Exciton Model by Stepan Mashnik, A. J. Sierk,
! Konstantin Gudima et al.
! Version CEM2KPH of July 2003
!
! cem03.f, fermi2k2.f, gadd90.f, gem2fit.f, photo.f:
! modified August 2003 by NVMokhov
!
!   "Last" change: 14-AUG-2003 by NVMokhov
!
! Modifications by NVMokhov: 12-14 August 2003:

! Conversion to Global DOUBLE PRECISION

! ALOG  replaced with LOG
! Similarly: DABS->ABS, DSQRT ->SQRT, DCOS->COS, DSIN->SIN, DEXP->EXP
! AMAX1 replaced with MAX
! AMIN1 replaced with MIN
! IFIX  replaced with INT
! FLOAT replaced with DBLE
! REAL  replaced with DBLE
! DNINT replaced with ANINT
! IDNINT replaced with NINT

! functions RAN1, RAN3 and RNDM(-1) removed,
! double precision MARS'S RNDM(-1.) is used instead
!   ****  See 3) above ****

! Modifications by NVMokhov: 9-15 June1998, update 08/14/03:

! subroutine TRANS  renamed to TRANS8
! subroutine DIRECT renamed to DIRECT8
! function   COLOMB renamed to COLOMB8
! function   SIGMAT renamed to SIGMAT8
! function   GEOM   renamed to GEOM8

! common/O/      has been aligned
! common/VUL/    has been aligned
! common/ADBF/   has been aligned
! common/MENU2/  has been aligned
! common/FISCEM/ has been aligned

! Modifications by NVMokhov: 01-DEC-1998, update 08/14/03:

! function FIS   -> subroutine sfis
! function FINT  -> subroutine sfint
! function FINT1 -> subroutine sfint1
! function FINT2 -> subroutine sfint2 - commented (unused)
! function fname -> subroutine sfname
! function s1    -> subroutine ss1
! function S2    -> subroutine ss2

! Modification by NVMokhov: 09-Oct-1999, update 08/14/03:
! common /INDEX/ -> /BINDEX/
! common /FISS/ replaced with /FISCEM/

!   "Last" change: 14-AUG-2003 by NVMokhov
!   "Last" change: 17-SEP-2004 by KKGudima
!
! ======================================================================
!
!   Changes made by S. G. Mashnik, LANL T-16, 2000-2001 to get a
!        preliminary version of CEM2k:
!
!    1. Changed the method for transition from the CASCADE stage of a
!       reaction to the preequilibrium and the time of equilibration.
!       For incident energies above 150 Mev, the transition from the
!       cascade stage of a reaction to the preequiliblium is determined
!       comparing the energy of the cascade nucleons to a cutoff energy
!       of 1 MeV (still a parameter): all cascade nucleons with energies
!       less than 1 MeV above the Fermi surface are absorbed (as excitons),
!       to continue the following equilibration of the preequilibrium nucleus
!       and emit particles at the preequilibrium stage; the cascade nucleons
!       with higher energies are considered as cascade particles with
!       further interactions with intranuclear nucleons (spectators), if the
!       Pauli principle allows this.
!
!    2. At incident energies above 150 Mev, only transitions with Delta_n = +2
!       (in the direction of equilibration) are allowed at the preequilibrium
!       part of a reaction. This makes the preequilibrium stage shorter
!       and the evaporation stage longer, that in couple with 1) result in
!       a longer evaporation part of a reaction, with higher excitation
!       energy, and, as a result, more evaporated neutrons. Both points
!       1) and 2) were inspired by the recent GSI measurements (Phys. Rev.
!       Lett. 84 (2000) 5736): Incorporating 1) and 2) in CEM97 allows us
!       to describe well the GSI data, while without these modifications
!       we got with CEM97 too few neutrons and a bad agreement with the
!       data, just as we got with other 7 different codes we tested.
!
!    3. To be able to describe correctly backward emitted particles at
!       incident energies bellow 150, MeV we keep in CEM2k the "Proximity"
!       parameter, P=0.3, used previously in CEM97 and CEM95, and take into
!       account all transitions, Delta_n = +2, 0, and -2 at the
!       preequilibrium stage of a reaction. Some "physical" speculations
!       may be made about this (changing the time of the cascade and
!       preequilibrium stages of reactions with increasing the energy of the
!       bombarding particle), but we do not have now a good understanding of
!       how quantitatively and qualitatively different mechanisms of nuclear
!       reactions change with the bombarding energy, and the mentioned
!       changes here are of a more pragmati! type: to describe better all
!       the available data.
!
!    4. The real binding energies were incorporated for nucleons at the
!       cascade stage of reactions instead of the approximation of a constant
!       separation energy of 7 MeV used in previous versions of the CEM. This
!       required an additional routine, BindNUC, and improved significantly
!       the agreement with the data for some excitation functions, like
!       (p,n), (n,p), (p,2p), etc., where the threshold of reactions is
!       important and should be calculated accurately.
!
!    5. The momentum-energy conservation for the each cascade simulated
!       event was imposed using real masses of the projectile and target
!       and of emitted cascade particles and residual excited nucleus after
!       the cascade: The excitation energy of a residual nucleus is
!       redefined according to the momentum-energy conservation law, or
!       this cascade event is considered false and resimulated, if the
!       conservation law does not provide a positive excitation energy for
!       the simulated already emitted particles and residual nucleus.
!       A new routine, RENORM, was written for this, and included in CEM2k.
!       The last points 4) and 5) have solved the problem #1 in the list
!       of "deficients" written above in 1997 by A.J. Sierk, and allow to
!       get a better agreement with the data. Note, than the "deficient #6
!       of the mentioned list was solved in 1999 by A.J. Sierk: the present
!       version of the code is several times (up to 6-7, for heavy targets)
!       faster than the initial version of CEM97.
!
!    6. The reduced masses of particles were included in calculation particle
!       widths at preequilibrium and evaporation stages of a reaction instead
!       of using the approximation of no-recoil used in previous versions of
!       CEM. This point is especially important when calculating emission of
!       complex particles from light nuclei, and allowed to get a much better
!       agreement with the data for the yields of He4, He3, t, and d (gas
!       production) from light nuclei. But to solve completely the problem
!       of complex particle emission (gas production) from arbitrary targets,
!       better inverse cross sections have to be incorporated in CEM2k, as
!       mentioned above in "deficient" #2. We derived already our own
!       universal systematics for inverse cross sections of complex particles
!       and light fragments and plan to incorporate them in CEM2k in FY2002,
!       solving the "deficients" #2 and 4 from the above list.
!
!    7. To solve the problem of using reliable total reaction cross sections
!       (for the absolute normalization of calculated yields) at incident
!       energies bellow 100 Mev, we have incorporated into CEM2k the NASA
!       systematics by Tripathi et al. (niM b117 (1996) 347) for all incident
!       protons and neutrons with energies above the maximum in the reaction
!       cross section, and the Kalbach systematics (J. Phys. G: Nucl. Phys.,
!       24 (1998) 847) for neutrons of lower energies. All nucleon-induced
!       reaction results by CEM2k are now normalized to these reaction cross
!       sections (printed in the output file as "Used here inelastic cross
!       section"), while the old inelastic cross section calculated by the code
!       using the geometrical cross section and the numbers of inelastic and
!       elastic simulated events are printed nearby in the output file as
!       ("Monte Carlo inelastic cross section", so a renormalization of results
!       to the old cross section is possible as well. Note, that the new
!       normalization provides a much better agreement with practically all
!       available nucleon-induced reaction data.
!
!  with KKG and MIB (Baznat) changes for interpolations of the afMultiplier
!  and czMultiplier
!
!  last changes(gb, isotope printing,ang. momenta etc.) by KKG 07/12/02
!  sigpre is included by KKG  07/15/02
!  Fermi decay of excited nuclei (A <= 12) is included by KKG 08/06/02
!  last corrections from 10.08/02 by KKG
!  Coalescence of the cascade nucleons into d,t,he3,he4 is included by
!  KKG 02/06/2003
!  last corrections by KKG&MIB 03/20/03 (see below and vlob) and
!  23/06/03 (see cascad)
!
!   "Last" change: 13-AUG-2003 by NVMokhov
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by K. K. Gudima, LANL T-16, Feb.-Mar., 2004.
!
! ======================================================================
!
!*******************************************************************
!
!   December, 1997
!
!   This program is a modified version of the code CEM95 written by
!   S. G. Mashnik.
!
!   Many of the historical options have been removed; so that only the
!   current best physics options are retained.  There are two basic
!   options:
!     1).  The CEM95 program as released by Mashnik in 1995; with its
!          preferred options for level densities (iljinov, et al.,#3,
!          1992), microscopic ground state corrections (Truran, Cameron,
!          and Hilf, 1970), and macroscopic fission barrier heights
!          (krappe, Nix, & Sierk, 1979). This version of the code
!          uses old-style nuclear masses (M(O-16) = 16.0000),
!          and the Cameron (1957) macroscopic energy calculation; it
!          also assumes all nucleon masses are 0.94 Gev, and all pion
!          masses are 0.14 GeV.
!          This version uses the new elementary N-N and pi-N
!          cross sections in the cascade, so is not identical to the
!          original CEM95; differences are small in the cases so far
!          checked.
!
!     2).  The partially upgraded 1997 version which includes the
!          mass tables of Moller, Nix, Myers, & Swiatecki, with
!          the corresponding shell corrections and pairing gaps,
!          the level densities as calculated by Sierk in 1996 to be
!          consistent with the Moller-Nix microscopic corrections using
!          identical methods to those of Iljinov, et al., new fits
!          to the elementary N-N and pi-N cross sections used in the
!          cascade portion of the model, actual free nucleon and pion
!          masses used in the cascade, and the Sierk (1986) angular
!          momentum dependnt macroscopic fission barriers, which are
!          consistent with the macroscopic model used in the Moller-
!          Nix mass calculations.
!
!    NOTE: This version still has many deficiencies which will be
!          addressed in the coming months.  This preliminary release
!          is being made because of the desire of the MCNPX team to
!          have a version of the code available.
!
!    Among the deficiencies are:
!
!     1).  The cascade does not conserve energy properly. The energy
!          of the Fermi surface is assumed to be -7.0 Mev, irrespective
!          of the actual nuclide. The excitation energy of the
!          residual nucleus after emission of a cascade particle is
!          not calculated with the actual ground state energy or
!          particle separation energy.
!
!     2).  The inverse cross sections used in the calculation of
!          preequilibrium and compound particle decay widths are not
!          a good representation of real cross sections.
!
!     3).  While calculating fission, the code does not have a model
!          for selecting fission fragments and following their
!          subsequent decay. When fission is a significant channel,
!          the emission of low energy nucleons will not be properly
!          modeled.
!
!     4).  The code does not have a model for fragmentation (emission
!          of particles heavier than alphas), either in the pre-
!          equilibrium or the compound decay parts of the code.
!
!     5).  Nuclei near the proton or neutron drip line are allowed to
!          to exist, although their particle decay will occur in much
!          less than the compound decay times implicitly assumed in
!          the model.
!
!     6).  No serious attempt has been made to improve the execution
!          speed of the code by analyzing where the majority of time
!          is being spent.
!
!     ....
!
! ======================================================================
!
! NOTE:   In the derivation of the level density parameters by Iljinov,
!         et al., pairing of 11.0/sqrt(A) was used, although the value
!         of 12.0 was stated in the paper.  This code uses 12.0, which
!         is inconsistent with the parameter sets.
!
!  Changes made at LANL, July-August 1995 (D. W. Muir):
!
!     1.  3 redundant variables removed (prinp:zapran,franz and partn:v)
!     2.  hollerith constants replaced with single-quote limited strings
!     3.  program card activated
!     4.  all Uc alphabetic characters outside of strings lowered.
!     5.  comments added to functions ranf1 and rndm to explain purpose
!     6.  save statement added to function ranf1 (not actually required)
!     7.  added cpu timing call and timing print
!     8.  most variables read from a single record for easier updating
!
!   Changes made by A. J. Sierk,  LANL  T-2, February-May, 1996; July-
!        December, 1997:
!
!    Moved program card to beginning of file.
!    Changed arithmetic if and go to (a,b,c...) statements to
!        if..then..else constructions; thereby removing many numbers
!        from statements.
!    Added indentation to if..then..else  and  do loops.
!    Inserted real and int functions to remove implicit typing where
!        noticed.
!    Explicitly open unit 15; cem97n.inp .
!    Change unit 16 to cem97n.info; for error and diagnostic messages.
!    Changed comment lines to be more easily readable.
!    Added extensive new comments.
!    Changed output formats to enhance readability.
!    Suppressed printout of descriptions of calculation after the
!        first energy; originally this was repeated for each energy
!        of the projectile.
!    Reduced maximum length of lines to 72 columns.
!    Arranged common blocks alphabetically in each subroutine.
!    Converted to generic Fortran implicit functions (alog -> log, etc.)
!   NOTE:  This RNDM was a different generator and has been discarded;
!   AJS  1/07/05.
!    Modified rndm to remove a problem with random numbers [rndm(i-1)]
!        greater than 0.999 . [The succeeding random number [rndm(i)]
!        was constrained to be in the neighborhood of 0.80, the spread
!        depending on the magnitude of 1.0 - rndm(i-1)].
!    Removed the ground-state shell correction calculation from bf,
!        putting it into a call to shell.
!    Changed the calling arguments to shell, to implement above^.
!    Added the parameter ipair, which allows selection between the
!        original approximate pairing gap of 12.0/sqrt(A) * chi and
!        the tabulated pairing gaps of Moller, Nix & Kratz (1997);
!        [Atomic Data Nucl. data Tables 66, #2 131-343].
!    Added the character variable heading(nh) to allow for printing a
!        descriptive text at the beginning of cem97n.res.
!    Added the Moller-Nix shell correction option to the ground-state
!        microscopic correction to the barrier height.
!    Put in tabulated (where available) or calculated (Moller-Nix)
!        mass excesses, in lieu of the formula in deltam.
!        For nuclei outside the table, use the macroscopic RLDM
!        formula from Moller, Nix, Myers, & Swiatecki.
!    Derived additional semiempirical level density parameters using
!        Moller-Nix ground-state microscopic corrections, both with
!        and without Moller-Nix-Kratz pairing gaps.
!        [Atomic Data Nucl. data Tables, 59, 185 (1995)]
!        (DISCOVERD THAT ILJINOV, ET AL. USED 11.0/sqrt(A) in DERIVING
!        THEIR LEVEL DENSITY SYSTEMATICS INSTEAD OF THE STATED VALUE
!        OF 12.0/sqrt(A)!!)
!    Corrected the energetics in SIMULATEDECAY to use a more accurate mass
!        for the compound or residual nucleus.
!    Removed the use of the ground-state microscopic correction from
!        the calculation of the level density parameter at the
!        saddle point.
!    Extended the number of channels printed out for mchy=1.
!    Added the mpyld flag to print out a summary table of particles
!        (n, p, d, t, He-3, alpha, pi +,-,0) emitted without the detail
!        of mchy=1.
!    Print out the time, the energy, and the reaction at the end of
!        each energy of a multi-energy run.
!    Made several changes to eliminate infinite loops.
!
! ======================================================================
!
!                    Program input variables (2016 version and older):
!**cm0     = projectile mass in GeV (enter a negative number to end run)
!  projNucleus%kinEnergy   = initial projectile kinetic energy (in GeV)
!  anucl   = target mass number.
!  znucl   = target proton number.
!**me0     = projectile charge (in units of proton charge)
!**mb0     = projectile baryon number
!  Above projectile quantities cm0, me0, mb0, are OBSOLETE!
!  limc    = total number of inelastic events, normally 2000-500000.
!             => Migrated to GSMOutput%numInelasticEvents
!**nnnp    = number of initial cascade steps with detailed diagnostics.
!  nnnp REMOVED!
!**wam     = ratio of saddle-point fission level density to compound
!            level density; af/ac.  (Usually called af/an).
!**idel    = do evaporation calculations with or without fission.
!          = 0 do calculations without fission.
!          = 1 allow competition between fission and particle emission
!            at compound stage.
!          = 2 to accumulate particle spectra, etc. ONLY for fission
!            events.
! NOTE:  idel removed from input; idel=2 function moved to mpyld!
!  dt0     = projectile kinetic energy step size in MeV.
!  t0max   = final projectile kinetic energy in MeV (last energy run
!            will be less than or equal to t0max).
!     Note:  To run a single energy, input a negative dt0 (-5.0) and
!             a value of t0max larger than 1000.*projNucleus%kinEnergy
!  dteta   = step size (degrees) in ejectile angular distributions.
!  mspec   = (0/1/2) if ejectile energy spectra (not/are) needed.
!             (See variables ipar1 and ipar2.)
!            mspec = 2 gives a breakdown of the source of evaporated
!            particles (prefiss, fiss frags, evap residue).
!  mpyld   = (0/1/2) if particle yield summary table (is not/is)
!            desired.  Prints cross sections and mean energies for
!            particles n -- pi+ broken into total, cascade,
!            preequilibrium, and compound production.
!            mpyld = 2 gives multiplicities ONLY for fission events.
!  mchy    = (0/1) if excitation functions (are not/are) needed.
!            If mchy = 1, 192 "channel" cross sections  with values of
!            at least 1 mb and up to 26 summed channel cross sections
!            in mb  (not isotope yields, see misy) are calculated.
!            To save cpu time in routine calculations, use mchy=0.
!  misy    = (0/1 or 2, or 3) If isotope yields (are not/are) needed.
!            If misy=1, the primary yield for all nuclides with non-zero
!            production will be printed.
!            If misy=2, the primary yield for all nuclides with non-zero
!            production and forward/backward correlations will be printed.
!            If misy=3, the residual (after cascade and preeq.) nuclear
!            information, and fission-fragment opening angle distribution,
!            and the neutron multiplicity distribution will be printed.
!  mdubl   = (0/1/2) if ddx spectra (are not/are) needed.
!            (See variables ipar1 and ipar2.)
!            mdubl = 2 gives a breakdown of the source of evaporated
!            particles (prefiss, fiss frags, evap residue).
!  mang    = (0/1/2) if angular distributions (are not/are) needed.
!            (See variables ipar1 and ipar2.)
!            mang = 2 gives a breakdown of the source of evaporated
!            particles (prefiss, fiss frags, evap residue).
!  ipar1,ipar2 = range of ejectile types for spectrum calcs.
!            See variables mspec, mang, mdubl.
!            Types are defined as follows:
!            1 n          4 t          7 pi-
!            2 p          5 He3        8 pi0
!            3 d          6 He4        9 pi+
!
!  tet1(j),tet2(j) = up to 10 angle bins for angular distributions.
!            used if mdubl=1.  Terminate with a negative tet1 value.
!  tmin(j),tmax(j),dt(j) = 4 energy regions and bin sizes (MeV).
!            used if mdubl>=1 or mspec>=1
!
!**The following 5 flags are set in INITIAL!!
!**noprec  = (0/1) for (using/not using) preequilibrium decay.
!**sigpre  = value of "smearing" parameter for transition from
!            preq. to compound decay. Usually 0.4 for CEM03.
!**nevtype = number of particle types evaporated in GEM2; 6 for
!            p-->alpha; 66 for all isotopes up to Mg-28. See table
!            in Block Data BDJEC.
!**irefrac = (0/1) for (NOT/changing directions) when passing between
!            zones of the potential during the cascade.
!**icostan = (0/1) for using (OLD/NEW) approximations for the
!            angular distributions of N-N elastic scattering and
!            gamma-N [Duarte approx for T <= 2 GeV; extrapolation
!            to avoid unphysical deviations in old approximations for
!            angles near 0 and pi for other energies.
!
!  nh      = Number of lines of text (up to 10) to be printed as
!            a header describing the calculation at the beginning of
!            the cem97n.res file.
!  header()  Lines of text as described above (up to 72 characters per
!            line).
!
! ======================================================================
!
!    Edited by A. J. Sierk,  LANL  T-2  February, 1996.
!    Edited by A. J. Sierk,  LANL  T-2  November-December, 1997.
!    Modified by AJS, March, 1999.
!    Modified by SGM at 06/16/2000
!    Modified by KKG at 10/29/2001
!    Modified by A. J. Sierk, LANL T-16, October, 2003.
!    Modified by KKG , Feb, 2004
!    Comments modified by AJS, May, 2004.
!    Modified by KKG , June, October, 2004
!    Edited by AJS, January, 2005.
!    Modified by AJS, January, 2005.
!    Modified by AJS, July-September, 2005.
!    Updated copyright notice, AJS, February, 2006.
!    Edited by AJS, LANL T-2, December, 2011 - February, 2012.
!    Modified by LMK, XCP-3, June 2012-July 2013 (expand preeq).
!    Edited by LMK, XCP-3, July 2013 (included error protection).
!    Edited by CMJ, XCP-3, 2016-2017 (LAQGSM expansion)
!
! ======================================================================

    use, intrinsic:: iso_fortran_env, only: int32, int64, real32, real64
    use gsm_params, only: zro, parName, numParticleNames

    ! For tallying stuff:
    use hist_mod, only: hist_print, hist_init

! ======================================================================

    implicit none
    class(GSM),           intent(inout) :: gsmObj
    class(GSMProjectile), intent(inout), target :: projNucleus
    class(GSMTarget),     intent(inout), target :: targNucleus
    class(GSMOutput),     intent(inout), target :: output
    integer(int32),       intent(in   ) :: maxProgeny

    integer(int32) :: icase  = 0_int32
    integer(int32) :: iclose = 0_int32
    real(real32)   :: dtime, utime(2)
    real(real64)   :: thour, tsec, ttmin, utim2, utim3, xxdum

! ======================================================================

    integer(int32) :: incFlag = defaultINCFlag
    real(real64)   :: initialEnergy = defaultKinEnergy

    ! Reaction information (to be passed in by client):
    type(gsmprogeny), dimension( maxProgeny ) :: progenybnk
    type(GSMResults) :: results

    ! Create containers for reaction-specific data and physics models:
    type(GSMReaction)       :: gsmRxn

! ======================================================================

    results = newGSMResults( progenyBnk )
    results%tallySim = .TRUE.

    ! CHECK FOR ERRORS WITH OBJECT, DATA, AND RESULTS CONSTRUCTION:
    ! (Ensure object was constructed:)
    call gsmObj%validateGSMState(results)
    if ( results%simState /= successfulSingleEvent ) return

    ! Validate the output object
    call output%validateOutputObj()

    ! Stops the Calculation
    if ( trim(projNucleus%particleName) == stopName ) then
       if (iclose.ne.1) then
          write(gsmObj%io%message,1900) ! empty line
          call gsmObj%io%print(2, 5, gsmObj%io%message)
          xxdum = dtime (utime)
          thour = dble(int (utime(1)/3600.d0))
          utim2 = (dble(utime(1)) - 3600.d0*thour)
          ttmin = dble(int (utim2/60.d0))
          utim3 = (utim2 - 60.d0*ttmin)
          tsec = utim3
          if (thour > zro) then
             write (31, 1000) thour, ttmin, tsec
             write (gsmObj%io%message, 1000) thour, ttmin, tsec
          else
             write (31, 1100) ttmin, tsec
             write (gsmObj%io%message, 1100) ttmin, tsec
          endif
          call gsmObj%io%print(2, 5, gsmObj%io%message)
          write(gsmObj%io%message,1900) ! empty line
          call gsmObj%io%print(2, 5, gsmObj%io%message)
          close (31)
          iclose = 1
       endif
       close (16)
       return
    endif

    ! Fill in blanks for the projectile/target nucleus:
    call gsmObj%formNuclei(projNucleus, targNucleus)
    initialEnergy = projNucleus%kinEnergy
    incFlag = gsmObj%inquireINCModel(projNucleus, targNucleus)

    select case( incFlag )
       case( sDCMFlagged )
          write(gsmObj%io%message, 4900) "Standard DCM"
       case( mDCMFlagged )
          write(gsmObj%io%message, 4900) "Modified DCM"
       case default
          incFlag = sDCMFlagged
          write(gsmObj%io%message, 4900) "Standard DCM"
    end select
    call gsmObj%io%print(3, 4, gsmObj%io%message)
    if (projNucleus%numBaryons <= 1) then
       if ( gsmObj%options%smoothTransition ) then
          write(gsmObj%io%message, 4950) "will"
       else
          write(gsmObj%io%message, 4950) "will not"
       end if
    end if
    call gsmObj%io%print(3, 4, gsmObj%io%message)

    !>>> ASSIGN TO AN OPENMP TASK
    ! Establish all modified DCM arguments:
    call gsmObj%setMDCMReaction(projNucleus, targNucleus)

    ! Point to appropriate projectile/target pairs for the reaction:
    results%initialProj => projNucleus
    results%initialTarg => targNucleus

    ! Set the frequency the event number is printed to the client:
    if ( gsmObj%options%printIncrement <= 0 ) then
       if(output%numInelasticEvents > 20000) then
          gsmObj%options%printIncrement = 10000.d0
       elseif(output%numInelasticEvents > 2000 .and. output%numInelasticEvents <= 20000) then
          gsmObj%options%printIncrement = 1000.d0
       elseif(output%numInelasticEvents > 200  .and. output%numInelasticEvents <= 2000) then
          gsmObj%options%printIncrement = 100.d0
       elseif(output%numInelasticEvents > 20   .and. output%numInelasticEvents <= 200) then
          gsmObj%options%printIncrement = 10.d0
       elseif(output%numInelasticEvents >= 1    .and. output%numInelasticEvents <= 20) then
          gsmObj%options%printIncrement = 1.d0
       else
          ! output%numInelasticEvents <= 0
          output%numInelasticEvents = 10000.d0   ! default limit if incorrect value entered
          gsmObj%options%printIncrement = 1000.d0   ! default scale set for specified limit
       endif
       write(gsmObj%io%message, 3000) gsmObj%options%printIncrement
       call gsmObj%io%print(2, 4, gsmObj%io%message)
    end if

    write(gsmObj%io%message, 1200) "Setting up data for calculation..."
    call gsmObj%io%print(5, 4, gsmObj%io%message)


!>>> ASSIGN AN OPENMP TASK
! Initialize various constants and arrays:
    call gsmObj%initial(projNucleus, output)
    icase = icase + 1
    output%maxEventAttempts = 10 * output%numInelasticEvents
    if (projNucleus%particleFlag == photonProjFlag .or. &
         & projNucleus%particleFlag == bremsProjFlag) &
         & output%maxEventAttempts = 600 * output%numInelasticEvents

20  continue

    !>>> ASSIGN THE IMMEDIATE BELOW CALLS AN OPENMP TASK (IF TIME)
    ! Resets initial values/arrays for starting/new calculation
    if(output%printPISA) call gsmObj%pisaInit()
    if(output%printHIST) call hist_init()
    if(abs(initialEnergy - projNucleus%kinEnergy) > div0Lim) then
       call gsmObj%initial(projNucleus, output)

       ! Establish all modified DCM arguments:
       call gsmObj%setMDCMReaction(projNucleus, targNucleus)
    endif


    ! Geometrical cross section using radius at which density is
    ! 0.01 x central density:
    if (projNucleus%kinEnergy > projNucleus%kinEnergyMax) return
    iclose = 0

    ! Construct reaction specific data and models:
    ! NOTE: This must be done AFTER the incident energy is incremented for
    !       incident energy spectrum calculations
    gsmRxn    = gsmObj%setupReaction(projNucleus, targNucleus, results, output)
    gsmRxn%incFlag = incFlag


    ! Print level density multiplier information:
    write(gsmObj%io%message, 1900)
    ! Target af/an ratios
    call gsmObj%io%print(3, 5, gsmObj%io%message)
    write(gsmObj%io%message, 2600) "target"
    call gsmObj%io%print(3, 4, gsmObj%io%message)
    write(gsmObj%io%message, 2650) targNucleus%numBaryons, targNucleus%numProtons, &
         & targNucleus%afMultiplier, targNucleus%czMultiplier
    call gsmObj%io%print(3, 4, gsmObj%io%message)
    ! Projectile af/an ratios
    if ( projNucleus%numBaryons > 1 ) then
       write(gsmObj%io%message,2600) "projectile"
       call gsmObj%io%print(3, 4, gsmObj%io%message)
       write(gsmObj%io%message, 2650) dble(projNucleus%numBaryons), dble(projNucleus%numProtons), &
            & projNucleus%afMultiplier, projNucleus%czMultiplier
       call gsmObj%io%print(3, 4, gsmObj%io%message)
    end if

    ! Print setup time for calculation:
    write(gsmObj%io%message, 1900)
    call gsmObj%io%print(2, 5, gsmObj%io%message)
    xxdum = dtime (utime)
    thour = dble(int (utime(1)/3600.d0))
    utim2 = (dble(utime(1)) - 3600.d0*thour)
    ttmin = dble(int (utim2/60.d0))
    utim3 = (utim2 - 60.d0*ttmin)
    tsec = utim3
    if (thour > zro) then
       write(gsmObj%io%message, 1050) thour, ttmin, tsec
    else
       write(gsmObj%io%message, 1150) ttmin, tsec
    endif
    call gsmObj%io%print(2, 5, gsmObj%io%message)
    write(gsmObj%io%message, 1900)
    call gsmObj%io%print(2, 5, gsmObj%io%message)


    ! Write that actual calculation is starting:
    write(gsmObj%io%message, 1900)
    call gsmObj%io%print(5, 5, gsmObj%io%message)
    write(gsmObj%io%message, 1200) "Starting simulation of spallation physics..."
    call gsmObj%io%print(5, 4, gsmObj%io%message)

!  =====================================================================

    ! Perform simulation
    call gsmObj%simulateEvents( gsmRxn, output%numInelasticEvents, output%maxEventAttempts )

    ! Storing data into output file
    call fdate (output%date)
    call gsmObj%prinp (projNucleus, targNucleus, output, gsmRxn%outData, icase)   ! Print license
    call gsmObj%typeout(projNucleus, targNucleus, output, gsmRxn%outData)   ! Print results
    if (output%printHist) call hist_print(gsmRxn%outData%ncas)   ! Print histogram information

    ! Check if energy needs incremented
    if (projNucleus%dKinEnergy > 0.0_real64) then

        ! Incrementing energy, re-run calculation for that incident energy
        projNucleus%kinEnergy = projNucleus%kinEnergy + projNucleus%dKinEnergy/1000.d0

        ! Delete auxiliary file
        call flush (31)
        go to 20
    end if

    return
! ======================================================================
1000 format (' Elapsed CPU Time (computation) = ',f4.0,' hr, ',f3.0,' min, and ', &
         & f6.3,' sec.')
1050 format (' Elapsed CPU Time (setup) = ',f4.0,' hr, ',f3.0,' min, and ', &
         & f6.3,' sec.')
1100 format (' Elapsed CPU Time (computation) = ',f3.0,' min and ', &
         & f6.3,' sec.')
1150 format (' Elapsed CPU Time (setup) = ',f3.0,' min and ',f6.3,' sec.')
1200 format (A)
1900 format ("")
2600 format("The af/an fission ratios for the ", a, " particle")
2650 format("   (A=", f6.2, ", Z=", f6.2, ") are ", f5.2, " and ", f6.2, &
          & ", respectively.")
3000 format ("The event number will be printed in increments of ", i7, ".")
4900 format ("The ", A, " will be used for GSM's INC model.")
4950 format ("   A smooth INC transition energy ", A, " be used for the simulation.")
! ======================================================================
  end subroutine gsmMain
