
! ==============================================================================
!
! This file contains all of the data types used by the standard Dubna Cascade
! Model's Data Class
!
!
! Written by CMJ, XCP-3, 12/2018
!
! ==============================================================================

  ! Create an IO type
  type, private :: sDCMDataIO
     private
     procedure(IOHandler), private, nopass, pointer :: print   => printSDCMData
     character(LEN=512),   private                  :: message =  ""
  end type sDCMDataIO


  ! Options for the data class
  type, public :: sDCMDataOptions
     private
     integer(int32), public  :: numZones = numZonesDefault
     
     ! Note: Target data is based on the nucleon density formula: rho = rho0{ 1 - exp[(r-r0 A**(1/3) ) / expDenom] }
     real(real64),   public  :: expDenom = expDenomDefault
     real(real64),   public  :: r0       = r0Default
     real(real64),   public  :: maxRad   = maxRadDefault
     real(real64),   public, dimension(maxZonesAllowed) :: aveDen = aveDenDefault

  end type sDCMDataOptions



  ! Create a "target" object
  type, private :: sDCMDataTarget
     private
     ! Target nucleus composition
     real(real64),   private :: numBaryons = 0.0_real64     ! Number of baryons (i.nucleons) in the nucleus
     real(real64),   private :: numProtons = 0.0_real64     ! Number of protons in the nucleus
     real(real64),   private :: aTargThrd  = 0.0_real64     ! A**(1/3) for target

     ! Regarding pion potentials and separation/binding energies
     real(real64),   private :: pionPote  = pionPoteDefault
     real(real64),   private :: sepEnergy = sepEnergyDefault

     ! Geometric Cross Section
     real(real64),   private :: geomCrossSection = 0.0_real64   ! Geometric cross section [fm^2]


     ! Regarding different nuclear "zones"
     ! NOTE: The trawling effect (depletion) is NEGLECTED within the Standard DCM. The below values will NOT change!
     ! rsm / rbig
     real(real64), private, dimension(maxZonesAllowed) :: zoneBoundR      = 0.0_real64   ! Bound of a nuclear zone
     real(real64), private, dimension(maxZonesAllowed) :: zoneBRFrac      = 0.0_real64   ! Normalized (to zoneBoundR's outer used bound) radius
     ! rhop / rhon
     real(real64), private, dimension(maxZonesAllowed) :: protonDensity   = 0.0_real64   ! Average nuclear density of protons  in each zone
     real(real64), private, dimension(maxZonesAllowed) :: neutronDensity  = 0.0_real64   ! Average nuclear density of neutrons in each zone
     ! af
     real(real64), private, dimension(maxZonesAllowed) :: coulombPote     = 0.0_real64   ! Average Coulomb potential           in each zone
     ! tfp / tfn
     real(real64), private, dimension(maxZonesAllowed) :: protFermiMom    = 0.0_real64   ! Average Fermi momenta of protons    in each zone
     real(real64), private, dimension(maxZonesAllowed) :: neutFermiMom    = 0.0_real64   ! Average Fermi momenta of neutrons   in each zone
   contains
  end type sDCMDataTarget
