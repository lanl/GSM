# ============================================================================ #
#
# CMakeLists.txt file for the main source directory
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# General Fortran objects used:
add_subdirectories(
  Harness
  Utilities
  RandomMC
  )


# Primary physics:
add_subdirectories(
  General
  PhotonEventGenerator
  Molnix
  FissionBarrier
  )


# GSM sub-models:
add_subdirectories(
  StandardDCM
  ModifiedDCM
  Coalescence
  FermiBreakup
  Preequilibrium
  Evaporation
  )


# Create the main GSM object:
add_subdirectories(
  GSM_Setup
  GeneralizedSpallation
  )


# Clients of GSM:
add_subdirectories(
  DriverGSM
  )

# ============================================================================ #

