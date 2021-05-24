# ============================================================================ #
#
# CMakeLists.txt for the gsmMain sub-model library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "GSM")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  # Modules for various objects
  ${CMAKE_CURRENT_LIST_DIR}/objects/BremsPhoton/BremsPhoton.f90
  ${CMAKE_CURRENT_LIST_DIR}/objects/Particles/Particles.f90
  ${CMAKE_CURRENT_LIST_DIR}/objects/EventData/EventData.f90
  ${CMAKE_CURRENT_LIST_DIR}/objects/OutputData/OutputData.f90
  ${CMAKE_CURRENT_LIST_DIR}/gsm_params.f90
  ${CMAKE_CURRENT_LIST_DIR}/generalizedSpallationData.f90
  ${CMAKE_CURRENT_LIST_DIR}/generalizedSpallationClass.f90
  )

# List the sub-model's dependencies
set(dependencies
   # General physics:
   General
   Molnix
   Fission_Barrier
   Photon_Event_Generator
   Random_Number_Generator
   # Physics models:
   Standard_DCM
   Modified_DCM
   Coalescence
   Fermi_Breakup
   Preequilibrium
   Evaporation
   # GSM Specific:
   GSM_Deprecated
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #

