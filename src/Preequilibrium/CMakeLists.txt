# ============================================================================ #
#
# CMakeLists.txt for the Preequilibrium sub-model library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "Preequilibrium")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/preequilibriumParams.f90
  ${CMAKE_CURRENT_LIST_DIR}/preequilibriumData.f90
  ${CMAKE_CURRENT_LIST_DIR}/lambda_j.f90
  ${CMAKE_CURRENT_LIST_DIR}/ModelData/gammaJClass.f90
  ${CMAKE_CURRENT_LIST_DIR}/ModelData/preequilibriumDataClass.f90
  ${CMAKE_CURRENT_LIST_DIR}/preequilibriumClass.f90
  )

# List the sub-model's dependencies
set(dependencies
  General
  Molnix
  Fission_Barrier
  Fermi_Breakup
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #
