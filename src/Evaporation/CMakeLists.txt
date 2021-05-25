# ============================================================================ #
#
# CMakeLists.txt for the Evaporation sub-model library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "Evaporation")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/Parameters/Parameters.f90
  ${CMAKE_CURRENT_LIST_DIR}/Data/PhysicsData.f90
  ${CMAKE_CURRENT_LIST_DIR}/ModelData/evaporationDataClass.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaporationClass.f90
  )

# List the sub-model's dependencies
set(dependencies
  Molnix
  Fermi_Breakup
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #
