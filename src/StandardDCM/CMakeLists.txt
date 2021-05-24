# ============================================================================ #
#
# CMakeLists.txt for the Standard_DCM sub-model library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "Standard_DCM")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/Parameters/Parameters.f90
  ${CMAKE_CURRENT_LIST_DIR}/Parameters/DataClassParameters.f90
  ${CMAKE_CURRENT_LIST_DIR}/Data/PhysicsData.f90
  ${CMAKE_CURRENT_LIST_DIR}/ModelData/standardDCMDataClass.f90
  ${CMAKE_CURRENT_LIST_DIR}/standardDCMClass.f90
  )

# List the sub-model's dependencies
set(dependencies
  Molnix
  Photon_Event_Generator
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #
