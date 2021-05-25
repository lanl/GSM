# ============================================================================ #
#
# CMakeLists.txt for the Fermi_Breakup sub-model library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "Fermi_Breakup")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/Parameters/Parameters.f90
  ${CMAKE_CURRENT_LIST_DIR}/Data/PhysicsData.f90
  ${CMAKE_CURRENT_LIST_DIR}/fermiBreakupClass.f90
  )

# List the sub-model's dependencies
set(dependencies
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #
