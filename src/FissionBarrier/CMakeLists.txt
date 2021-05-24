# ============================================================================ #
#
# CMakeLists.txt for the Fission_Barrier sub-model library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "Fission_Barrier")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/fissionBarrierParams.f90
  ${CMAKE_CURRENT_LIST_DIR}/fissionBarrierClass.f90
  )

# List the sub-model's dependencies
set(dependencies
  Molnix
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #
