# ============================================================================ #
#
# CMakeLists.txt for the GSM Driver sub-model library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "GSM_Driver")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/standaloneGSM.f90
  )

# List the sub-model's dependencies
set(dependencies
  Contracts
  Parallel
  GSM
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #
