# ============================================================================ #
#
# CMakeLists.txt for the abstracting parallelization
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "Parallel")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/OpenMPAbstraction.f90
  # ${CMAKE_CURRENT_LIST_DIR}/MPI.f90
  )

# List the sub-model's dependencies
set(dependencies
  Contracts
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #
