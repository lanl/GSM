# ============================================================================ #
#
# CMakeLists.txt for the General sub-model library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "General")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/numbers.f90
  ${CMAKE_CURRENT_LIST_DIR}/fund_data.f90
  ${CMAKE_CURRENT_LIST_DIR}/coulomb_barrier.f90
  ${CMAKE_CURRENT_LIST_DIR}/inverse_x_section.f90
  ${CMAKE_CURRENT_LIST_DIR}/hist_mod.f90
  )

# List the sub-model's dependencies
set(dependencies
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #
