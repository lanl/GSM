# ============================================================================ #
#
# CMakeLists.txt for the gsmOther sub-model library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

# Set sub-model's name
set(model_name "GSM_Deprecated")

# ============================================================================ #

# Load the sources used for the sub-model
set(sources
  ${CMAKE_CURRENT_LIST_DIR}/../GeneralizedSpallation/gsm_params.f90
  ${CMAKE_CURRENT_LIST_DIR}/gsm_derived_types.f90
  ${CMAKE_CURRENT_LIST_DIR}/gsm_options.f90
  ${CMAKE_CURRENT_LIST_DIR}/defaults_mod.f90
  )

# List the sub-model's dependencies
set(dependencies
  )

# ============================================================================ #

# Create the sub-model's library
add_sublibrary(${model_name} "${sources}" "${dependencies}")

# ============================================================================ #
