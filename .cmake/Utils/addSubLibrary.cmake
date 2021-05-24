# ============================================================================ #
#
# Defines a macro that creates a sub-library
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
include (Prepend)

macro (add_sublibrary target sources dependencies)

  message(STATUS
    "-- Configuring the ${Yellow}${target}${ColorReset} sub-library."
    )

  # Perform general setup of a library
  set(subtarget "${target}")

  # Obtain library properties
  set(lib_type "STATIC")
  if ( ${PROJECT_NAME}.SHARED )
    set(lib_type "SHARED")
  endif()

  # Create the general structure for the library
  add_library(${subtarget} "${libtype}" "")
  add_library(${PROJECT_NAME}::${target} ALIAS ${subtarget})

  # Set library sources
  target_sources(${subtarget} PRIVATE
    ${sources}
    )

  # Link to all given dependencies
  if (dependencies)
    target_link_libraries(${subtarget}
      PRIVATE
      ${dependencies}
      )
  endif()

  # Add a test directory for the sub-model
  add_test_directory()


endmacro()

# ============================================================================ #

