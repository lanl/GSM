# ============================================================================ #
#
# Adds a test to the project given an executable name
# (test will have the same name)
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

macro (add_project_test target)

  # message(STATUS "      Adding the ${Dim}${target}${FmtReset} test...")

  add_test(
    NAME ${target}
    COMMAND ${target}
    # --exe $<TARGET_FILE:${target}>
    # --exe ${CMAKE_BINARY_DIR}/test/bin
    ${ARGN}
    # CONFIGURATIONS ${SUPPORTED_CONFIGURATIONS}
    )

  # Now add the test executable to the build_and_test target
  add_dependencies(build_and_test ${target})

endmacro()

# ============================================================================ #

