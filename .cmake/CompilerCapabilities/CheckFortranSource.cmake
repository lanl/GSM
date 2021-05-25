# ============================================================================ #
#
# Tests if a compiler feature is supported by the compiler
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# ============================================================================ #
# Checks if the source file can be compiled
# ============================================================================ #
function(Check_Fortran_Source_Compile test_file compiler_feature)

  # Print that test will be performed
  if (NOT CMAKE_REQUIRED_QUIET)
    message (${envVERBOSE} "Testing ${compiler_feature}")
  endif ()

  # Compile the source code for the feature being tested
  try_compile (
    ${compiler_feature}
    ${CMAKE_BINARY_DIR}
    ${test_file}
    COMPILE_DEFINITIONS ${CMAKE_Fortran_FLAGS_RELEASE}
    OUTPUT_VARIABLE build_output
    )

  # Print if test was successful or not
  if (${compiler_feature})
    if (NOT CMAKE_REQUIRED_QUIET)
      message(${envVERBOSE} "Testing ${compiler_feature}: SUCCESSS")
    endif ()
    add_definitions(-D${compiler_feature})
  else ()
    if (NOT CMAKE_REQUIRED_QUIET)
      message(${envVERBOSE} "Testing ${compiler_feature}: FAILURE")
      message(${envVERBOSE} "The build messages were:\n${build_output}")
    endif ()
  endif ()
endfunction()


# ============================================================================ #
# Checks if the source code can be compiled and ran
# ============================================================================ #
function(CHECK_FORTRAN_SOURCE_RUN test_file compiler_feature)

  if (NOT CMAKE_REQUIRED_QUIET)
    message (${envVERBOSE} "Testing ${compiler_feature}")
  endif ()

  try_run (
    code_runs
    code_compiles
    ${CMAKE_BINARY_DIR}
    ${test_file}
    )

  if (${code_compiles})
    if (${code_runs} EQUAL 0)
      if (NOT CMAKE_REQUIRED_QUIET)
        message (${envVERBOSE} "Testing ${compiler_feature}: SUCCESS")
      endif ()

      add_definitions(-D${compiler_feature})
      set (${compiler_feature} 1)

    else ()

      if (NOT CMAKE_REQUIRED_QUIET)
        message (${envVERBOSE} "Performing Test ${compiler_feature}: RUN FAILURE")
      endif ()

    endif ()
  else ()
      if (NOT CMAKE_REQUIRED_QUIET)
        message (${envVERBOSE} "Performing Test ${compiler_feature}: BUILD FAILURE")
      endif ()
  endif()
endfunction()

# ============================================================================ #

