# ============================================================================ #
#
# Defines a macro for to add a test directory
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

set(SUPPORTED_TEST_DIR_NAMES
  test
  tests
  )

macro (add_test_directory)
  if (${PROJECT_NAME}.TEST)
    foreach(test in ${SUPPORTED_TEST_DIR_NAMES})
      if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/${test}")
        add_subdirectory(${test})
      endif()
    endforeach()
  endif()
endmacro()

# ============================================================================ #

