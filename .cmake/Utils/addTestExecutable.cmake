# ============================================================================ #
#
# Adds an executable to the /bin/test/ directory
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
include (Prepend)
include (addProjectTest)
PreventInSourceBuild()

macro(add_test_executable target sources dependencies)

  # Add the exectable, link, and install it:
  set(execTarget "tst${target}")
  add_executable(${execTarget} EXCLUDE_FROM_ALL ${sources})
  # target_compile_options(${execTarget} BEFORE PRIVATE "${CMAKE_Fortran_FLAGS_DEBUG}")

  # If dependencies exist, add them:
  if (dependencies)
    Prepend(mod_dependencies "${PROJECT_NAME}_" ${dependencies})
    target_link_libraries(${execTarget} ${mod_dependencies})
  endif()

  if (FALSE)
  # Install the executable
  install(TARGETS ${execTarget}
    CONFIGURATIONS ${ALLOWED_CONFIGURATIONS}
    COMPONENT "${execTarget} (test executible)"
    RUNTIME DESTINATION "${CMAKE_BINARY_DIR}/test/bin"
    LIBRARY DESTINATION "${CMAKE_BINARY_DIR}/test/lib"
    ARCHIVE DESTINATION "${CMAKE_BINARY_DIR}/test/lib"
    INCLUDES DESTINATION "${CMAKE_BINARY_DIR}/test/include"
    ${ARGN}
    )
  endif ()

  # Add a test named according to the executable
  add_project_test(${execTarget})

endmacro()

# ============================================================================ #

