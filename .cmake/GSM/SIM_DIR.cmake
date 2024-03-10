# ============================================================================ #
#
# Copies the executable and all required data files/tables to a single directory
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

if (${PROJECT_NAME}.SIM_DIR)

  # Create a label for the directory:
  set(simDir "my_${PROJECT_NAME}")

  # Add a Makefile command
  add_custom_target(${simDir} ALL

    # Create the directory
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/${simDir}"

    # Copy the data files/tables
    COMMAND ${CMAKE_COMMAND} -E copy
       "${CMAKE_SOURCE_DIR}/src/dataTbl/*"
       "${CMAKE_SOURCE_DIR}/${simDir}/"

    # Copy the executable
    COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_BINARY_DIR}/bin/x${PROJECT_NAME}*"
       "${CMAKE_SOURCE_DIR}/${simDir}/"
    )

  # The simulation directory will depend on any files mentioned above.
  # Ensure we add dependencies to those.
    add_dependencies(${simDir}
      "x${PROJECT_NAME}${PROJECT_VERSION_MAJOR}"
      )

endif()


