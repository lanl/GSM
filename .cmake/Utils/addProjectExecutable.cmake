# ============================================================================ #
#
# Adds an executable to the project in the /bin/ directory
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
include (Prepend)
PreventInSourceBuild()

macro(add_project_executable target sources dependencies)

  # Inform user of configuration:
  message(${envDEBUG}
    "${Green}Configuring the ${target} executable...${ColorReset}")

  # Add the exectable, link, and install it:
  add_executable(${target})
  add_executable(${PROJECT_NAME}::${target} ALIAS ${target})

  target_sources(${target} PRIVATE
    ${sources}
    )

  target_compile_definitions(
    ${target}
    PUBLIC
    ${${PROJECT_NAME}_compile_definitions}
    )

  # If dependencies exist, add them:
  if (dependencies)
    target_link_libraries(${target}
      PRIVATE
      ${dependencies}
      )
  endif()

  # Install the executable
  install(TARGETS ${target}
    CONFIGURATIONS ${SUPPORTED_CONFIGURATIONS}
    COMPONENT "${target} (executable)"
    RUNTIME DESTINATION "${CMAKE_BINARY_DIR}/bin"
    LIBRARY DESTINATION "${CMAKE_BINARY_DIR}/lib"
    ARCHIVE DESTINATION "${CMAKE_BINARY_DIR}/lib"
    INCLUDES DESTINATION "${CMAKE_BINARY_DIR}/include"
    ${ARGN}
    )

endmacro()

# ============================================================================ #

