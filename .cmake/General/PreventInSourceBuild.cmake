# ============================================================================ #
#
# Prints an error message on an attempt to build inside the source directory
# tree
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
include (colors)

macro (PreventInSourceBuild)

  # Determine if in the parent directory or a sub-directory of the parent:
  string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insrc)
  get_filename_component(parent_dir ${CMAKE_SOURCE_DIR} DIRECTORY)
  string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${parent_dir}" insrc_subdir)

  # Check for CMakeLists.txt file as well
  file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)

  if (${insrc} OR ${insrc_subdir} OR EXISTS "${LOC_PATH}")
    message(FATAL_ERROR
      "\n${Red}ERROR!${ColorReset}\n"
      "     CMAKE_CURRENT_SOURCE_DIR="
      "${LightBlue}${CMAKE_CURRENT_SOURCE_DIR}${ColorReset}\n"
      "     CMAKE_CURRENT_BINARY_DIR="
      "${LightMagenta}${CMAKE_CURRENT_BINARY_DIR}${ColorReset}"

      "\n${PROJECT_NAME} does not support in-source builds!\n"
      "     You must now delete the CMakeCache.txt file and the CMakeFiles/\n"
      "     directory under the top ${PROJECT_NAME} directory or you will\n"
      "     not be able to configure correctly!"

      "\nYou must now run something like...${Cyan}\n"
      "   $ rm -r CMakeCache.txt CMakeFiles/${ColorReset}"

      "\nPlease create a different directory and configure ${PROJECT_NAME}"
      "     under that different directory, such as:${Cyan}\n"
      "   $ mkdir MY_BUILD\n"
      "   $ cd MY_BUILD\n"
      "   $ cmake [OPTIONS] ../\n${ColorReset}"
      ""
      )
  endif()
endmacro()

macro(macro_ensure_out_of_source_build project)
endmacro()
