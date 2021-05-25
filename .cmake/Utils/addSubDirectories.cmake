# ============================================================================ #
#
# Defines a macro for add several subdirectories at once
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

macro (add_subdirectories)
  foreach(dir ${ARGV})
    add_subdirectory(${dir})
  endforeach()
endmacro()

# ============================================================================ #

