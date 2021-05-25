# ============================================================================ #
#
# Defines a function to prepend text to a variable or variable list
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

function(Prepend var prefix)

  # Establish and create a new list - "list_var" - from any additional arguments
  set(list_var "")
  foreach(f ${ARGN})
    list(APPEND list_var "${prefix}${f}")
  endforeach()

  # Update the provided variable - var - in the calling file
  set(${var} "${list_var}" PARENT_SCOPE)

endfunction()

# ============================================================================ #

