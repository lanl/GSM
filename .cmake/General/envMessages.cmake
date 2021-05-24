# ============================================================================ #
#
# Defines some global variables for project messaging
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

if (${CMAKE_MAJOR_VERSION} VERSION_LESS "3.15.3")
  set(envNOTICE "" CACHE STRING "Notice message type" FORCE)
  set(envVERBOSE "STATUS" CACHE STRING "Verbose message type" FORCE)
  set(envDEBUG "STATUS" CACHE STRING "Debug message type" FORCE)
else ()
  set(envNOTICE "NOTICE" CACHE STRING "Notice message type" FORCE)
  set(envVERBOSE "VERBOSE" CACHE STRING "Verbose message type" FORCE)
  set(envDEBUG "DEBUG" CACHE STRING "Debug message type" FORCE)
endif ()

# Prevent variables from appearing in the list:
mark_as_advanced(envNOTICE)
mark_as_advanced(envVERBOSE)
mark_as_advanced(envDEBUG)

# ============================================================================ #
#
