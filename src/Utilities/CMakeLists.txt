# ============================================================================ #
#
# CMakeLists.txt for the Utilities used by GSM
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)
PreventInSourceBuild()

add_subdirectories(
  Contract
  Parallel
  Logger
  LogManager
  Attributes
  )


