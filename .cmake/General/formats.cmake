# ============================================================================ #
#
# Defines various text formats available for use
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Specify escape character
STRING(ASCII 27 Esc)

# Specify standard colors
SET(FmtReset   "${Esc}[0m")
SET(Bold       "${Esc}[1m")
SET(Dim        "${Esc}[2m")
SET(UnderLine  "${Esc}[4m")
SET(Blink      "${Esc}[5m")
SET(Reverse    "${Esc}[7m")
SET(Hidden     "${Esc}[8m")

