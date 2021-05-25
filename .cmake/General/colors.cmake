# ============================================================================ #
#
# Defines various colors available for use
#
# ============================================================================ #
cmake_minimum_required (VERSION 3.8.0)

# Specify escape character
STRING(ASCII 27 Esc)

# Specify standard colors
SET(ColorReset "${Esc}[m")
SET(White      "${Esc}[97")
SET(Black      "${Esc}[30m")
SET(Red        "${Esc}[31m")
SET(Green      "${Esc}[32m")
SET(Yellow     "${Esc}[33m")
SET(Blue       "${Esc}[34m")
SET(Magenta    "${Esc}[35m")
SET(Cyan       "${Esc}[36m")

# Specify light colors
SET(LightRed     "${Esc}[91m")
SET(LightGreen   "${Esc}[92m")
SET(LightYellow  "${Esc}[93m")
SET(LightBlue    "${Esc}[94m")
SET(LightMagenta "${Esc}[95m")
SET(LightCyan    "${Esc}[96m")
