#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Initialization stuff I don't understand.
AC_PREREQ([2.68])
AC_INIT(OpenMOC, 0.0, lululi@mit.edu)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_MAKE_SET
PKG_PROG_PKG_CONFIG

# Checks for libraries.
AC_CHECK_LIB([expat], [XML_ParserCreate])
AC_CHECK_LIB([silo], [silo_Init])
AC_CHECK_LIB([Magick++], [MagickPlusPlusDestroyMagick])

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h string.h sys/time.h])

# Uses pkg-config to find alternate library install paths
PKG_CHECK_MODULES(MAGICKXX, [Magick++])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([clock_gettime floor gettimeofday mkdir sqrt strdup strndup])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
