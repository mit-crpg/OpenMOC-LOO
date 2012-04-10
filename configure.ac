#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Initialization (note the explicit autoconf dependency)
AC_PREREQ([2.68])
AC_INIT(OpenMOC, 0.0, lululi@mit.edu)
AM_INIT_AUTOMAKE

# Be sure to add any additional subdirs here
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 src/Makefile])

# Checks for programs.
AC_PROG_CXX
AC_PROG_MAKE_SET
PKG_PROG_PKG_CONFIG

# All of the OpenMOC code is C++, and some libraries won't link with C
AC_LANG(C++)

# Searches the system for the libraries needed to build OpenMOC
AC_SEARCH_LIBS([XML_ParserCreate], [expat], [], [
  AC_MSG_ERROR([unable to find expat])
], [])

AC_SEARCH_LIBS([silo_version], [silo siloh5], [], [
  AC_MSG_ERROR([unable to find silo])
], [-lz])

# Requires some non-standard header files
AC_CHECK_HEADERS([limits.h stdlib.h string.h sys/time.h])

# Uses pkg-config to find alternate library install paths
PKG_CHECK_MODULES(MAGICKXX, [Magick++], [], [
  AC_MSG_ERROR([unable to find magick++])
])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Ensures that we have all the library functions we use in OpenMOC
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([clock_gettime floor gettimeofday mkdir sqrt strdup strndup])

AC_OUTPUT
