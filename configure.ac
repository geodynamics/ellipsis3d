dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.
dnl
dnl Copyright (C) 1995 The GeoFramework Consortium
dnl
dnl This file is part of Ellipsis3D.
dnl
dnl Ellipsis3D is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License, version 2,
dnl as published by the Free Software Foundation.
dnl
dnl Ellipsis3D is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl Author:
dnl   Louis Moresi <louis.moresi@sci.monash.edu>
dnl

# $Id$

AC_PREREQ(2.59)
AC_INIT([Ellipsis3D], [1.0.0], [cig-mc@geodynamics.org], [Ellipsis3D])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_SRCDIR([Ellipsis3D.c])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([foreign])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h rpc/rpc.h stdlib.h string.h strings.h sys/resource.h sys/time.h sys/times.h time.h unistd.h])
AC_CHECK_HEADERS([rpc/xdr.h], [], [],
[#if HAVE_RPC_RPC_H
#include <rpc/rpc.h>
#endif
])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_TYPE_SIGNAL

# Checks for library functions.
AC_SEARCH_LIBS([sqrt], [m])
AC_CHECK_FUNCS([floor pow sqrt], [], [AC_MSG_FAILURE([missing math function '$ac_func'])])
AC_CHECK_FUNCS([getrusage times])
AC_SEARCH_LIBS([xdrstdio_create], [], [
    AC_DEFINE([HAVE_LIBXDR], [1], [Define to 1 if you have library routines for external data representation.])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

dnl end of configure.ac
