#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pTurbFields], [0.0], [jgraha8@gmail.com])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/fields.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_FC

AC_FC_SRCEXT([f90])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])
AC_CHECK_HEADERS([autofd.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_LANG_PUSH([Fortran])
#AC_CHECK_LIB([matsolv], [matsolv_lu_decomp],[],[AC_MSG_ERROR([matsolv library not found])])
AC_LINK_IFELSE(	
		[AC_LANG_PROGRAM([],
		[use matsolv; real(matsolv_rp) LU(3,3); call matsolv_lu_decomp(LU)])],
		[],
		[AC_MSG_ERROR([matsolv library not found])])
AC_LANG_POP([Fortran])

AC_CHECK_LIB([autofd],[autofd_stencil_coefs_c],[],[AC_MSG_ERROR([autofd library not found])])

AC_OUTPUT
