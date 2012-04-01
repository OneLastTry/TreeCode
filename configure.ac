#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([Treecode], [0.1], [stefans.mezulis@gmail.com])
AC_CONFIG_SRCDIR([TreeCode2/Particle.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
LT_INIT

AC_LANG([C++])


AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for header files.
AC_CHECK_HEADERS([cstdlib cmath])
BOOST_REQUIRE([1.49])
# We also want boost test, and program options
BOOST_TEST
BOOST_PROGRAM_OPTIONS

#Read eigen dir from command line if supplied
AC_ARG_WITH([eigen-dir],
		#AS_HELP_STRING is just a pretty formatter
        AS_HELP_STRING([--with-eigen-dir=EIGEN_DIR],
        [Force location of Eigen]),
        [
		 if test -d "$withval"; then
			 CPPFLAGS="$CPPFLAGS -I$withval"
         else
        	AC_MSG_ERROR(--with-eigen-dir expected directory name)
         fi
        ],
        [CPPFLAGS="$CPPFLAGS"]
)
AC_CHECK_HEADER([Eigen/Dense], [], [AC_MSG_ERROR([Eigen is required. Use --with-eigen-dir to set location])])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([floor sqrt])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
