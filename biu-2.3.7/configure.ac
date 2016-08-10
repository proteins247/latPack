#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
AC_PREREQ([2.59])
AC_INIT([BIU], [2.3.7], [http://www.bioinf.uni-freiburg.de/SW/BIU/])
AC_CONFIG_SRCDIR([src/biu/OptionParser.hh])
AC_CONFIG_AUX_DIR([.])
AC_CANONICAL_HOST

# handle windows specific gnu linker flags
lt_enable_auto_import=""
case "$host_os" in 
	cygwin* | mingw* | cegcc*)
		AM_LDFLAGS="-Wl,--enable-auto-import $AM_LDFLAGS"
esac


# Checks for programs.
#AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_RANLIB
AC_CHECK_PROG([HAVE_PKG_CONFIG], [pkg-config], [yes],[no])

#AM_INIT_AUTOMAKE([subdir-objects])

##########################################################################
# Checks for programs needed for tests
##########################################################################

# Check for 'diff' and get full path.
AC_ARG_VAR([DIFF],[the 'diff' program to use for test output comparison])
AC_PATH_PROG([DIFF],[diff],[])
if test "x$DIFF" = "x"; then
	AC_MSG_NOTICE([==> diff command not found!])
	AC_MSG_NOTICE([==> Set DIFF variable if present in non-standard path!])
	AC_MSG_ERROR([diff is mandatory to run the tests : will stop here!])
fi

# Check for 'cat' and get full path.
AC_ARG_VAR([CAT],[the 'cat' program used for printing test output files])
AC_PATH_PROG([CAT],[cat],[])
if test "x$CAT" = "x"; then
	AC_MSG_NOTICE([==> cat command not found!])
	AC_MSG_NOTICE([==> Set CAT variable if present in non-standard path!])
	AC_MSG_ERROR([cat is mandatory to run the tests : will stop here!])
fi

# Check for 'printf' and get full path.
AC_ARG_VAR([PRINTF],[the 'printf' program used to print test information])
AC_PATH_PROG([PRINTF],[printf],[])
if test "x$PRINTF" = "x"; then
	AC_MSG_NOTICE([==> printf command not found!])
	AC_MSG_NOTICE([==> Set PRINTF variable if present in non-standard path!])
	AC_MSG_ERROR([printf is mandatory to run the tests : will stop here!])
fi

##########################################################################

AC_MSG_CHECKING(whether to build with debug information)
debuger=no
AC_ARG_ENABLE(debug,
	AC_HELP_STRING(
	    [--enable-debug],
	    [enable debug data generation (def=no)]
	),
  	debuger="$enableval"
)
AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
	AC_DEFINE(DEBUG)
	CXXFLAGS="$CXXFLAGS -g -Wall -Werror -Wno-uninitialized -O0 -ansi -pedantic -march=native "
else
	AC_DEFINE(NDEBUG)
	CXXFLAGS="$CXXFLAGS -O2 -march=native "
fi

AC_MSG_CHECKING([wether to build with gcov code coverage information (GCC only)])
gcov=no
AC_ARG_ENABLE(gcov,
	AC_HELP_STRING(
	    [--enable-gcov],
	    [enable gcov code coverage data generation (GCC only) (def=no)]
	),
  	gcov="$enableval"
)
AC_MSG_RESULT($gcov)
if test x"$gcov" = x"yes"; then
	CXXFLAGS="$CXXFLAGS -fprofile-arcs -ftest-coverage"
	LDFLAGS="$LDFLAGS -fprofile-arcs -ftest-coverage"
	CCACHE_DISABLE=1
fi

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(BIU, misc/doxygen/doxygen.config, doc)
#
AM_INIT_AUTOMAKE([subdir-objects])

AC_LANG([C++])


#
# Checks for libraries.
#
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])
#
# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

################################################################################
# CHECK FOR AVAILABLE HASH IMPLEMENTATIONS
################################################################################

# check for unordered_map
AC_MSG_CHECKING([for std::unordered_map])
AC_COMPILE_IFELSE(
	[
		AC_LANG_PROGRAM([[#include <unordered_map>]],
			[[std::unordered_map< int, int > myMap; return 0;]])
	], 
	[	
		AC_MSG_RESULT(yes) 
		# define if the unorderd_map is available
		STD_UNORDERED_MAP_DEFINE="#define HAVE_UNORDERED_MAP 1"
	],
	[ 	
		AC_MSG_RESULT(no)
		STD_UNORDERED_MAP_DEFINE="#define HAVE_UNORDERED_MAP 0"
	]) # unordered_map
# distribute result
AC_SUBST(STD_UNORDERED_MAP_DEFINE, ${STD_UNORDERED_MAP_DEFINE})
	
# check for tr1::unordered_map
AC_MSG_CHECKING([for std::tr1::unordered_map])
AC_COMPILE_IFELSE(
	[
		AC_LANG_PROGRAM([[#include <tr1/unordered_map>]],
			[[std::tr1::unordered_map< int, int > myMap; return 0;]])
	], 
	[	
		AC_MSG_RESULT(yes) 
		# define if the unorderd_map is available
		STD_TR1_UNORDERED_MAP_DEFINE="#define HAVE_TR1_UNORDERED_MAP 1"
	],
	[ 	
		AC_MSG_RESULT(no)
		# define if the unorderd_map is available
		STD_TR1_UNORDERED_MAP_DEFINE="#define HAVE_TR1_UNORDERED_MAP 0"
	]) # tr1::unordered_map
# distribute result
AC_SUBST(STD_TR1_UNORDERED_MAP_DEFINE, ${STD_TR1_UNORDERED_MAP_DEFINE})
	
# check for gnu hash_map
AC_MSG_CHECKING([for GNU hash_map])
AC_COMPILE_IFELSE(
	[
		AC_LANG_PROGRAM([[#include <ext/hash_map>]],
			[[__gnu_cxx::hash_map< int, int > myMap; return 0;]])
	], 
	[	
		AC_MSG_RESULT(yes) 
		# define if the GNU hash_map is available
		GNU_HASH_MAP_DEFINE="#define HAVE_GNU_HASH_MAP 1"
	],
	[ 	
		AC_MSG_RESULT(no)
		# define if the GNU hash_map is available
		GNU_HASH_MAP_DEFINE="#define HAVE_GNU_HASH_MAP 0"
	]) # gnu::hash_map
	
# distribute prefix
AC_SUBST(GNU_HASH_MAP_DEFINE, ${GNU_HASH_MAP_DEFINE})

################################################################################


# distribute prefix
AC_SUBST(PREFIX_PATH, ${prefix})
AC_SUBST(HASH_MAP_DEFINE, ${HASH_MAP_DEFINE})

# distribute the configure based compiler and linker setups
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([ Makefile ])
AC_CONFIG_FILES([ src/doxygen.hh])
AC_CONFIG_FILES([ src/biu/HashMap.hh ])
AC_CONFIG_FILES([ src/Makefile ])
AC_CONFIG_FILES([ misc/Makefile misc/doxygen/Makefile ])
AC_CONFIG_FILES([ tests/Makefile ])

# set up package information if pkg-config is present
if test "x$HAVE_PKG_CONFIG" = "xyes"; then
AC_CONFIG_FILES([ biu.pc ])
fi

AC_OUTPUT
