#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
AC_PREREQ(2.59)
AC_INIT([LatPack],[1.9.1], [http://www.bioinf.uni-freiburg.de/Software/LatPack/])
AC_CONFIG_SRCDIR([src/latFit.cc])
AC_CONFIG_AUX_DIR([.])

# minimal reqired version of the BIU library
BIU_REQUIRED_VERSION=2.3.7
ELL_REQUIRED_VERSION=3.2.4

# biu version distributed with this package
LOCAL_BIU_FOLDER="biu-2.3.7"
AC_CONFIG_SUBDIRS([biu-2.3.7])

# Checks for programs.
AC_PROG_CXX
AC_CHECK_PROG([HAVE_PKG_CONFIG], [pkg-config], [yes],[no])

AC_MSG_CHECKING(whether to build in debug mode with assertions)
debuger=no
AC_ARG_ENABLE(debug,
	AC_HELP_STRING(
	    [--enable-debug],
	    [enable debug mode with assertions (def=no)]
	),
  	debuger="$enableval"
)
AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
	AC_DEFINE(DEBUG)
	CXXFLAGS="$CXXFLAGS -g -Wall -Werror -Wno-uninitialized -O0  -ansi -pedantic -march=native "
else
	AC_DEFINE(NDEBUG)
	CXXFLAGS="$CXXFLAGS -O2 -march=native "
fi


AM_INIT_AUTOMAKE([subdir-objects])

AC_LANG_CPLUSPLUS
AC_DEFINE(__cplusplus)


###############################################################################
# BEGIN DOCUMENTATION BUILD CHECK
###############################################################################
dnl Find latex and friends

AC_MSG_CHECKING([whether to build docs])
AC_ARG_WITH(docs,
        AS_HELP_STRING([--without-docs],
                           [do not try to build docs]),
           ,with_docs=yes)
AC_MSG_RESULT($with_docs)

if test "$with_docs" = "yes"; then
        AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
        if test -z "$PDFLATEX"; then
        	AC_MSG_WARN([Cannot find pdflatex in your path!])
        fi
        if test ! '(' -z "$PDFLATEX" ')'; then
        	BUILDDOC="yes"
        fi
else
        BUILDDOC="no"
fi

AM_CONDITIONAL(build_docs, [test "$BUILDDOC" = "yes"])


###############################################################################
# END DOCUMENTATION BUILD CHECK
###############################################################################


###############################################################################
# BEGIN STATIC LINKING CHECK
###############################################################################

AC_ARG_ENABLE(static-linking,
    AC_HELP_STRING(
    	[--enable-static-linking],
	    [enables static linking of executables (default=disabled)]
	),
    [enable_static_linking=$enableval],
    [enable_static_linking=no]
)

if test "$enable_static_linking" = "yes"; then
	LDFLAGS="-static $LDFLAGS"
fi
AC_MSG_CHECKING([whether to link only with static libraries])
AC_MSG_RESULT($enable_static_linking)


###############################################################################
# END STATIC LINKING CHECK
###############################################################################



###############################################################################
# BEGIN BIU LIBRARY PATH SUPPORT
###############################################################################

LIBS="-lbiu $LIBS"

# BIU library path support, if not installed in usual directories
AC_ARG_WITH(BIU,
	AC_HELP_STRING(
	    [--with-BIU=prefix],
	    [try this for the BIU-library prefix install directory]
	  ),
	  BIUPATHSET=1,
	  BIUPATHSET=0
)
if test  $BIUPATHSET = 1 ; then
	CXXFLAGS="$CXXFLAGS -I$with_BIU/include"
	LDFLAGS="-L$with_BIU/lib $LDFLAGS"
	export PKG_CONFIG_PATH="$with_BIU/lib/pkgconfig/:$PKG_CONFIG_PATH"
	LOCAL_BIU_FOLDER=""
else
	# biu version distributed with this package
	AC_SUBST(LOCAL_BIU_FOLDER)
	LOCAL_BIU_ROOT=`pwd`/$LOCAL_BIU_FOLDER;
	CXXFLAGS="-I$LOCAL_BIU_ROOT/src/ $CXXFLAGS"
	LDFLAGS="-L$LOCAL_BIU_ROOT/src/ $LDFLAGS"
fi

AM_CONDITIONAL(build_local_biu, [test "$BIUPATHSET" = 0])

###############################################################################
# END BIU LIBRARY PATH SUPPORT
###############################################################################


###############################################################################
# BEGIN ELL LIBRARY PATH SUPPORT
###############################################################################

# ELL library path support, if not installed in usual directories
AC_ARG_WITH(ELL,
	AC_HELP_STRING(
	    [--with-ELL=prefix],
	    [try this for the ELL-library prefix install directory. Mandatory for latFold program.]
	  ),
	  ELLPATHSET=1,
	  ELLPATHSET=0
)
if test  $ELLPATHSET = 1 ; then
	CXXFLAGS="$CXXFLAGS -I$with_ELL/include"
	LDFLAGS="-L$with_ELL/lib $LDFLAGS"
	export PKG_CONFIG_PATH="$with_ELL/lib/pkgconfig/:$PKG_CONFIG_PATH"
fi

###############################################################################
# END ELL LIBRARY PATH SUPPORT
###############################################################################


###############################################################################
# BEGIN GSL LIBRARY PATH SUPPORT
###############################################################################

# GSL library path support, if not installed in usual directories
AC_ARG_WITH(GSL,
	AC_HELP_STRING(
	    [--with-GSL=prefix],
	    [try this for a non-standard GNU scientific library (GSL) install prefix. Mandatory for latFit program.]
	  ),
	  GSLPATHSET=1,
	  GSLPATHSET=0
)
if test  $GSLPATHSET = 1 ; then
	CXXFLAGS="-I$with_GSL/include $CXXFLAGS"
	LDFLAGS="-L$with_GSL/lib $LDFLAGS"
fi

###############################################################################
# END GSL LIBRARY PATH SUPPORT
###############################################################################



###############################################################################
# BEGIN FOLDING MODULE CHECK
###############################################################################

AC_ARG_ENABLE(latFold,
    AC_HELP_STRING(
    	[--enable-latFold],
	    [enable compilation of latFold program (default=disabled)]
	),
[enable_latfold=$enableval],
[enable_latfold=no]
)

AM_CONDITIONAL(with_latfold, 
	[test "$enable_latfold" = "yes"])
	
if test "$enable_latfold" = "yes"; then
	LIBS="-lell $LIBS"
fi

###############################################################################
# END FOLDING MODULE CHECK
###############################################################################

# Checks for libraries and header files.


# check for stl headers
AC_CHECK_HEADERS([vector])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

###############################################################################
# BEGIN DEPENDENCY CHECK
###############################################################################

###############################################################################
# BEGIN BIU CHECK
###############################################################################
# check for biu headers
AC_MSG_CHECKING([for the BIU C++ library])
	AC_COMPILE_IFELSE([AC_LANG_SOURCE([[#include <biu/VirtualList.hh>]])], 
     [ 
       AC_MSG_RESULT(yes)
       BIUNOTFOUND=0;
     ],
	 [
       AC_MSG_RESULT(no)
       BIUNOTFOUND=1;
     ]
)


if test "$enable_latfold" = "yes"; then
	# check for ELL headers
	AC_MSG_CHECKING([for the ELL C++ library headers])
	AC_COMPILE_IFELSE(
	     [AC_LANG_PROGRAM([[ #include <ell/State.hh> ]])],
	     [ 
	       AC_MSG_RESULT(yes)
	       ELLNOTFOUND=0;
	     ],
	     [
	       AC_MSG_RESULT(no)
	       ELLNOTFOUND=1;
	     ]
	)
else
   ELLNOTFOUND=0;
fi

# store current *FLAGS and merge with AM_*FLAGS for compilation and linker check   
OLD_CXXFLAGS=$CXXFLAGS;
OLD_CPPFLAGS=$CPPFLAGS;
OLD_LDFLAGS=$LDFLAGS;
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
CPPFLAGS="$AM_CXXFLAGS $CPPFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"

# ensure the library to check for is covered by the LIBS variable
OLD_LIBS=$LIBS
LIBS="$LIBS -lgsl -lgslcblas -lm"
	AC_MSG_NOTICE([])
	AC_MSG_NOTICE([current CXXFLAGS = $CXXFLAGS])
	AC_MSG_NOTICE([current LDFLAGS  = $LDFLAGS])
	AC_MSG_NOTICE([current LIBS     = $LIBS])
	AC_MSG_NOTICE([])

# try to link a function of GSL library
AC_MSG_CHECKING([for the GSL C++ library headers])
AC_COMPILE_IFELSE(
     [AC_LANG_PROGRAM([[ #include <gsl/gsl_matrix_double.h> ]])],
     [ 
       AC_MSG_RESULT(yes)
       FOUND_GSL_HEADER=1;
       FOUND_GSL=1;
     ],
     [
       AC_MSG_RESULT(no)
       FOUND_GSL_HEADER=0;
       FOUND_GSL=0;
     ]
)
#AC_CHECK_LIB([m],[main])
#AC_CHECK_LIB([gslcblas],[main])
#AC_CHECK_LIB([gsl],[main])
#AC_MSG_CHECKING([whether the GNU scientific library (GSL) can be linked])
#AC_LINK_IFELSE(
#    [AC_LANG_PROGRAM([[ #include <gsl/gsl_matrix_double.h> ]],
#        [[ gsl_matrix_alloc(3, 3); ]])],
#    [AC_MSG_RESULT([yes])
#        FOUND_GSL=1;],
#    [AC_MSG_RESULT([no])
#        LIBS=$OLD_LIBS; 
#        FOUND_GSL=0;])

# reset original *FLAGS
CXXFLAGS=$OLD_CXXFLAGS
CPPFLAGS=$OLD_CPPFLAGS
LDFLAGS=$OLD_LDFLAGS


DEPENDENCYNOTFOUND=0;


# error output if BIU-headers not found
if test $BIUNOTFOUND = 1; then
	AC_MSG_NOTICE([])
	AC_MSG_NOTICE([the BIU C++ library is required for building the tools.])
	AC_MSG_NOTICE([it can be obtained from http://www.bioinf.uni-freiburg.de/SW/BIU/.])
	AC_MSG_NOTICE([])
 	if test $BIUPATHSET = 1; then
		AC_MSG_NOTICE([cannot find BIU in given path '$with_BIU'.])
 	else
		AC_MSG_NOTICE([if BIU is installed elsewhere use --with-BIU=prefix.])
 	fi
 	DEPENDENCYNOTFOUND=1;
fi

# error output if ELL-headers not found
if test $ELLNOTFOUND = 1; then
	AC_MSG_NOTICE([])
	AC_MSG_NOTICE([the ELL C++ library is required for building the tools.])
	AC_MSG_NOTICE([it can be obtained from http://www.bioinf.uni-freiburg.de/SW/ELL/.])
	AC_MSG_NOTICE([])
 	if test $ELLPATHSET = 1; then
		AC_MSG_NOTICE([cannot find ELL in given path '$with_ELL'.])
 	else
		AC_MSG_NOTICE([if ELL is installed elsewhere use --with-ELL=prefix.])
 	fi
 	DEPENDENCYNOTFOUND=1;
fi

# error output if GSL-headers not found
if test $FOUND_GSL = 0; then
	AC_MSG_NOTICE([])
	AC_MSG_NOTICE([the GNU scientific library (GSL) is required for building the tools.])
	AC_MSG_NOTICE([it can be obtained from www.gnu.org/software/gsl/])
	AC_MSG_NOTICE([and is usually available for linux distributions.])
	AC_MSG_NOTICE([])
 	if test $GSLPATHSET = 1; then
		AC_MSG_NOTICE([cannot find GSL in given path '$with_GSL'.])
 	else
		AC_MSG_NOTICE([if GSL is installed elsewhere use --with-GSL=prefix.])
 	fi
 	DEPENDENCYNOTFOUND=1;
fi


# error ABORT if one of the libraries was not found
if test $DEPENDENCYNOTFOUND = 1; then
	AC_MSG_NOTICE([])
	AC_MSG_NOTICE([current CXXFLAGS = $CXXFLAGS])
	AC_MSG_NOTICE([current LDFLAGS  = $LDFLAGS])
	AC_MSG_NOTICE([])
	AC_MSG_ERROR([library not found! Try './configure --help'.])
fi

if test "x$HAVE_PKG_CONFIG" = "xyes"; then

PKG_PROG_PKG_CONFIG
if test  $BIUPATHSET = 1 ; then
	AC_MSG_NOTICE([check if BIU version >= $BIU_REQUIRED_VERSION])
	PKG_CHECK_MODULES([BIU], [biu >= $BIU_REQUIRED_VERSION])
fi
if test  $ELLPATHSET = 1 ; then
	AC_MSG_NOTICE([check if ELL version >= $ELL_REQUIRED_VERSION])
	PKG_CHECK_MODULES([ELL], [ell >= $ELL_REQUIRED_VERSION])
fi

fi


###############################################################################
# END DEPENDENCY CHECK
###############################################################################

AC_CONFIG_FILES([ Makefile ])
AC_CONFIG_FILES([ src/version.hh ])
AC_CONFIG_FILES([ src/Makefile ])
AC_CONFIG_FILES([ doc/Makefile doc/version-info.tex ])
AC_CONFIG_FILES([ energy/Makefile ])

AC_OUTPUT
