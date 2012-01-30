AC_PREREQ([2.61])
AC_INIT([pepStat], [0.0.1], [http://www.rglab.org])

AC_CONFIG_SRCDIR([DESCRIPTION])

# determine RHOME
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

## pick up compiler as will be used by R CMD INSTALL/SHLIB
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP


## BLAS.
AC_ARG_WITH([blas],
[AC_HELP_STRING([--with-blas],
                [use system gsl BLAS library (if available)])],
[BLAS_LIBS=$with_blas],
[BLAS_LIBS="-lgslcblas"])

PKG_CHECK_MODULES(GSL, [gsl >= 1.2], [HAVEGSL="yes"], [HAVEGSL="no"])
if test "$HAVEGSL" = "yes"; then
  LIBS=$GSL_LIBS
else 
  AM_PATH_GSL(1.2,[LIBS="${LIBS} -lgsl ${BLAS_LIBS} -lm"],[AC_MSG_ERROR([Cannot find Gnu Scientific Library >=1.6])])
fi

AC_SUBST(OPENMP_CFLAGS)
# Checks for header files.
AC_CHECK_HEADERS(dispatch/dispatch.h,
  :,
  [
    AC_CHECK_HEADERS(omp.h,
      [
        LIBS+=-lgomp
        OPENMP_CFLAGS=-fopenmp
      ],
      [
        AC_MSG_NOTICE([no parallelism])
      ])
  ]
)

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT()
