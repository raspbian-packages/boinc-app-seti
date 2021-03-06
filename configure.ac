#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id: $


AC_PREREQ(2.57)
AC_INIT(setiathome_v8, 8.00, korpela@ssl.berkeley.edu)
AC_CONFIG_MACRO_DIR([m4])
svnrev="`svn info | grep Revision:`"

if test -n "${svnrev}" ; then
  AC_REVISION("${svnrev}")
  AC_DEFINE_UNQUOTED(SVN_REV,"${svnrev}",
    [Define to be the subversion revision number])
else
  AC_REVISION([$Revision: 1146 $])
  AC_DEFINE(SVN_REV,"$Revision: 1146 $",
    [Define to be the subversion revision number])
fi

echo "WARNING! Changes have been made which are not compatible with SETI@home"
echo "v7.  Do not use this version with the projects until v8 is released"

AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE(dist-zip)
AC_CONFIG_SRCDIR([client/analyze.h])

AC_CACHE_SAVE

# Set up the prefix early.
if test "x$prefix" = "xNONE" ; then
  prefix="${ac_default_prefix}"
fi
if test "x$exec_prefix" = "xNONE" ; then
  exec_prefix="${ac_default_prefix}"
fi


# Make sure we use the same initial CFLAGS and CXXFLAGS unless otherwise told.
if test -z "${CXXFLAGS}"
then
  CXXFLAGS="${CFLAGS}"
elif test -z "${CFLAGS}"
then
  CFLAGS="${CXXFLAGS}"
fi

# make sure $includedir is in CFLAGS and $libdir is in LDFLAGS.
tmpinc="`eval echo ${includedir}`"
tmplib="`eval echo ${libdir}`"
if test -z "`echo ${CFLAGS} | grep ${tmpinc}`" ; then
  CFLAGS="${CFLAGS} -I${tmpinc}"
fi
if test -z "`echo ${CXXFLAGS} | grep ${tmpinc}`" ; then
  CXXFLAGS="${CXXFLAGS} -I${tmpinc}"
fi
if test -z "`echo ${LDFLAGS} | grep ${tmplib}`" ; then
  LDFLAGS="${LDFLAGS} -L${tmplib}"
fi

AC_ARG_ENABLE(intrinsics,
    AS_HELP_STRING([--disable-intrinsics],
                   [disable use of intrinsics in SIMD code]))


if test "x$enable_intrinsics" != xno ; then
  enable_intrinsics=yes
fi

AC_ARG_ENABLE(graphics,
    AC_HELP_STRING([--disable-graphics],
                   [disable building the client graphics]))

if test "x$enable_graphics" != xno ; then
  enable_graphics=yes
else
  enable_graphics=no
fi


AC_ARG_ENABLE(tests,
    AC_HELP_STRING([--disable-tests],
                   [disable tests (fakedata and hires_timer_test)]),
    [],
    [enable_tests=yes])

AC_ARG_ENABLE(server,
    AC_HELP_STRING([--disable-server],
                   [disable building the seti@home server components]),
    [],
    [enable_server=yes])


AC_ARG_ENABLE(client,
    AC_HELP_STRING([--disable-client],
                   [disable building the client]),
    [],
    [enable_client=yes])

AC_ARG_ENABLE(static-client,
    AC_HELP_STRING([--disable-static-client],
                   [allow dynamic libraries to be used in the client]),
    [enable_static_client=${enableval}],
    [enable_static_client=yes])

m4_divert_once([HELP_ENABLE],
  AC_HELP_STRING([], [Default: --enable-server --enable-client:
                      builds both server and client]))

if test "${enable_server}" = yes ; then
    if test "${enable_client}" = yes ; then
        configured_to_build='server & client'
    else
        configured_to_build='server only'
    fi
else
    if test "${enable_client}" = yes ; then
        configured_to_build='client only'
    else
        AC_MSG_WARN([Youve disabled both the server and the client. Nothing useful will be built])
    fi
fi

AC_MSG_NOTICE(["--- Configuring SETI_BOINC AC_PACKAGE_VERSION (${configured_to_build}) ---"])


AC_COPYRIGHT([
Copyright (C) 2004 Regents of the University of California

SETI_BOINC is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2, or (at your option) any later
version.

SETI_BOINC is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with SETI_BOINC; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

Contributor(s):
])

major_version=`echo AC_PACKAGE_VERSION | sed 's/\..*//'`
minor_version=`echo AC_PACKAGE_VERSION | sed 's/.*\.//' | sed 's/^0//'`
AC_SUBST([MAJOR_VERSION], [$major_version])
AC_SUBST([MINOR_VERSION], [$minor_version])
AC_DEFINE_UNQUOTED([VERSION_MAJOR],$major_version, [SETI@home major version number])
AC_DEFINE_UNQUOTED([VERSION_MINOR],$minor_version, [SETI@home minor version number])
AC_DEFINE_UNQUOTED([SAH_APP_NAME],["$PACKAGE_NAME"], [Define to the BOINC application name for setiathome])
SAH_TOP_DIR=`pwd`
AC_SUBST([SAH_TOP_DIR])

AM_MAINTAINER_MODE

# Checks for programs.
PATH="/usr/xpg4/bin:${PATH}"
AC_PROG_CXX
AC_PROG_CXXCPP
AM_PROG_CC_C_O
AM_PROG_AS
AC_PROG_AWK
AC_PROG_LN_S
AC_PATH_PROG(TR,[tr])
AC_PATH_PROGS(AR,[ar lib])
AC_PATH_PROG(GREP,[grep])
AC_PATH_PROG(AUTOCONF,[autoconf])
AC_PATH_PROG(AUTOHEADER,[autoheader])
AC_PATH_PROGS(INDENT,[astyle indent])
AC_PATH_PROGS(CP,[cp copy])
AC_PATH_PROGS(LN,[ln cp copy])
AC_PATH_PROG(SORT,[sort])
AC_PATH_PROG(UNIQ,[uniq])
AC_PATH_PROG(CAT,[cat type])
AC_PATH_PROG(MV,[mv])
AC_PATH_PROGS(RM,[rm Rm del erase delete])

AC_LIBTOOL_DLOPEN
AC_LIBTOOL_WIN32_DLL
AM_PROG_LIBTOOL
AC_SUBST(PICFLAGS,${lt_prog_compiler_pic})
AC_SUBST(SED)
AC_C_BIGENDIAN
AX_C_FLOAT_WORDS_BIGENDIAN

SAH_OPTION_BITNESS

SAH_REQUIRES([boinc],[SAH_CHECK_BOINC],["${no_boinc}" = yes],
[
ERROR: trying to build the seti_boinc client or server but BOINC was not
found.
You can get boinc at http://boinc.ssl.berkeley.edu/
],
[exit 1])
AC_SUBST([BOINCDIR])

dnl Determine the BOINC platform given the target arch-platform-os.
BOINC_PLATFORM

if test -n `echo $INDENT | grep astyle`
then
  AC_SUBST([INDENT_FLAGS],["--c --indent-classes --indent-switches --brackets=attach --convert-tabs"])
else
  AC_SUBST([INDENT_FLAGS],["-kr"])
fi
AC_PROG_MAKE_SET
SAH_DLLEXT
SAH_LIBEXT
if test -n "$EXEEXT"
then
  DOTEXEEXT=".$EXEEXT"
fi
AC_SUBST(DOTEXEEXT)
AC_SYS_LARGEFILE
AX_PTHREAD([AC_DEFINE(HAVE_PTHREAD,1, [Have pthread])])

if test "x$enable_client" = xyes ; then

# find the right flag for building a static app
if test "x$enable_static_client" != xno ; then
  svldflags="$LDFLAGS"
  SAH_CHECK_LDFLAG([-static],[LDSTATIC=-static])
  LDFLAGS="$svldflags"
fi

# Check libraries required or desired for seti_boinc app
# Checks for libraries.
APP_LDFLAGS="$LDSTATIC $LDFLAGS"
APP_CFLAGS="$CFLAGS"
AC_CHECK_LIB([log],[__android_log_print],
    [APP_LIBS="-llog ${APP_LIBS}"])

dnl AC_CHECK_LIB([corkscrew],[unwind_backtrace_signal_arch],
dnl    [APP_LIBS="-lcorkscrew ${APP_LIBS}"])

AC_CHECK_LIB([m], [sin],[
    AC_DEFINE([HAVE_LIBM],[1],[Define to 1 if you have the math library])
    APP_LIBS="-lm ${APP_LIBS}"
])
SAH_CHECK_MATH_FUNCS([sqrt floor sinf cosf sincos sincosf atanf isnan _isnan __isnan __isnanf _isnanf isnanf])
AC_CHECK_LIB([fftw3f], [fftwf_plan_dft_1d], [
    AC_DEFINE([HAVE_LIBFFTW3F],[1],[Define to 1 if you have the fftw3f library])
    APP_LIBS="-lfftw3f ${APP_LIBS}"
],
[],
[-lm])
AC_CHECK_LIB([fftw3], [fftwf_plan_dft_1d], [
    AC_DEFINE([HAVE_LIBFFTW3F],[1],[Define to 1 if you have the fftw3f library])
    APP_LIBS="-lfftw3 ${APP_LIBS}"
],
[],
[-lm])
AC_CHECK_LIB([cpufeatures], [android_getCpuFeatures], [
    AC_DEFINE([HAVE_LIBCPUFEATURES],[1],[Define to 1 if you have the cpufeatures library])
    APP_LIBS="-lcpufeatures ${APP_LIBS}"
],
[],
[-lm])
AC_ARG_ENABLE(neon,
    AC_HELP_STRING([--disable-neon],
                   [Disable neon instruction set on ARM devices]),
    [enable_neon=${enableval}],
    [enable_neon=yes])
if test "x${enable_neon}" = xyes ; then
  AC_DEFINE([USE_NEON],1,[Define to 1 to enable neon instructions on ARM devices])
fi

SAH_CHECK_ASMLIB
APP_LIBS="${ASMLIB_LIBS} ${APP_LIBS}"
APP_CFLAGS="${ASMLIB_CFLAGS} ${APP_CFLAGS}"
APP_LDFLAGS="${ASMLIB_LDFLAGS} ${APP_LDFLAGS}"

# Check for additional optimization flags
SAH_OPTIMIZATIONS

# Check for intel performance primatives
SAH_CHECK_IPP
if test x$found_ipp = xyes ; then
  APP_LDFLAGS="-L${IPPDIR}/lib ${APP_LDFLAGS}"
  APP_LIBS="-lippcore -lippsmerged ${APP_LIBS}"
  APP_CFLAGS="-I${IPPDIR}/include ${APP_CFLAGS}"
fi

AC_SUBST([APP_LIBS])
AC_SUBST([APP_LDFLAGS])
AC_SUBST([APP_CFLAGS])

if test "$enable_graphics" = yes ; then
  SAH_GRX_LIBS
  AX_CHECK_GL
  AX_CHECK_GLU
  AX_CHECK_GLUT
  if test "$no_gl" = yes -o "$no_glu" = yes -o "$no_glut" = yes ; then
     have_glut=no
     enable_graphics=no
     AC_MSG_WARN([
================================================================================
WARNING: GL/GLU/GLUT not found.

The GL, GLU, and GLUT libraries are required in order to build the graphical
parts of the BOINC application API library.

==> only building non-graphical parts of the BOINC API Library for now.

================================================================================
       ])
else
     AC_DEFINE([BOINC_APP_GRAPHICS],1,[Define to 1 to build a graphical application])
     SAH_GRX_INCLUDES
     GRAPHICS_CFLAGS="${GRX_CFLAGS} ${GL_CFLAGS} ${GLU_CFLAGS} ${GLUT_CFLAGS}"
     GRAPHICS_LIBS="${GL_LIBS} ${GLU_LIBS} ${GLUT_LIBS} ${GRX_LIBS} -ljpeg"

case ${target} in
        powerpc-apple*)       GRAPHICS_LIBS="${GRAPHICS_LIBS} -framework AppKit -framework CoreServices"
        ;;
esac
      GRAPHICS_LIBS_RAW="${PTHREAD_CFLAGS} -L${BOINCDIR}/api -L${BOINCDIR}/lib -lboinc_graphics2 -lboinc ${GRAPHICS_LIBS}"
dnl     SAH_RAW_LDFLAGS(["${PTHREAD_CFLAGS} -L${BOINCDIR}/api -L${BOINCDIR}/lib -lboinc_graphics2 -lboinc ${GRAPHICS_LIBS}"],[GRAPHICS_LIBS_RAW])
     have_glut=yes
     enable_graphics=yes
  fi
fi
AC_SUBST([LDSTATIC])
AC_SUBST([GRAPHICS_CFLAGS])
AC_SUBST([GRAPHICS_LIBS])
AC_SUBST([GRAPHICS_LIBS_RAW])

AC_CACHE_SAVE

fi

dnl SAH_CHECK_LIB([c],[atexit],
dnl     AC_DEFINE([HAVE_LIBC],[1],[Define to 1 if you have the c library]))
dnl SAH_CHECK_LIB([pthread],[pthread_join],
dnl     AC_DEFINE([HAVE_LIBPTHREAD],[1],[Define to 1 if you have the pthread library]))
dnl AC_CHECK_LIB([gcc],[sscanf],
dnl     AC_DEFINE([HAVE_LIBGCC],[1],[Define to 1 if you have the gcc library]))
dnl SAH_CHECK_LIB([gcc_eh],[_Unwind_Resume],
dnl     AC_DEFINE([HAVE_LIBGCC_EH],[1],[Define to 1 if you have the gcc_eh library]))
dnl SAH_CHECK_LIB([gcc_s],[sscanf],
dnl     AC_DEFINE([HAVE_LIBGCC_S],[1],[Define to 1 if you have the gcc_s library]))
dnl AC_CHECK_LIB([stdc++],[sscanf],
dnl     AC_DEFINE([HAVE_LIBSTDC__],[1],[Define to 1 if you have the stdc++ library]))
dnl AC_CHECK_LIB([z], [uncompress],
dnl     AC_DEFINE([HAVE_LIBZ],[1],[Define to 1 if you have the z library]))
dnl AC_CHECK_LIB([socket], [bind],
dnl     AC_DEFINE([HAVE_LIBSOCKET],[1],[Define to 1 if you have the socket library]))
dnl AC_CHECK_LIB([rt], [sched_get_priority_min],
dnl     AC_DEFINE([HAVE_LIBRT],[1],[Define to 1 if you have the rt library]))
dnl AC_CHECK_LIB([nsl], [gethostbyname],
dnl     AC_DEFINE([HAVE_LIBNSL],[1],[Define to 1 if you have the nsl library]))
dnl AC_CHECK_LIB([elf], [elf_hash],
dnl     AC_DEFINE([HAVE_LIBELF],[1],[Define to 1 if you have the elf library]))
dnl AC_CHECK_LIB([aio], [aio_fork],
dnl     AC_DEFINE([HAVE_LIBAIO],[1],[Define to 1 if you have the aio library]))
CHECK_SSL
AC_CHECK_LIB([ssl], [fopen],
    [BOINC_EXTRA_LIBS="${BOINC_EXTRA_LIBS} -L${SSLDIR}/lib ${sah_lib_last}"])
AC_CHECK_LIB([crypto], [RSA_new],
    [BOINC_EXTRA_LIBS="${BOINC_EXTRA_LIBS} -L${SSLDIR}/lib ${sah_lib_last}"])

dnl AC_CHECK_LIB([dl], [dlopen],
dnl   AC_DEFINE([HAVE_LIBDL],[1],[Define to 1 if you have the dl library]))
AC_LANG(C++)

case ${target} in
        *-apple-darwin*)       CLIENT_EXTRA_LIBS="${CLIENT_EXTRA_LIBS} -framework Accelerate -framework CoreServices"
        ;;
esac


if test "${ac_cv_cxx_compiler_gnu}" = "yes" ; then
  gcc_version_string="`${CXX} --version | head -1`"
  AC_DEFINE_UNQUOTED([COMPILER_STRING],"$gcc_version_string",[Define to a string identifying your compiler])
else if test -n "${CXX}" ; then
  AC_DEFINE_UNQUOTED([COMPILER_STRING],"$CXX",[Define to a string identifying your compiler])
else if test -n "${CC}" ; then
  AC_DEFINE_UNQUOTED([COMPILER_STRING],"$CC",[Define to a string identifying your compiler])
fi
fi
fi



SAH_SERVER_REQUIRES([mysql],[SAH_CHECK_MYSQL],["${no_mysql}" = yes])
SAH_SERVER_REQUIRES([informix],[SAH_CHECK_INFORMIX],["${no_informix}" = yes])
SAH_SERVER_REQUIRES([setilib],[SAH_FIND_SETILIB],["${no_setilib}" = yes])
SAH_SERVER_REQUIRES([healpix],[SAH_CHECK_HEALPIX],["${no_healpix}" = yes])
SAH_SERVER_REQUIRES([cfitsio],[SAH_CHECK_CFITSIO],["${no_cfitsio}" = yes])
SAH_SERVER_REQUIRES([seti_gbt],[SAH_CHECK_SETI_GBT],["${no_seti_gbt}" = yes])

# Checks for header files.
AC_HEADER_STDC
SAH_LARGEFILE_BREAKS_CXX
SAH_HEADER_STDCXX
AC_HEADER_SYS_WAIT
AC_HEADER_TIME
AC_CHECK_HEADERS([mmap.h mman.h sys/mman.h fcntl.h inttypes.h limits.h memory.h alloca.h malloc.h stdlib.h string.h strings.h errno.h sys/time.h sys/types.h sys/ioctl.h sys/statvfs.h unistd.h dirent.h math.h float.h ieeefp.h floatingpoint.h complex.h fftw3.h setjmp.h signal.h mach/mach_time.h sys/param.h sys/systm.h machine/cpu.h cpu-features.h sys/cpu-features.h machine/cpu-features.h])
save_cxxflags="${CXXFLAGS}"
save_cppflags="${CPPFLAGS}"
CXXFLAGS="${CXXFLAGS} -msse3"
CPPFLAGS="${CPPFLAGS} -msse3"
AC_CHECK_HEADERS([intrin.h x86intrin.h pmmintrin.h xmmintrin.h emmintrin.h])
CXXFLAGS="${save_cxxflags} -mavx"
CPPFLAGS="${save_cppflags} -mavx"
AC_CHECK_HEADERS([immintrin.h avxintrin.h])
CXXFLAGS="${save_cxxflags}"
CPPFLAGS="${save_cppflags}"
if test "x$enable_intrinsics" != "xno" ; then
  if test "x$ac_cv_header_xmmintrin_h" != "xno" -o "x$ac_cv_header_emmintrin_h" != "xno" ; then
    AC_DEFINE([USE_INTRINSICS],[1],[Define to 1 to use SIMD intrinsics rather than inline assembly])
  fi
fi
AC_CACHE_SAVE

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_C_LONG_DOUBLE
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_STRUCT_ST_BLOCKS
AC_STRUCT_TM
AC_CHECK_SIZEOF([short])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long int])
AC_CHECK_SIZEOF([float])
AC_CHECK_SIZEOF([double])
AC_CHECK_SIZEOF([long double])
AC_CHECK_TYPES([long long,_int64,int64_t,bool,_int32, _uint32, int32_t, uint32_t, hrtime_t, uint_fast64_t, uint64_t, _uint64, u_int64_t, ptrdiff_t, ssize_t, off64_t])
AC_CHECK_TYPES([sigjmp_buf, jmp_buf], [], [], [[#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
# include <sys/stat.h>
#endif
#ifdef STDC_HEADERS
# include <stdlib.h>
# include <stddef.h>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# endif
#endif
#ifdef HAVE_STRING_H
# if !defined STDC_HEADERS && defined HAVE_MEMORY_H
#  include <memory.h>
# endif
# include <string.h>
#endif
#ifdef HAVE_STRINGS_H
# include <strings.h>
#endif
#ifdef HAVE_INTTYPES_H
# include <inttypes.h>
#endif
#ifdef HAVE_STDINT_H
# include <stdint.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif
#ifdef HAVE_SIGNAL_H
# include <signal.h>
#endif
#ifdef HAVE_SETJMP_H
# include <setjmp.h>
#endif
]])
CXXFLAGS="${save_cxxflags} -mavx"
CPPFLAGS="${save_cppflags} -mavx"
AC_CHECK_TYPES([__m256,_m256,__m256d,_m256d,__m256i,_m256i],[],[],[[#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
# include <sys/stat.h>
#endif
#ifdef STDC_HEADERS
# include <stdlib.h>
# include <stddef.h>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# endif
#endif
#ifdef HAVE_STRING_H
# if !defined STDC_HEADERS && defined HAVE_MEMORY_H
#  include <memory.h>
# endif
# include <string.h>
#endif
#ifdef HAVE_STRINGS_H
# include <strings.h>
#endif
#ifdef HAVE_INTTYPES_H
# include <inttypes.h>
#endif
#ifdef HAVE_STDINT_H
# include <stdint.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif
#define __AVX__ 1
#if defined(HAVE_IMMINTRIN_H)
#include <immintrin.h>
#elif defined(HAVE_AVXINTRIN_H)
#include <avxintrin.h>
#elif defined(HAVE_X86INTRIN_H)
#include <x86intrin.h>
#elif defined(HAVE_INTRIN_H)
#include <intrin.h>
#endif
]])
CXXFLAGS="${save_cxxflags}"
CPPFLAGS="${save_cppflags}"

AC_CACHE_SAVE

# Checks for library functions.
AC_LANG(C)
AC_FUNC_FORK
AC_HEADER_MAJOR
AC_FUNC_ALLOCA
if test "${ac_cv_func_alloca_works}" = "yes" ; then
  ac_cv_func_alloca="yes"
fi
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_CHECK_DECLS([__builtin_alloca,__builtin_popcount,__builtin_prefetch,__builtin_clz,__builtin_ffs,__builtin_isnan],
   [],[],[[
#include <stdio.h>
#include <stdlib.h>
]])
if test "${ac_cv_have_decl___builtin_alloca}" = "yes" ; then
  ac_cv_func___builtin_alloca="yes"
fi
if test "${ac_cv_have_decl___builtin_popcount}" = "yes" ; then
  ac_cv_func___builtin_popcount="yes"
fi
if test "${ac_cv_have_decl___builtin_prefetch}" = "yes" ; then
  ac_cv_func___builtin_prefetch="yes"
fi
if test "${ac_cv_have_decl___builtin_clz}" = "yes" ; then
  ac_cv_func___builtin_clz="yes"
fi
if test "${ac_cv_have_decl___builtin_ffs}" = "yes" ; then
  ac_cv_func___builtin_ffs="yes"
fi
if test "${ac_cv_have_decl___builtin_isnan}" = "yes" ; then
  ac_cv_func___builtin_isnan="yes"
fi

AC_CHECK_FUNCS([ _aligned_malloc alloca _alloca __builtin_alloca __builtin_popcount __builtin_prefetch __builtin_clz __builtin_ffs posix_memalign memalign atexit exit _exit getcwd memset munmap putenv strchr strstr strcasestr atoll dlopen gethrtime mach_absolute_time get_cyclecount clock_gettime nanotime microtime gettimeofday sqrt floor sinf cosf sincos sincosf atanf isnan _isnan __isnan __builtin_isnan isnanf _isnanf __isnanf siglongjmp sigsetjmp sigaction sysv_signal bsd_signal ])
SAH_CHECK_NAMESPACES
AH_TEMPLATE([HAVE_STD_MIN],[Define to 1 if min is in namespace std::])
AH_TEMPLATE([HAVE_STD_MAX],[Define to 1 if max is in namespace std::])
AH_TEMPLATE([HAVE_STD_TRANSFORM],[Define to 1 if transform is in namespace std::])
SAH_FUNCS_IN_NAMESPACE([['min(0,0)'] ['max(0,0)'] ['transform((char *)0,(char *)0,(char *)0,(int(*)(int))malloc)']],std)

#SAH_AVX
have_avx=no

AM_CONDITIONAL(USE_MY_LIBJPEG, [false])
AM_CONDITIONAL(USE_MY_IMAGELIBS, [false])
AM_CONDITIONAL(USE_MY_GLUT, [false])
AM_CONDITIONAL(ENABLE_TESTS, [test "${enable_tests}" = yes])
AM_CONDITIONAL(ENABLE_SERVER, [test "${enable_server}" = yes])
AM_CONDITIONAL(ENABLE_CLIENT, [test "${enable_client}" = yes])
AM_CONDITIONAL(ENABLE_GUI, [test "${enable_graphics}" = yes])
AM_CONDITIONAL(LINUX, [test -n "`echo ${target} | grep linux`"])
AM_CONDITIONAL(I386, [test -n "`echo ${target} | $EGREP 'i.86|x86_64|amd64'`" -a $COMPILER_MODEL_BITS -eq 32])
AM_CONDITIONAL(X86_64, [test -n "`echo ${target} | $EGREP 'i.86|x86_64|amd64'`" -a $COMPILER_MODEL_BITS -eq 64 ])
AM_CONDITIONAL(PPC, [test -n "`echo ${target} | grep 'powerpc'`" -o -n "`echo ${target} | grep 'ppc'`" ])
AM_CONDITIONAL(AVX, [test "x${have_avx}" = "xyes"])
if test "x${avx_type}" != "x" ; then
  AC_DEFINE_UNQUOTED([AVX_MASKSTORE_TYPECAST(x)],[reinterpret_cast<$avx_type>(x)],
    [Define to the typecast required for arg 2 of _mm256_maskstore_ps()])
else
  AC_DEFINE_UNQUOTED([AVX_MASKSTORE_TYPECAST(x)],[(x)],
    [Define to the typecast required for arg 2 of _mm256_maskstore_ps()])
fi

AC_ARG_WITH([float-abi],
    AC_HELP_STRING([--with-float-abi=<abi>],
                   [set -mfloat-abi parameter for gcc]),
    [FP_ABI_FLAGS="-mfloat-abi=$withval"],[
    if test "x$FP_ABI_FLAGS" = "x" && test "x`echo ${CFLAGS} ${CPPFLAGS} ${CXXFLAGS} | grep mfloat-ab`" = "x" ; then
        case $host in
            arm-*-armhf|arm-*-gnueabihf)
                    FP_ABI_FLAGS="-mfloat-abi=hard"
                    ;;
            arm-*-armel|arm-*-androideabi|arm-*-gnueabi)
                    FP_ABI_FLAGS="-mfloat-abi=softfp"
                    ;;
            *)      FP_ABI_FLAGS=
                    ;;
        esac
    fi
])

AC_ARG_WITH([fpu],
    AC_HELP_STRING([--with-fpu=<fpu>],
                   [set -mfpu parameter for gcc]),
    [FP_FLAGS="-mfpu=$withval"],[
    if test "x$FP_ABI_FLAGS" != "x-mfloat-abi=soft" && ! ( echo ${CFLAGS} ${CPPFLAGS} ${CXXFLAGS} | grep "mfloat-abi=soft[[ $]]" >/dev/null ) ; then
        if test "x$FP_FLAGS" = "x" && ! ( echo ${CFLAGS} ${CPPFLAGS} ${CXXFLAGS} | grep mfpu >/dev/null  ) ; then
            case $host in
                arm-*-armhf|arm-*-armel|arm-*-androideabi)
                        FP_FLAGS="-mfpu=vfp"
                        ;;
                *)      FP_FLAGS=
                        ;;
            esac
        fi
    fi
])


AC_SUBST([FP_ABI_FLAGS])
AC_SUBST([FP_FLAGS])

AC_CACHE_SAVE

AH_TOP([
#ifndef _SAH_CONFIG_H_
#define _SAH_CONFIG_H_

#if defined(__linux__)
#define _POSIX_C_SOURCE 1
#endif

#ifdef _WIN32
#include "win-sah_config.h"
#else
])


AH_BOTTOM([

/* Define USE_NAMESPACES if you may access more than one database from the
 * same program
 */

#endif

/*
 * Use fftw if we have the library
 */
#if defined(HAVE_LIBFFTW3F) && defined(HAVE_FFTW3_H)
#define USE_FFTWF
#endif

#if defined(USE_INFORMIX) && defined(USE_MYSQL) && defined(HAVE_NAMESPACES)
#define USE_NAMESPACES
#endif

#if !defined(CUSTOM_STRING) && defined(COMPILER_STRING)
#define CUSTOM_STRING PACKAGE_STRING" "SVN_REV" "COMPILER_STRING
#endif

#include "std_fixes.h"

#endif
])

AC_SUBST([CLIENT_EXTRA_LIBS])

AC_CONFIG_FILES([Makefile
                 image_libs/Makefile
                 client/Makefile
                 tools/Makefile
                 db/Makefile
                 db/tools/Makefile
                 db/schema_to_class
                 assimilator/Makefile
                 validate/Makefile
                 splitter_pfb/Makefile
                 splitter_fft/Makefile
                ])
AC_CONFIG_HEADERS([sah_config.h])

AC_OUTPUT
chmod +x db/schema_to_class
