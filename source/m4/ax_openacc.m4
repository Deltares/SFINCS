# ===========================================================================
#        https://www.gnu.org/software/autoconf-archive/ax_openacc.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_OPENACC([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use OpenACC a
#   standard API and set of compiler directives for gpu computing
#   (see https://www.openacc.org/). This macro is based on the ax_openmp macro.
#
#   On success, it sets the OPENACC_CFLAGS/OPENACC_CXXFLAGS/OPENACC_F77FLAGS
#   output variable to the flag (e.g. -acc) used both to compile *and* link
#   OpenACC programs in the current language.
#
#   NOTE: You are assumed to not only compile your program with these flags,
#   but also link it with them as well.
#
#   If you want to compile everything with OpenACC, you should set:
#
#     CFLAGS="$CFLAGS $OPENACC_CFLAGS"
#     #OR#  CXXFLAGS="$CXXFLAGS $OPENACC_CXXFLAGS"
#     #OR#  FFLAGS="$FFLAGS $OPENACC_FFLAGS"
#
#   (depending on the selected language).
#
#   The user can override the default choice by setting the corresponding
#   environment variable (e.g. OPENACC_CFLAGS).
#
#   ACTION-IF-FOUND is a list of shell commands to run if an OpenACC flag is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_OPENACC.
#
#
# LICENSE
#
#   Copyright (c) 2020 Fedor Baart <fedor.baart@deltares.nl>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 13


# _AX_LANG_OPENACC
# ---------------
# Expands to some language dependent source code for testing the presence of
# OpenACC.
AC_DEFUN([_AX_LANG_OPENACC],
[AC_LANG_SOURCE([_AC_LANG_DISPATCH([$0], _AC_LANG, $@)])])

# _AX_LANG_OPENACC(C)
# ------------------
m4_define([_AX_LANG_OPENACC(C)],
[[
#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] )
{

    // Minimal example
    // input
    double a[1];
    double b[1];
    // Output vector
    double c[1];

    a[0] = 0;
    b[0] = 1;

    #pragma acc kernels copyin(a[0],b[0]), copyout(c[0])
    c[0] = a[0] + b[0];
    printf("%f\n", c[0]);

    return 0;
}
]])

# _AX_LANG_OPENACC(C++)
# --------------------
m4_copy([_AX_LANG_OPENACC(C)], [_AX_LANG_OPENACC(C++)])

# _AX_LANG_OPENACC(Fortran 77)
# ---------------------------
m4_define([_AX_LANG_OPENACC(Fortran 77)],
[
        program main
        implicit none
        integer a, b, c
        a = 0
        b = 1
        !$acc data copy(a, b) copyout(c)
        c = a + b
        !$acc end data
        end
])

# _AX_LANG_OPENACC(Fortran)
# ------------------------
m4_copy([_AX_LANG_OPENACC(Fortran 77)], [_AX_LANG_OPENACC(Fortran)])

# AX_OPENACC
# ---------
# Check which options need to be passed to the C compiler to support OpenACC.
# Set the OPENACC_CFLAGS / OPENACC_CXXFLAGS / OPENACC_FFLAGS variable to these
# options.
# The options are necessary at compile time (so the #pragmas are understood)
# and at link time (so the appropriate library is linked with).
# This macro takes care to not produce redundant options if $CC $CFLAGS already
# supports OpenACC. It also is careful to not pass options to compilers that
# misinterpret them; for example, most compilers accept "-openacc" and create
# an output file called 'penmp' rather than activating OpenACC support.
AC_DEFUN([AX_OPENACC],
[
  OPENACC_[]_AC_LANG_PREFIX[]FLAGS=
  AC_ARG_ENABLE([openacc],
    [AS_HELP_STRING([--disable-openacc], [do not use OpenACC])])
  if test "$enable_openacc" != no; then
    AC_CACHE_CHECK([for $[]_AC_CC[] option to support OpenACC],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openacc],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openacc='unsupported'
      dnl Try these flags:
      dnl   GCC >= 4.2           -fopenacc
      dnl   pg, nvfortran            -acc
      dnl If in this loop a compiler is passed an option that it doesn't
      dnl understand or that it misinterprets, the AC_LINK_IFELSE test
      dnl will fail (since we know that it failed without the option),
      dnl therefore the loop will continue searching for an option, and
      dnl no output file called 'penmp' or 'mp' is created.
      for ac_option in -fopenacc -acc; do
        ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
        _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
        AC_LINK_IFELSE([_AX_LANG_OPENACC],
          [ac_cv_prog_[]_AC_LANG_ABBREV[]_openacc=$ac_option])
        _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
        if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openacc" != unsupported; then
          break
        fi
      done])
    case $ac_cv_prog_[]_AC_LANG_ABBREV[]_openacc in #(
      "none needed" | unsupported)
        ;; #(
      *)
        OPENACC_[]_AC_LANG_PREFIX[]FLAGS=$ac_cv_prog_[]_AC_LANG_ABBREV[]_openacc ;;
    esac
  fi
  AC_SUBST([OPENACC_]_AC_LANG_PREFIX[FLAGS])
])
