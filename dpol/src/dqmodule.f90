module dqmodule

!  DQFUN: A thread-safe double-quad precision computation package

!  This software requires IEEE 128-bit floating-point arithmetic, in hardware or
!  software (it is currently provided, for instance, by the gfortran compiler).

!  High-level language interface module (DQMODULE).

!  Revision date:  5 Sep 2019

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2017 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to double-quad precision (approximately 66 digits), by making only
!    relatively minor changes to existing Fortran-90 programs.  All basic arithmetic
!    operations and transcendental functions are supported, together with several
!    special functions.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.  There are NO global shared variables (except
!    static compile-time data), and NO initialization is necessary.

!  DOCUMENTATION:
!    Full documentation is not yet available, but will be provided.  In the
!    meantime, see the brief summary in README.txt in the main DQFUN directory.

!  DESCRIPTION OF THIS MODULE (DQMODULE):
!    This module contains all high-level Fortran-90 language interfaces.

use dqfunmod
implicit none
type dq_real
  sequence
  real (dqknd) dqr(2)
end type
type dq_complex
  sequence
  real (dqknd) dqc(4)
end type
type (dq_real), public:: dqeps

private &
  dq_eqqq, dq_eqqz, dq_eqdq, dq_eqqd, dq_eqxq, dq_eqqx, &
  dq_addqq, dq_addqz, dq_adddq, dq_addqd, dq_addxq, dq_addqx, &
  dq_subqq, dq_subqz, dq_subdq, dq_subqd, dq_subxq, dq_subqx, dq_negq, &
  dq_mulqq, dq_mulqz, dq_muldq, dq_mulqd, dq_mulxq, dq_mulqx, &
  dq_divqq, dq_divqz, dq_divdq, dq_divqd, dq_divxq, dq_divqx, &
  dq_expqq, dq_expqi, dq_expdq, dq_expqd, &
  dq_eqtqq, dq_eqtqz, dq_eqtdq, dq_eqtqd, dq_eqtxq, dq_eqtqx, &
  dq_netqq, dq_netqz, dq_netdq, dq_netqd, dq_netxq, dq_netqx, &
  dq_letqq, dq_letdq, dq_letqd, dq_getqq, dq_getdq, dq_getqd, &
  dq_lttqq, dq_lttdq, dq_lttqd, dq_gttqq, dq_gttdq, dq_gttqd

private &
  dq_eqzq, dq_eqzz, dq_eqdz, dq_eqzd, dq_eqxz, dq_eqzx, &
  dq_addzq, dq_addzz, dq_adddz, dq_addzd, dq_addxz, dq_addzx, &
  dq_subzq, dq_subzz, dq_subdz, dq_subzd, dq_subxz, dq_subzx, dq_negz, &
  dq_mulzq, dq_mulzz, dq_muldz, dq_mulzd, dq_mulxz, dq_mulzx, &
  dq_divzq, dq_divzz, dq_divdz, dq_divzd, dq_divxz, dq_divzx, &
  dq_expzi, dq_expzz, dq_exprz, dq_expzr, &
  dq_eqtzq, dq_eqtzz, dq_eqtdz, dq_eqtzd, dq_eqtxz, dq_eqtzx, &
  dq_netzq, dq_netzz, dq_netdz, dq_netzd, dq_netxz, dq_netzx

private &
  dq_absq, dq_absz, dq_acos, dq_imag, dq_aint, dq_anint, dq_asin, dq_atan, &
  dq_atan2, dq_conjg, dq_cos, dq_cosz, dq_cosh, dq_qtod, dq_ztod, dq_qtox, &
  dq_ztox, dq_exp, dq_expz, dq_log, dq_logz, dq_log10, dq_maxq, dq_maxq3, &
  dq_minq, dq_minq3, dq_modq, dq_qtoz, dq_dtoz, dq_xtoz, dq_qqtoz, &
  dq_ddtoz, dq_cssh, dq_cssn, dq_nrt, dq_poly, dq_inpq, dq_qtoq, dq_ztoz, &
  dq_inpz, dq_ztoq, dq_dtoq, dq_xtoq, dq_atoq, dq_itoq, dq_outq, dq_outz, &
  dq_log2, dq_pi, dq_signq, dq_sin, dq_sinz, dq_sinh, dq_sqrtq, dq_sqrtz, &
  dq_tan, dq_tanh

!  Operator extension interface blocks.

interface assignment (=)
  module procedure dq_eqqq
  module procedure dq_eqqz
  module procedure dq_eqdq
  module procedure dq_eqqd
  module procedure dq_eqxq
  module procedure dq_eqqx
  module procedure dq_eqqa
  module procedure dq_eqzq
  module procedure dq_eqzz
  module procedure dq_eqdz
  module procedure dq_eqzd
  module procedure dq_eqxz
  module procedure dq_eqzx
end interface

interface operator (+)
  module procedure dq_addqq
  module procedure dq_addqz
  module procedure dq_adddq
  module procedure dq_addqd
  module procedure dq_addxq
  module procedure dq_addqx
  module procedure dq_addzq
  module procedure dq_addzz
  module procedure dq_adddz
  module procedure dq_addzd
  module procedure dq_addxz
  module procedure dq_addzx
end interface

interface operator (-)
  module procedure dq_subqq
  module procedure dq_subqz
  module procedure dq_subdq
  module procedure dq_subqd
  module procedure dq_subxq
  module procedure dq_subqx
  module procedure dq_subzq
  module procedure dq_subzz
  module procedure dq_subdz
  module procedure dq_subzd
  module procedure dq_subxz
  module procedure dq_subzx

  module procedure dq_negq
  module procedure dq_negz
end interface

interface operator (*)
  module procedure dq_mulqq
  module procedure dq_mulqz
  module procedure dq_muldq
  module procedure dq_mulqd
  module procedure dq_mulxq
  module procedure dq_mulqx
  module procedure dq_mulzq
  module procedure dq_mulzz
  module procedure dq_muldz
  module procedure dq_mulzd
  module procedure dq_mulxz
  module procedure dq_mulzx
end interface

interface operator (/)
  module procedure dq_divqq
  module procedure dq_divqz
  module procedure dq_divdq
  module procedure dq_divqd
  module procedure dq_divxq
  module procedure dq_divqx
  module procedure dq_divzq
  module procedure dq_divzz
  module procedure dq_divdz
  module procedure dq_divzd
  module procedure dq_divxz
  module procedure dq_divzx
end interface

interface operator (**)
  module procedure dq_expqq
  module procedure dq_expqi
  module procedure dq_expdq
  module procedure dq_expqd
  module procedure dq_expzi
  module procedure dq_expzz
  module procedure dq_exprz
  module procedure dq_expzr
end interface

interface operator (.eq.)
  module procedure dq_eqtqq
  module procedure dq_eqtqz
  module procedure dq_eqtdq
  module procedure dq_eqtqd
  module procedure dq_eqtxq
  module procedure dq_eqtqx
  module procedure dq_eqtzq
  module procedure dq_eqtzz
  module procedure dq_eqtdz
  module procedure dq_eqtzd
  module procedure dq_eqtxz
  module procedure dq_eqtzx
end interface

interface operator (.ne.)
  module procedure dq_netqq
  module procedure dq_netqz
  module procedure dq_netdq
  module procedure dq_netqd
  module procedure dq_netxq
  module procedure dq_netqx
  module procedure dq_netzq
  module procedure dq_netzz
  module procedure dq_netdz
  module procedure dq_netzd
  module procedure dq_netxz
  module procedure dq_netzx
end interface

interface operator (.le.)
  module procedure dq_letqq
  module procedure dq_letdq
  module procedure dq_letqd
end interface

interface operator (.ge.)
  module procedure dq_getqq
  module procedure dq_getdq
  module procedure dq_getqd
end interface

interface operator (.lt.)
  module procedure dq_lttqq
  module procedure dq_lttdq
  module procedure dq_lttqd
end interface

interface operator (.gt.)
  module procedure dq_gttqq
  module procedure dq_gttdq
  module procedure dq_gttqd
end interface

interface abs
  module procedure dq_absq
  module procedure dq_absz
end interface

interface acos
  module procedure dq_acos
end interface

interface aimag
  module procedure dq_imag
end interface

interface aint
  module procedure dq_aint
end interface

interface anint
  module procedure dq_anint
end interface

interface asin
  module procedure dq_asin
end interface

interface atan
  module procedure dq_atan
end interface

interface atan2
  module procedure dq_atan2
end interface

interface conjg
  module procedure dq_conjg
end interface

interface cos
  module procedure dq_cos
  module procedure dq_cosz
end interface

interface cosh
  module procedure dq_cosh
end interface

interface exp
  module procedure dq_exp
  module procedure dq_expz
end interface

interface log
  module procedure dq_log
  module procedure dq_logz
end interface

interface log10
  module procedure dq_log10
end interface

interface max
  module procedure dq_maxq
  module procedure dq_maxq3
end interface

interface min
  module procedure dq_minq
  module procedure dq_minq3
end interface

interface mod
  module procedure dq_modq
end interface

interface dqcmplx
  module procedure dq_ztoz
  module procedure dq_qtoz
  module procedure dq_dtoz
  module procedure dq_xtoz
  module procedure dq_qqtoz
  module procedure dq_ddtoz
end interface

interface dqcssh
  module procedure dq_cssh
end interface

interface dqcssn
  module procedure dq_cssn
end interface

interface dqnrt
  module procedure dq_nrt
end interface

interface dqlog2
  module procedure dq_log2
end interface

interface dqpi
  module procedure dq_pi
end interface

interface dqpolyr
  module procedure dq_poly
end interface

interface dqread
  module procedure dq_inpq
  module procedure dq_inpz
end interface

interface dqreal
  module procedure dq_qtoq
  module procedure dq_ztoq
  module procedure dq_dtoq
  module procedure dq_xtoq
  module procedure dq_atoq
  module procedure dq_itoq
end interface

interface dqwrite
  module procedure dq_outq
  module procedure dq_outz
end interface

interface qreal
  module procedure dq_qtod
  module procedure dq_ztod
end interface

interface qcmplx
  module procedure dq_qtox
  module procedure dq_ztox
end interface

interface sign
  module procedure dq_signq
end interface

interface sin
  module procedure dq_sin
  module procedure dq_sinz
end interface

interface sinh
  module procedure dq_sinh
end interface

interface sqrt
  module procedure dq_sqrtq
  module procedure dq_sqrtz
end interface

interface tan
  module procedure dq_tan
end interface

interface tanh
  module procedure dq_tanh
end interface

contains

  subroutine dqxzc (a, b)

!  This converts the DC variable A to the dqC variable B.
!  This routine is not intended to be called directly by the user.

	implicit none
    complex (dqknd) a
    real (dqknd) b(4)

    b(1) = a
    b(2) = 0.q0
    b(3) = aimag (a)
    b(4) = 0.q0
    return
  end subroutine

  subroutine dqmzc (a, b)

!  This converts the dq real variable A to the dqC variable B.
!  This routine is not intended to be called directly by the user.

    implicit none
    real (dqknd) a(2), b(4)

    b(1) = a(1)
    b(2) = a(2)
    b(3) = 0.q0
    b(4) = 0.q0
    return
  end subroutine

!  Assignment routines.

  subroutine dq_eqqq (qa, qb)
    implicit none
    type (dq_real), intent (out):: qa
    type (dq_real), intent (in):: qb
    qa%dqr(1) = qb%dqr(1)
    qa%dqr(2) = qb%dqr(2)
    return
  end subroutine

  subroutine dq_eqqz (qa, zb)
    implicit none
    type (dq_real), intent (out):: qa
    type (dq_complex), intent (in):: zb
    qa%dqr(1) = zb%dqc(1)
    qa%dqr(2) = zb%dqc(2)
    return
  end subroutine

  subroutine dq_eqdq (da, qb)
    implicit none
    real (dqknd), intent (out):: da
    type (dq_real), intent (in):: qb
    da = qb%dqr(1)
    return
  end subroutine

  subroutine dq_eqqd (qa, db)
    implicit none
    type (dq_real), intent (out):: qa
    real (dqknd), intent (in):: db
    qa%dqr(1) = db
    qa%dqr(2) = 0.q0
    return
  end subroutine

  subroutine dq_eqxq (xa, qb)
    implicit none
    complex (dqknd), intent (out):: xa
    type (dq_real), intent (in):: qb
    xa = qb%dqr(1)
    return
  end subroutine

  subroutine dq_eqqx (qa, xb)
    implicit none
    type (dq_real), intent (out):: qa
    complex (dqknd), intent (in):: xb
    qa%dqr(1) = xb
    qa%dqr(2) = 0.q0
    return
  end subroutine

  subroutine dq_eqqa (qa, ab)
    implicit none
    type (dq_real), intent (out):: qa
    character(*), intent (in):: ab
    character*80 cht
    cht = ab
    call dqinpc (cht, qa%dqr)
    return
  end subroutine

  subroutine dq_eqzq (za, qb)
    implicit none
    type (dq_complex), intent (out):: za
    type (dq_real), intent (in):: qb
    call dqmzc (qb%dqr, za%dqc)
    return
  end subroutine

  subroutine dq_eqzz (za, zb)
    implicit none
    type (dq_complex), intent (out):: za
    type (dq_complex), intent (in):: zb
    call dqceq (zb%dqc, za%dqc)
    return
  end subroutine

  subroutine dq_eqdz (da, zb)
    implicit none
    real (dqknd), intent (out):: da
    type (dq_complex), intent (in):: zb
    da = zb%dqc(1)
    return
  end subroutine

  subroutine dq_eqzd (za, db)
    implicit none
    type (dq_complex), intent (out):: za
    real (dqknd), intent (in):: db
    complex (dqknd) xb
    xb = db
    call dqxzc (xb, za%dqc)
    return
  end subroutine

  subroutine dq_eqxz (xa, zb)
    implicit none
    complex (dqknd), intent (out):: xa
    type (dq_complex), intent (in):: zb
    real (dqknd) db, dc
    db = zb%dqc(1)
    dc = zb%dqc(3)
    xa = cmplx (db, dc, dqknd)
    return
  end subroutine

  subroutine dq_eqzx (za, xb)
    implicit none
    type (dq_complex), intent (out):: za
    complex (dqknd), intent (in):: xb
    call dqxzc (xb, za%dqc)
    return
  end subroutine

!  Addition routines.

  function dq_addqq (qa, qb)
    implicit none
    type (dq_real) dq_addqq
    type (dq_real), intent (in):: qa, qb
    call dqadd (qa%dqr, qb%dqr, dq_addqq%dqr)
    return
  end function

  function dq_addqz (qa, zb)
    implicit none
    type (dq_complex) dq_addqz
    type (dq_real), intent (in):: qa
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    call dqmzc (qa%dqr, qt1)
    call dqcadd (qt1, zb%dqc, dq_addqz%dqc)
    return
  end function

  function dq_adddq (da, qb)
    implicit none
    type (dq_real):: dq_adddq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqadd (qt1, qb%dqr, dq_adddq%dqr)
    return
  end function

  function dq_addqd (qa, db)
    implicit none
    type (dq_real):: dq_addqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    real (dqknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.q0
    call dqadd (qa%dqr, qt1, dq_addqd%dqr)
    return
  end function

  function dq_addxq (xa, qb)
    implicit none
    type (dq_complex):: dq_addxq
    complex (dqknd), intent (in):: xa
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4), qt2(4)
    call dqxzc (xa, qt1)
    call dqmzc (qb%dqr, qt2)
    call dqcadd (qt1, qt2, dq_addxq%dqc)
    return
  end function

  function dq_addqx (qa, xb)
    implicit none
    type (dq_complex):: dq_addqx
    type (dq_real), intent (in):: qa
    complex (dqknd), intent (in):: xb
    real (dqknd) qt1(4), qt2(4)
    call dqmzc (qa%dqr, qt1)
    call dqxzc (xb, qt2)
    call dqcadd (qt1, qt2, dq_addqx%dqc)
    return
  end function

  function dq_addzq (za, qb)
    implicit none
    type (dq_complex):: dq_addzq
    type (dq_complex), intent (in):: za
    type (dq_real), intent (in)::qb
    real (dqknd) qt1(4)
    call dqmzc (qb%dqr, qt1)
    call dqcadd (za%dqc, qt1, dq_addzq%dqc)
    return
  end function

  function dq_addzz (za, zb)
    implicit none
    type (dq_complex):: dq_addzz
    type (dq_complex), intent (in):: za, zb
    call dqcadd (za%dqc, zb%dqc, dq_addzz%dqc)
    return
  end function

  function dq_adddz (da, zb)
    implicit none
    type (dq_complex):: dq_adddz
    real (dqknd), intent (in):: da
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    complex (dqknd) xa
    xa = da
    call dqxzc (xa, qt1)
    call dqcadd (qt1, zb%dqc, dq_adddz%dqc)
    return
  end function

  function dq_addzd (za, db)
    implicit none
    type (dq_complex):: dq_addzd
    type (dq_complex), intent (in):: za
    real (dqknd), intent (in):: db
    real (dqknd) qt1(4)
    complex (dqknd) xb
    xb = db
    call dqxzc (xb, qt1)
    call dqcadd (za%dqc, qt1, dq_addzd%dqc)
    return
  end function

  function dq_addxz (xa, zb)
    implicit none
    type (dq_complex):: dq_addxz
    complex (dqknd), intent (in):: xa
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    call dqxzc (xa, qt1)
    call dqcadd (qt1, zb%dqc, dq_addxz%dqc)
    return
  end function

  function dq_addzx (za, xb)
    implicit none
    type (dq_complex):: dq_addzx
    type (dq_complex), intent (in):: za
    complex (dqknd), intent (in):: xb
    real (dqknd) qt1(4)
    call dqxzc (xb, qt1)
    call dqcadd (za%dqc, qt1, dq_addzx%dqc)
    return
  end function

!  Subtraction routines.

  function dq_subqq (qa, qb)
    implicit none
    type (dq_real):: dq_subqq
    type (dq_real), intent (in):: qa, qb
    call dqsub (qa%dqr, qb%dqr, dq_subqq%dqr)
    return
  end function

  function dq_subqz (qa, zb)
    implicit none
    type (dq_complex):: dq_subqz
    type (dq_real), intent (in):: qa
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    call dqmzc (qa%dqr, qt1)
    call dqcsub (qt1, zb%dqc, dq_subqz%dqc)
    return
  end function

  function dq_subdq (da, qb)
    implicit none
    type (dq_real):: dq_subdq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqsub (qt1, qb%dqr, dq_subdq%dqr)
    return
  end function

  function dq_subqd (qa, db)
    implicit none
    type (dq_real):: dq_subqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    real (dqknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.q0
    call dqsub (qa%dqr, qt1, dq_subqd%dqr)
    return
  end function

  function dq_subxq (xa, qb)
    implicit none
    type (dq_complex):: dq_subxq
    complex (dqknd), intent (in):: xa
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4), qt2(4)
    call dqxzc (xa, qt1)
    call dqmzc (qb%dqr, qt2)
    call dqcsub (qt1, qt2, dq_subxq%dqc)
    return
  end function

  function dq_subqx (qa, xb)
    implicit none
    type (dq_complex):: dq_subqx
    type (dq_real), intent (in):: qa
    complex (dqknd), intent (in):: xb
    real (dqknd) qt1(4), qt2(4)
    call dqmzc (qa%dqr, qt1)
    call dqxzc (xb, qt2)
    call dqcsub (qt1, qt2, dq_subqx%dqc)
    return
  end function

  function dq_subzq (za, qb)
    implicit none
    type (dq_complex):: dq_subzq
    type (dq_complex), intent (in):: za
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4)
    call dqmzc (qb%dqr, qt1)
    call dqcsub (za%dqc, qt1, dq_subzq%dqc)
    return
  end function

  function dq_subzz (za, zb)
    implicit none
    type (dq_complex):: dq_subzz
    type (dq_complex), intent (in):: za, zb
    call dqcsub (za%dqc, zb%dqc, dq_subzz%dqc)
    return
  end function

  function dq_subdz (da, zb)
    implicit none
    type (dq_complex):: dq_subdz
    real (dqknd), intent (in):: da
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    complex (dqknd) xa
    xa = da
    call dqxzc (xa, qt1)
    call dqcsub (qt1, zb%dqc, dq_subdz%dqc)
    return
  end function

  function dq_subzd (za, db)
    implicit none
    type (dq_complex):: dq_subzd
    type (dq_complex), intent (in):: za
    real (dqknd), intent (in):: db
    real (dqknd) qt1(4)
    complex (dqknd) xb
    xb = db
    call dqxzc (xb, qt1)
    call dqcsub (za%dqc, qt1, dq_subzd%dqc)
    return
  end function

  function dq_subxz (xa, zb)
    implicit none
    type (dq_complex):: dq_subxz
    complex (dqknd), intent (in):: xa
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    call dqxzc (xa, qt1)
    call dqcsub (qt1, zb%dqc, dq_subxz%dqc)
    return
  end function

  function dq_subzx (za, xb)
    implicit none
    type (dq_complex):: dq_subzx
    type (dq_complex), intent (in):: za
    complex (dqknd), intent (in):: xb
    real (dqknd) qt1(4)
    call dqxzc (xb, qt1)
    call dqcsub (za%dqc, qt1, dq_subzx%dqc)
    return
  end function
  
!  dqR negation routines.

  function dq_negq (qa)
    implicit none
    type (dq_real):: dq_negq
    type (dq_real), intent (in):: qa
    call dqeq (qa%dqr, dq_negq%dqr)
    dq_negq%dqr(1) = - qa%dqr(1)
    dq_negq%dqr(2) = - qa%dqr(2)
    return
  end function

  function dq_negz (za)
    implicit none
    type (dq_complex):: dq_negz
    type (dq_complex), intent (in):: za
    call dqceq (za%dqc, dq_negz%dqc)
    dq_negz%dqc(1) = - za%dqc(1)
    dq_negz%dqc(2) = - za%dqc(2)
    dq_negz%dqc(3) = - za%dqc(3)
    dq_negz%dqc(4) = - za%dqc(4)
    return
  end function

!  dqR multiply routines.

  function dq_mulqq (qa, qb)
    implicit none
    type (dq_real):: dq_mulqq
    type (dq_real), intent (in):: qa, qb
    call dqmul (qa%dqr, qb%dqr, dq_mulqq%dqr)
    return
  end function

  function dq_mulqz (qa, zb)
    implicit none
    type (dq_complex):: dq_mulqz
    type (dq_real), intent (in):: qa
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    call dqmzc (qa%dqr, qt1)
    call dqcmul (qt1, zb%dqc, dq_mulqz%dqc)
    return
  end function

  function dq_muldq (da, qb)
    implicit none
    type (dq_real):: dq_muldq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    call dqmuld (qb%dqr, da, dq_muldq%dqr)
    return
  end function

  function dq_mulqd (qa, db)
    implicit none
    type (dq_real):: dq_mulqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    call dqmuld (qa%dqr, db, dq_mulqd%dqr)
    return
  end function

  function dq_mulxq (xa, qb)
    implicit none
    type (dq_complex):: dq_mulxq
    complex (dqknd), intent (in):: xa
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4), qt2(4)
    call dqxzc (xa, qt1)
    call dqmzc (qb%dqr, qt2)
    call dqcmul (qt1, qt2, dq_mulxq%dqc)
    return
  end function

  function dq_mulqx (qa, xb)
    implicit none
    type (dq_complex):: dq_mulqx
    type (dq_real), intent (in):: qa
    complex (dqknd), intent (in):: xb
    real (dqknd) qt1(4), qt2(4)
    call dqmzc (qa%dqr, qt1)
    call dqxzc (xb, qt2)
    call dqcmul (qt1, qt2, dq_mulqx%dqc)
    return
  end function

  function dq_mulzq (za, qb)
    implicit none
    type (dq_complex):: dq_mulzq
    type (dq_complex), intent (in):: za
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4)
    call dqmzc (qb%dqr, qt1)
    call dqcmul (za%dqc, qt1, dq_mulzq%dqc)
    return
  end function

  function dq_mulzz (za, zb)
    implicit none
    type (dq_complex):: dq_mulzz
    type (dq_complex), intent (in):: za, zb
    call dqcmul (za%dqc, zb%dqc, dq_mulzz%dqc)
    return
  end function

  function dq_muldz (da, zb)
    implicit none
    type (dq_complex):: dq_muldz
    real (dqknd), intent (in):: da
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    complex (dqknd) xa
    xa = da
    call dqxzc (xa, qt1)
    call dqcmul (qt1, zb%dqc, dq_muldz%dqc)
    return
  end function

  function dq_mulzd (za, db)
    implicit none
    type (dq_complex):: dq_mulzd
    type (dq_complex), intent (in):: za
    real (dqknd), intent (in):: db
    real (dqknd) qt1(4)
    complex (dqknd) xb
    xb = db
    call dqxzc (xb, qt1)
    call dqcmul (za%dqc, qt1, dq_mulzd%dqc)
    return
  end function

  function dq_mulxz (xa, zb)
    implicit none
    type (dq_complex):: dq_mulxz
    complex (dqknd), intent (in):: xa
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    call dqxzc (xa, qt1)
    call dqcmul (qt1, zb%dqc, dq_mulxz%dqc)
    return
  end function

  function dq_mulzx (za, xb)
    implicit none
    type (dq_complex):: dq_mulzx
    type (dq_complex), intent (in):: za
    complex (dqknd), intent (in):: xb
    real (dqknd) qt1(4)
    call dqxzc (xb, qt1)
    call dqcmul (za%dqc, qt1, dq_mulzx%dqc)
    return
  end function

!  dqR divide routines.

  function dq_divqq (qa, qb)
    implicit none
    type (dq_real):: dq_divqq
    type (dq_real), intent (in):: qa, qb
    call dqdiv (qa%dqr, qb%dqr, dq_divqq%dqr)
    return
  end function

  function dq_divqz (qa, zb)
    implicit none
    type (dq_complex):: dq_divqz
    type (dq_real), intent (in):: qa
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    call dqmzc (qa%dqr, qt1)
    call dqcdiv (qt1, zb%dqc, dq_divqz%dqc)
    return
  end function

  function dq_divdq (da, qb)
    implicit none
    type (dq_real):: dq_divdq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqdiv (qt1, qb%dqr, dq_divdq%dqr)
    return
  end function

  function dq_divqd (qa, db)
    implicit none
    type (dq_real):: dq_divqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    call dqdivd (qa%dqr, db, dq_divqd%dqr)
    return
  end function

  function dq_divxq (xa, qb)
    implicit none
    type (dq_complex):: dq_divxq
    complex (dqknd), intent (in):: xa
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4), qt2(4)
    call dqxzc (xa, qt1)
    call dqmzc (qb%dqr, qt2)
    call dqcdiv (qt1, qt2, dq_divxq%dqc)
    return
  end function

  function dq_divqx (qa, xb)
    implicit none
    type (dq_complex):: dq_divqx
    type (dq_real), intent (in):: qa
    complex (dqknd), intent (in):: xb
    real (dqknd) qt1(4), qt2(4)
    call dqmzc (qa%dqr, qt1)
    call dqxzc (xb, qt2)
    call dqcdiv (qt1, qt2, dq_divqx%dqc)
    return
  end function

  function dq_divzq (za, qb)
    implicit none
    type (dq_complex):: dq_divzq
    type (dq_complex), intent (in):: za
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4)
    call dqmzc (qb%dqr, qt1)
    call dqcdiv (za%dqc, qt1, dq_divzq%dqc)
    return
  end function

  function dq_divzz (za, zb)
    implicit none
    type (dq_complex):: dq_divzz
    type (dq_complex), intent (in):: za, zb
    call dqcdiv (za%dqc, zb%dqc, dq_divzz%dqc)
    return
  end function

  function dq_divdz (da, zb)
    implicit none
    type (dq_complex):: dq_divdz
    real (dqknd), intent (in):: da
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    complex (dqknd) xa
    xa = da
    call dqxzc (xa, qt1)
    call dqcdiv (qt1, zb%dqc, dq_divdz%dqc)
    return
  end function

  function dq_divzd (za, db)
    implicit none
    type (dq_complex):: dq_divzd
    type (dq_complex), intent (in):: za
    real (dqknd), intent (in):: db
    real (dqknd) qt1(4)
    complex (dqknd) xb
    xb = db
    call dqxzc (xb, qt1)
    call dqcdiv (za%dqc, qt1, dq_divzd%dqc)
    return
  end function

  function dq_divxz (xa, zb)
    implicit none
    type (dq_complex):: dq_divxz
    complex (dqknd), intent (in):: xa
    type (dq_complex), intent (in):: zb
    real (dqknd) qt1(4)
    call dqxzc (xa, qt1)
    call dqcdiv (qt1, zb%dqc, dq_divxz%dqc)
    return
  end function

  function dq_divzx (za, xb)
    implicit none
    type (dq_complex):: dq_divzx
    type (dq_complex), intent (in):: za
    complex (dqknd), intent (in):: xb
    real (dqknd) qt1(4)
    call dqxzc (xb, qt1)
    call dqcdiv (za%dqc, qt1, dq_divzx%dqc)
    return
  end function

!  dqR exponentiation routines.

  function dq_expqq (qa, qb)
    implicit none
    type (dq_real):: dq_expqq
    type (dq_real), intent (in):: qa, qb
    real (dqknd) qt1(4), qt2(4)
    call dqlog (qa%dqr, qt1)
    call dqmul (qt1, qb%dqr, qt2)
    call dqexp (qt2, dq_expqq%dqr)
    return
  end function

  function dq_expqi (qa, ib)
    implicit none
    type (dq_real):: dq_expqi
    type (dq_real), intent (in):: qa
    integer, intent (in):: ib
    call dqnpwr (qa%dqr, ib, dq_expqi%dqr)
    return
  end function

  function dq_expdq (da, qb)
    implicit none
    type (dq_real):: dq_expdq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    real (dqknd) qt1(4), qt2(4), qt3(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqlog (qt1, qt2)
    call dqmul (qt2, qb%dqr, qt3)
    call dqexp (qt3, dq_expdq%dqr)
    return
    end function

  function dq_expqd (qa, db)
    implicit none
    type (dq_real):: dq_expqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    real (dqknd) qt1(4), qt2(4)
    call dqlog (qa%dqr, qt1)
    call dqmuld (qt1, db, qt2)
    call dqexp (qt2, dq_expqd%dqr)
    return
  end function

  function dq_expzi (za, ib)
    implicit none
    type (dq_complex):: dq_expzi
    type (dq_complex), intent (in):: za
    integer, intent (in):: ib
    call dqcpwr (za%dqc, ib, dq_expzi%dqc)
    return
  end function

  function dq_expzz (za, zb)
    implicit none
    type (dq_complex):: dq_expzz
    type (dq_complex), intent (in):: za, zb
    real (dqknd) r1(4), r2(4), r3(4), r4(4), r5(4), r6(4)
    call dqmul (za%dqc(1), za%dqc(1), r1(1))
    call dqmul (za%dqc(3), za%dqc(3), r2(1))
    call dqadd (r1(1), r2(1), r3(1))
    call dqlog (r3(1), r4(1))
    call dqmuld (r4(1), 0.5q0, r5(1))
    call dqmul (zb%dqc(1), r5(1), r1(1))
    call dqang (za%dqc(1), za%dqc(3), r2(1))
    call dqmul (r2(1), zb%dqc(3), r3(1))
    call dqsub (r1(1), r3(1), r4(1))
    call dqexp (r4(1), r1(1))
    call dqmul (zb%dqc(3), r5(1), r3(1))
    call dqmul (zb%dqc(1), r2(1), r4(1))
    call dqadd (r3(1), r4(1), r6(1))
    call dqcssnf (r6(1), r3(1), r4(1))
    call dqmul (r1(1), r3(1), dq_expzz%dqc(1))
    call dqmul (r1(1), r4(1), dq_expzz%dqc(3))
    return
  end function

  function dq_exprz (qa, zb)
    implicit none
    type (dq_complex):: dq_exprz
    type (dq_real), intent (in):: qa
    type (dq_complex), intent (in):: zb
    real (dqknd) r1(4), r2(4), r3(4), r4(4), r5(4)
    call dqlog (qa%dqr(1), r2(1))
    call dqmul (r2(1), zb%dqc(1), r3(1))
    call dqexp (r3(1), r1(1))
    call dqlog (qa%dqr(1), r2(1))
    call dqmul (r2(1), zb%dqc(3), r3(1))
    call dqcssnf (r3(1), r4(1), r5(1))
    call dqmul (r1(1), r4(1), dq_exprz%dqc(1))
    call dqmul (r1(1), r5(1), dq_exprz%dqc(3))
    return
  end function

  function dq_expzr (za, qb)
    implicit none
    type (dq_complex):: dq_expzr
    type (dq_complex), intent (in):: za
    type (dq_real), intent (in):: qb
    real (dqknd) r1(4), r2(4), r3(4), r4(4), r5(4)
    call dqmul (za%dqc(1), za%dqc(1), r1(1))
    call dqmul (za%dqc(3), za%dqc(3), r2(1))
    call dqadd (r1(1), r2(1), r3(1))
    call dqlog (r3(1), r4(1))
    call dqmuld (r4(1), 0.5q0, r5(1))
    call dqmul (r5(1), qb%dqr(1), r1(1))
    call dqexp (r1(1), r2(1))
    call dqang (za%dqc(1), za%dqc(3), r3(1))
    call dqmul (qb%dqr(1), r3(1), r1(1))
    call dqcssnf (r1(1), r4(1), r5(1)) 
    call dqmul (r2(1), r4(1), dq_expzr%dqc(1))
    call dqmul (r2(1), r5(1), dq_expzr%dqc(3))
    return
  end function

!  .EQ. routines.

  function dq_eqtqq (qa, qb)
    implicit none
    logical dq_eqtqq
    type (dq_real), intent (in):: qa, qb
    integer ic
    call dqcpr (qa%dqr, qb%dqr, ic)
    if (ic .eq. 0) then
      dq_eqtqq = .true.
    else
      dq_eqtqq = .false.
    endif
    return
  end function

  function dq_eqtqz (qa, zb)
    implicit none
    logical dq_eqtqz
    type (dq_real), intent (in):: qa
    type (dq_complex), intent (in):: zb
    integer ic1, ic2
    real (dqknd) qt1(4)
    call dqmzc (qa%dqr, qt1)
    call dqcpr (qt1, zb%dqc, ic1)
    call dqcpr (qt1(3:4), zb%dqc(3:4), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      dq_eqtqz = .true.
    else
      dq_eqtqz = .false.
    endif
    return
  end function

  function dq_eqtdq (da, qb)
    implicit none
    logical dq_eqtdq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqcpr (qt1, qb%dqr, ic)
    if (ic .eq. 0) then
      dq_eqtdq = .true.
    else
      dq_eqtdq = .false.
    endif
    return
  end function

  function dq_eqtqd (qa, db)
    implicit none
    logical dq_eqtqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.q0
    call dqcpr (qa%dqr, qt1, ic)
    if (ic .eq. 0) then
      dq_eqtqd = .true.
    else
      dq_eqtqd = .false.
    endif
    return
  end function

  function dq_eqtxq (xa, qb)
    implicit none
    logical dq_eqtxq
    complex (dqknd), intent (in):: xa
    type (dq_real), intent (in):: qb
    integer ic1, ic2
    real (dqknd) qt1(4), qt2(4)
    call dqxzc (xa, qt1)
    call dqmzc (qb%dqr, qt2)
    call dqcpr (qt1, qt2, ic1)
    call dqcpr (qt1(3:4), qt2(3:4), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      dq_eqtxq = .true.
    else
      dq_eqtxq = .false.
    endif
    return
  end function

  function dq_eqtqx (qa, xb)
    implicit none
    logical dq_eqtqx
    type (dq_real), intent (in):: qa
    complex (dqknd), intent (in):: xb
    integer ic1, ic2
    real (dqknd) qt1(4), qt2(4)
    call dqmzc (qa%dqr, qt1)
    call dqxzc (xb, qt2)
    call dqcpr (qt1, qt2, ic1)
    call dqcpr (qt1(3:4), qt2(3:4), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      dq_eqtqx = .true.
    else
      dq_eqtqx = .false.
    endif
    return
  end function

  function dq_eqtzq (za, qb)
    implicit none
    logical dq_eqtzq
    type (dq_complex), intent (in):: za
    type (dq_real), intent (in):: qb
    integer ic1, ic2
    real (dqknd) qt1(4)
    call dqmzc (qb%dqr, qt1)
    call dqcpr (za%dqc, qt1, ic1)
    call dqcpr (za%dqc(3:4), qt1(3:4), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      dq_eqtzq = .true.
    else
      dq_eqtzq = .false.
    endif
    return
  end function

  function dq_eqtzz (za, zb)
    implicit none
    logical dq_eqtzz
    type (dq_complex), intent (in):: za, zb
    integer ic1, ic2
    call dqcpr (za%dqc, zb%dqc, ic1)
    call dqcpr (za%dqc(3:4), zb%dqc(3:4), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      dq_eqtzz = .true.
    else
      dq_eqtzz = .false.
    endif
    return
  end function

  function dq_eqtdz (da, zb)
    implicit none
    logical dq_eqtdz
    real (dqknd), intent (in):: da
    type (dq_complex), intent (in):: zb
    integer ic1, ic2
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqcpr (qt1, zb%dqc, ic1)
    call dqcpr (qt1(3:4), zb%dqc(3:4), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      dq_eqtdz = .true.
    else
      dq_eqtdz = .false.
    endif
    return
  end function

  function dq_eqtzd (za, db)
    implicit none
    logical dq_eqtzd
    type (dq_complex), intent (in):: za
    real (dqknd), intent (in):: db
    integer ic1, ic2
    real (dqknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.q0
    call dqcpr (za%dqc, qt1, ic1)
    call dqcpr (za%dqc(3:4), qt1(3:4), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      dq_eqtzd = .true.
    else
      dq_eqtzd = .false.
    endif
    return
  end function

  function dq_eqtxz (xa, zb)
    implicit none
    logical dq_eqtxz
    complex (dqknd), intent (in):: xa
    type (dq_complex), intent (in):: zb
    integer ic1, ic2
    real (dqknd) qt1(4)
    call dqxzc (xa, qt1)
    call dqcpr (qt1, zb%dqc, ic1)
    call dqcpr (qt1(3:4), zb%dqc(3:4), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      dq_eqtxz = .true.
    else
      dq_eqtxz = .false.
    endif
    return
  end function

  function dq_eqtzx (za, xb)
    implicit none
    logical dq_eqtzx
    type (dq_complex), intent (in):: za
    complex (dqknd), intent (in):: xb
    integer ic1, ic2
    real (dqknd) qt1(4)
    call dqxzc (xb, qt1)
    call dqcpr (za%dqc, qt1, ic1)
    call dqcpr (za%dqc(3:4), qt1(3:4), ic2)
    if (ic1 .eq. 0 .and. ic2 .eq. 0) then
      dq_eqtzx = .true.
    else
      dq_eqtzx = .false.
    endif
    return
  end function

!  .NE. routines.

  function dq_netqq (qa, qb)
    implicit none
    logical dq_netqq
    type (dq_real), intent (in):: qa, qb
    integer ic
    call dqcpr (qa%dqr, qb%dqr, ic)
    if (ic .ne. 0) then
      dq_netqq = .true.
    else
      dq_netqq = .false.
    endif
    return
  end function

  function dq_netqz (qa, zb)
    implicit none
    logical dq_netqz
    type (dq_real), intent (in):: qa
    type (dq_complex), intent (in):: zb
    integer ic1, ic2
    real (dqknd) qt1(4)
    call dqmzc (qa%dqr, qt1)
    call dqcpr (qt1, zb%dqc, ic1)
    call dqcpr (qt1(3:4), zb%dqc(3:4), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      dq_netqz = .true.
    else
      dq_netqz = .false.
    endif
    return
  end function

  function dq_netdq (da, qb)
    implicit none
    logical dq_netdq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqcpr (qt1, qb%dqr, ic)
    if (ic .ne. 0) then
      dq_netdq = .true.
    else
      dq_netdq = .false.
    endif
    return
  end function

  function dq_netqd (qa, db)
    implicit none
    logical dq_netqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.q0
    call dqcpr (qa%dqr, qt1, ic)
    if (ic .ne. 0) then
      dq_netqd = .true.
    else
      dq_netqd = .false.
    endif
    return
  end function

  function dq_netxq (xa, qb)
    implicit none
    logical dq_netxq
    complex (dqknd), intent (in):: xa
    type (dq_real), intent (in):: qb
    integer ic1, ic2
    real (dqknd) qt1(4), qt2(4)
    call dqxzc (xa, qt1)
    call dqmzc (qb%dqr, qt2)
    call dqcpr (qt1, qt2, ic1)
    call dqcpr (qt1(3:4), qt2(3:4), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      dq_netxq = .true.
    else
      dq_netxq = .false.
    endif
    return
  end function

  function dq_netqx (qa, xb)
    implicit none
    logical dq_netqx
    type (dq_real), intent (in):: qa
    complex (dqknd), intent (in):: xb
    integer ic1, ic2
    real (dqknd) qt1(4), qt2(4)
    call dqmzc (qa%dqr, qt1)
    call dqxzc (xb, qt2)
    call dqcpr (qt1, qt2, ic1)
    call dqcpr (qt1(3:4), qt2(3:4), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      dq_netqx = .true.
    else
      dq_netqx = .false.
    endif
    return
  end function

  function dq_netzq (za, qb)
    implicit none
    logical dq_netzq
    type (dq_complex), intent (in):: za
    type (dq_real), intent(in):: qb
    integer ic1, ic2
    real (dqknd) qt1(4)
    call dqmzc (qb%dqr, qt1)
    call dqcpr (za%dqc, qt1, ic1)
    call dqcpr (za%dqc(3:4), qt1(3:4), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      dq_netzq = .true.
    else
      dq_netzq = .false.
    endif
    return
  end function

  function dq_netzz (za, zb)
    implicit none
    logical dq_netzz
    type (dq_complex), intent (in):: za, zb
    integer ic1, ic2
    call dqcpr (za%dqc, zb%dqc, ic1)
    call dqcpr (za%dqc(3:4), zb%dqc(3:4), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      dq_netzz = .true.
    else
      dq_netzz = .false.
    endif
    return
  end function

  function dq_netdz (da, zb)
    implicit none
    logical dq_netdz
    real (dqknd), intent (in):: da
    type (dq_complex), intent (in):: zb
    integer ic1, ic2
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqcpr (qt1, zb%dqc, ic1)
    call dqcpr (qt1(3:4), zb%dqc(3:4), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      dq_netdz = .true.
    else
      dq_netdz = .false.
    endif
    return
  end function

  function dq_netzd (za, db)
    implicit none
    logical dq_netzd
    type (dq_complex), intent (in):: za
    real (dqknd), intent (in):: db
    integer ic1, ic2
    real (dqknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.q0
    call dqcpr (za%dqc, qt1, ic1)
    call dqcpr (za%dqc(3:4), qt1(3:4), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      dq_netzd = .true.
    else
      dq_netzd = .false.
    endif
    return
  end function

  function dq_netxz (xa, zb)
    implicit none
    logical dq_netxz
    complex (dqknd), intent (in):: xa
    type (dq_complex), intent (in):: zb
    integer ic1, ic2
    real (dqknd) qt1(4)
    call dqxzc (xa, qt1)
    call dqcpr (qt1, zb%dqc, ic1)
    call dqcpr (qt1(3:4), zb%dqc(3:4), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      dq_netxz = .true.
    else
      dq_netxz = .false.
    endif
    return
  end function

  function dq_netzx (za, xb)
    implicit none
    logical dq_netzx
    type (dq_complex), intent (in):: za
    complex (dqknd), intent (in):: xb
    integer ic1, ic2
    real (dqknd) qt1(4)
    call dqxzc (xb, qt1)
    call dqcpr (za%dqc, qt1, ic1)
    call dqcpr (za%dqc(3:4), qt1(3:4), ic2)
    if (ic1 .ne. 0 .or. ic2 .ne. 0) then
      dq_netzx = .true.
    else
      dq_netzx = .false.
    endif
    return
  end function

!  .LE. routines.

  function dq_letqq (qa, qb)
    implicit none
    logical dq_letqq
    type (dq_real), intent (in):: qa, qb
    integer ic
    call dqcpr (qa%dqr, qb%dqr, ic)
    if (ic .le. 0) then
      dq_letqq = .true.
    else
      dq_letqq = .false.
    endif
    return
  end function

  function dq_letdq (da, qb)
    implicit none
    logical dq_letdq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqcpr (qt1, qb%dqr, ic)
    if (ic .le. 0) then
      dq_letdq = .true.
    else
      dq_letdq = .false.
    endif
    return
  end function

  function dq_letqd (qa, db)
    implicit none
    logical dq_letqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.q0
    call dqcpr (qa%dqr, qt1, ic)
    if (ic .le. 0) then
      dq_letqd = .true.
    else
      dq_letqd = .false.
    endif
    return
  end function

!  .GE. routines.

  function dq_getqq (qa, qb)
    implicit none
    logical dq_getqq
    type (dq_real), intent (in):: qa, qb
    integer ic
    call dqcpr (qa%dqr, qb%dqr, ic)
    if (ic .ge. 0) then
      dq_getqq = .true.
    else
      dq_getqq = .false.
    endif
    return
  end function

  function dq_getdq (da, qb)
    implicit none
    logical dq_getdq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqcpr (qt1, qb%dqr, ic)
    if (ic .ge. 0) then
      dq_getdq = .true.
    else
      dq_getdq = .false.
    endif
    return
  end function

  function dq_getqd (qa, db)
    implicit none
    logical dq_getqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.q0
    call dqcpr (qa%dqr, qt1, ic)
    if (ic .ge. 0) then
      dq_getqd = .true.
    else
      dq_getqd = .false.
    endif
    return
  end function

!  .LT. routines.

  function dq_lttqq (qa, qb)
    implicit none
    logical dq_lttqq
    type (dq_real), intent (in):: qa, qb
    integer ic
    call dqcpr (qa%dqr, qb%dqr, ic)
    if (ic .lt. 0) then
      dq_lttqq = .true.
    else
      dq_lttqq = .false.
    endif
    return
  end function

  function dq_lttdq (da, qb)
    implicit none
    logical dq_lttdq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqcpr (qt1, qb%dqr, ic)
    if (ic .lt. 0) then
      dq_lttdq = .true.
    else
      dq_lttdq = .false.
    endif
    return
  end function

  function dq_lttqd (qa, db)
    implicit none
    logical dq_lttqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.q0
    call dqcpr (qa%dqr, qt1, ic)
    if (ic .lt. 0) then
      dq_lttqd = .true.
    else
      dq_lttqd = .false.
    endif
    return
  end function

!  .GT. routines.

  function dq_gttqq (qa, qb)
    implicit none
    logical dq_gttqq
    type (dq_real), intent (in):: qa, qb
    integer ic
    call dqcpr (qa%dqr, qb%dqr, ic)
    if (ic .gt. 0) then
      dq_gttqq = .true.
    else
      dq_gttqq = .false.
    endif
    return
  end function

  function dq_gttdq (da, qb)
    implicit none
    logical dq_gttdq
    real (dqknd), intent (in):: da
    type (dq_real), intent (in):: qb
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.q0
    call dqcpr (qt1, qb%dqr, ic)
    if (ic .gt. 0) then
      dq_gttdq = .true.
    else
      dq_gttdq = .false.
    endif
    return
  end function

  function dq_gttqd (qa, db)
    implicit none
    logical dq_gttqd
    type (dq_real), intent (in):: qa
    real (dqknd), intent (in):: db
    integer ic
    real (dqknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.q0
    call dqcpr (qa%dqr, qt1, ic)
    if (ic .gt. 0) then
      dq_gttqd = .true.
    else
      dq_gttqd = .false.
    endif
    return
  end function

!  Other dq subroutines and functions.

  function dq_absq (qa)
    implicit none
    type (dq_real):: dq_absq
    type (dq_real), intent (in):: qa
    call dqeq (qa%dqr, dq_absq%dqr)
    if (qa%dqr(1) .ge. 0.q0) then
      dq_absq%dqr(1) = qa%dqr(1)
      dq_absq%dqr(2) = qa%dqr(2)
    else
      dq_absq%dqr(1) = - qa%dqr(1)
      dq_absq%dqr(2) = - qa%dqr(2)
    endif
    return
  end function

  function dq_absz (za)
    implicit none
    type (dq_real):: dq_absz
    type (dq_complex), intent (in):: za
    real (dqknd) qt1(4), qt2(4), qt3(4)
    call dqmul (za%dqc, za%dqc, qt1)
    call dqmul (za%dqc(3:4), za%dqc(3:4), qt2)
    call dqadd (qt1, qt2, qt3)
    call dqsqrt (qt3, dq_absz%dqr)
    return
  end function

  function dq_acos (qa)
    implicit none
    type (dq_real):: dq_acos
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4), qt2(4), qt3(4)
    qt1(1) = 1.q0
    qt1(2) = 0.q0
    call dqmul (qa%dqr, qa%dqr, qt2)
    call dqsub (qt1, qt2, qt3)
    call dqsqrt (qt3, qt1)
    call dqang (qa%dqr, qt1, dq_acos%dqr)
    return
  end function

  function dq_aint (qa)
    implicit none
    type (dq_real):: dq_aint
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4)
    call dqinfr (qa%dqr, dq_aint%dqr, qt1)
    return
  end function

  function dq_anint (qa)
    implicit none
    type (dq_real):: dq_anint
    type (dq_real), intent (in):: qa
    call dqnint (qa%dqr, dq_anint%dqr)
    return
  end function

  function dq_asin (qa)
    implicit none
    type (dq_real):: dq_asin
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4), qt2(4), qt3(4)
    qt1(1) = 1.q0
    qt1(2) = 0.q0
    call dqmul (qa%dqr, qa%dqr, qt2)
    call dqsub (qt1, qt2, qt3)
    call dqsqrt (qt3, qt1)
    call dqang (qt1, qa%dqr, dq_asin%dqr)
    return
  end function

  function dq_atan (qa)
    implicit none
    type (dq_real):: dq_atan
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4)
    qt1(1) = 1.q0
    qt1(2) = 0.q0
    call dqang (qt1, qa%dqr, dq_atan%dqr)
    return
  end function

  function dq_atan2 (qa, qb)
    implicit none
    type (dq_real):: dq_atan2
    type (dq_real), intent (in):: qa, qb
    call dqang (qb%dqr, qa%dqr, dq_atan2%dqr)
    return
  end function

  function dq_conjg (za)
    implicit none
    type (dq_complex):: dq_conjg
    type (dq_complex), intent (in):: za
    call dqceq (za%dqc, dq_conjg%dqc)
    dq_conjg%dqc(3) = - za%dqc(3)
    dq_conjg%dqc(4) = - za%dqc(4)
    return
  end function

  function dq_cos (qa)
    implicit none
    type (dq_real):: dq_cos
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4)
    call dqcssnf (qa%dqr, dq_cos%dqr, qt1)
    return
  end function

  function dq_cosz (za)
    implicit none
    type (dq_complex):: dq_cosz
    type (dq_complex), intent (in):: za
    real (dqknd) qt1(4), qt2(4), qt3(4), qt4(4), qt5(4), qt6(4)
    call dqeq (za%dqc(3:4), qt2)
    qt2(1) = - qt2(1)
    qt2(2) = - qt2(2)
    call dqexp (qt2, qt1)
    qt3(1) = 1.q0
    qt3(2) = 0.q0
    call dqdiv (qt3, qt1, qt2)
    call dqcssnf (za%dqc, qt3, qt4)
    call dqadd (qt1, qt2, qt5)
    call dqmuld (qt5, 0.5q0, qt6)
    call dqmul (qt6, qt3, dq_cosz%dqc)
    call dqsub (qt1, qt2, qt5)
    call dqmuld (qt5, 0.5q0, qt6)
    call dqmul (qt6, qt4, dq_cosz%dqc(3:4))
    return
  end function

  function dq_cosh (qa)
    implicit none
    type (dq_real):: dq_cosh
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4)
    call dqcsshf (qa%dqr, dq_cosh%dqr, qt1)
    return
  end function

  function dq_qtod (qa)
    implicit none
    real (dqknd):: dq_qtod
    type (dq_real), intent (in):: qa
    dq_qtod = qa%dqr(1)
    return
  end function

  function dq_ztod (za)
    implicit none
    real (dqknd):: dq_ztod
    type (dq_complex), intent (in):: za
    dq_ztod = za%dqc(1)
    return
  end function

  function dq_qtox (qa, qb)
    implicit none
    complex (dqknd):: dq_qtox
    type (dq_real), intent (in):: qa, qb
    real (dqknd) da, db
    da = qa%dqr(1)
    db = qb%dqr(1)
    dq_qtox = cmplx (da, db, dqknd)
    return
  end function

  function dq_ztox (za)
    implicit none
    complex (dqknd):: dq_ztox
    type (dq_complex), intent (in):: za
    real (dqknd) da, db
    da = za%dqc(1)
    db = za%dqc(3)
    dq_ztox = cmplx (da, db, dqknd)
    return
  end function

  function dq_exp (qa)
    implicit none
    type (dq_real):: dq_exp
    type (dq_real), intent (in):: qa
    call dqexp (qa%dqr, dq_exp%dqr)
    return
  end function

  function dq_expz (za)
    implicit none
    type (dq_complex):: dq_expz
    type (dq_complex), intent (in):: za
    real (dqknd) qt1(4), qt2(4), qt3(4)
    call dqexp (za%dqc, qt1)
    call dqcssnf (za%dqc(3:4), qt2, qt3)
    call dqmul (qt1, qt2, dq_expz%dqc)
    call dqmul (qt1, qt3, dq_expz%dqc(3:4))
    return
  end function

  function dq_imag (za)
    implicit none
    type (dq_real):: dq_imag
    type (dq_complex), intent (in):: za
    call dqeq (za%dqc(3:4), dq_imag%dqr)
    return
  end function

  function dq_log (qa)
    implicit none
    type (dq_real):: dq_log
    type (dq_real), intent (in):: qa
    call dqlog (qa%dqr, dq_log%dqr)
    return
  end function

  function dq_logz (za)
    implicit none
    type (dq_complex):: dq_logz
    type (dq_complex), intent (in):: za
    real (dqknd) qt1(4), qt2(4), qt3(4), qt4(4)
    call dqmul (za%dqc, za%dqc, qt1)
    call dqmul (za%dqc(3:4), za%dqc(3:4), qt2)
    call dqadd (qt1, qt2, qt3)
    call dqlog (qt3, qt4)
    call dqmuld (qt4, 0.5q0, dq_logz%dqc)
    call dqang (za%dqc, za%dqc(3:4), dq_logz%dqc(3:4))
    return
  end function

  function dq_log10 (qa)
    implicit none
    type (dq_real):: dq_log10
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4), qt2(4), qt3(4)
    call dqlog (qa%dqr, qt1)
    qt2(1) = 10.q0
    qt2(2) = 0.q0
    call dqlog (qt2, qt3)
    call dqdiv (qt1, qt3, dq_log10%dqr)
    return
  end function

  function dq_maxq (qa, qb)
    implicit none
    type (dq_real):: dq_maxq
    type (dq_real), intent (in):: qa, qb
    integer ic
    call dqcpr (qa%dqr, qb%dqr, ic)
    if (ic .ge. 0) then
      call dqeq (qa%dqr, dq_maxq%dqr)
    else
      call dqeq (qb%dqr, dq_maxq%dqr)
    endif
    return
  end function


  function dq_maxq3 (qa, qb, qc)
    implicit none
    type (dq_real):: dq_maxq3
    type (dq_real), intent (in):: qa, qb, qc
    integer ic
    real (dqknd) qt0(4)
    call dqcpr (qa%dqr, qb%dqr, ic)
    if (ic .ge. 0) then
      call dqeq (qa%dqr, qt0)
    else
      call dqeq (qb%dqr, qt0)
    endif
    call dqcpr (qt0, qc%dqr, ic)
    if (ic .ge. 0) then
      call dqeq (qt0, dq_maxq3%dqr)
    else
      call dqeq (qc%dqr, dq_maxq3%dqr)
    endif
    return
  end function

  function dq_minq (qa, qb)
    implicit none
    type (dq_real):: dq_minq
    type (dq_real), intent (in):: qa, qb
    integer ic
    call dqcpr (qa%dqr, qb%dqr, ic)
    if (ic .lt. 0) then
      call dqeq (qa%dqr, dq_minq%dqr)
    else
      call dqeq (qb%dqr, dq_minq%dqr)
    endif
    return
  end function

  function dq_minq3 (qa, qb, qc)
    implicit none
    type (dq_real):: dq_minq3
    type (dq_real), intent (in):: qa, qb, qc
    integer ic
    real (dqknd) qt0(4)
    call dqcpr (qa%dqr, qb%dqr, ic)
    if (ic .lt. 0) then
      call dqeq (qa%dqr, qt0)
    else
      call dqeq (qb%dqr, qt0)
    endif
    call dqcpr (qt0, qc%dqr, ic)
    if (ic .lt. 0) then
      call dqeq (qt0, dq_minq3%dqr)
    else
      call dqeq (qc%dqr, dq_minq3%dqr)
    endif
    return
  end function

  function dq_modq (qa, qb)
    implicit none
    type (dq_real):: dq_modq
    type (dq_real), intent (in):: qa, qb
    real (dqknd) qt1(4), qt2(4), qt3(4)
    call dqdiv (qa%dqr, qb%dqr, qt1)
    call dqinfr (qt1, qt2, qt3)
    call dqmul (qb%dqr, qt2, qt1)
    call dqsub (qa%dqr, qt1, dq_modq%dqr)
    return
  end function

  function dq_qtoq (qa)
    implicit none
    type (dq_real):: dq_qtoq
    type (dq_real), intent (in):: qa
    call dqeq (qa%dqr, dq_qtoq%dqr)
    return
  end function

  function dq_ztoz (za)
    implicit none
    type (dq_complex):: dq_ztoz
    type (dq_complex), intent (in):: za
    call dqceq (za%dqc, dq_ztoz%dqc)
    return
  end function

  function dq_qtoz (qa)
    implicit none
    type (dq_complex):: dq_qtoz
    type (dq_real), intent (in):: qa
    call dqmzc (qa%dqr, dq_qtoz%dqc)
    return
  end function

  function dq_dtoz (da)
    implicit none
    type (dq_complex):: dq_dtoz
    real (dqknd), intent (in):: da
    complex (dqknd) xa
    xa = da
    call dqxzc (xa, dq_dtoz%dqc)
    return
  end function

  function dq_xtoz (xa)
    implicit none
    type (dq_complex):: dq_xtoz
    complex (dqknd), intent (in):: xa
    call dqxzc (xa, dq_xtoz%dqc)
    return
  end function

  function dq_qqtoz (qa, qb)
    implicit none
    type (dq_complex):: dq_qqtoz
    type (dq_real), intent (in):: qa, qb
    call dqqqc (qa%dqr, qb%dqr, dq_qqtoz%dqc)
    return
  end function

  function dq_ddtoz (da, db)
    implicit none
    type (dq_complex):: dq_ddtoz
    real (dqknd), intent (in):: da, db
    complex (dqknd) xa
    xa = cmplx (da, db, dqknd)
    call dqxzc (xa, dq_ddtoz%dqc)
    return
  end function

  subroutine dq_cssh (qa, qb, qc)
    implicit none
    type (dq_real), intent (in):: qa
    type (dq_real), intent (out):: qb, qc
    call dqcsshf (qa%dqr, qb%dqr, qc%dqr)
    return
  end subroutine

  subroutine dq_cssn (qa, qb, qc)
    implicit none
    type (dq_real), intent (in):: qa
    type (dq_real), intent (out):: qb, qc
    call dqcssnf (qa%dqr, qb%dqr, qc%dqr)
    return
  end subroutine

  function dq_nrt (qa, ib)
    implicit none
    type (dq_real):: dq_nrt
    type (dq_real), intent (in):: qa
    integer, intent (in):: ib
    call dqnrtf (qa%dqr, ib, dq_nrt%dqr)
    return
  end function

  function dq_log2 ()
    implicit none
    type (dq_real):: dq_log2
    call dqlog2c (dq_log2%dqr)
  end function    

  function dq_pi ()
    implicit none
    type (dq_real):: dq_pi
    call dqpic (dq_pi%dqr)
  end function    

  subroutine dq_poly (ia, qa, qb, qc)
    implicit none
    integer, intent (in):: ia
    type (dq_real), intent (in) :: qa(0:ia), qb
    type (dq_real), intent (out):: qc
    call dqpoly (ia, qa(0)%dqr, qb%dqr, qc%dqr)
    return
  end subroutine

  subroutine dq_inpq (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    implicit none
    integer, intent (in):: iu
    type (dq_real), intent (in):: q1, q2, q3, q4, q5, q6, q7, q8, q9
    optional:: q2, q3, q4, q5, q6, q7, q8, q9
    call dqinp (iu, q1%dqr)
    if (present (q2)) call dqinp (iu, q2%dqr)
    if (present (q3)) call dqinp (iu, q3%dqr)
    if (present (q4)) call dqinp (iu, q4%dqr)
    if (present (q5)) call dqinp (iu, q5%dqr)
    if (present (q6)) call dqinp (iu, q6%dqr)
    if (present (q7)) call dqinp (iu, q7%dqr)
    if (present (q8)) call dqinp (iu, q8%dqr)
    if (present (q9)) call dqinp (iu, q9%dqr)
    return
  end subroutine

  subroutine dq_inpz (iu, z1, z2, z3, z4, z5, z6, z7, z8, z9)
    implicit none
    integer, intent (in):: iu
    type (dq_complex), intent (in):: z1, z2, z3, z4, z5, z6, z7, z8, z9
    optional:: z2, z3, z4, z5, z6, z7, z8, z9
    call dqinp (iu, z1%dqc)
    call dqinp (iu, z1%dqc(3:4))
    if (present (z2)) call dqinp (iu, z2%dqc)
    if (present (z2)) call dqinp (iu, z2%dqc(3:4))
    if (present (z3)) call dqinp (iu, z3%dqc)
    if (present (z3)) call dqinp (iu, z3%dqc(3:4))
    if (present (z4)) call dqinp (iu, z4%dqc)
    if (present (z4)) call dqinp (iu, z4%dqc(3:4))
    if (present (z5)) call dqinp (iu, z5%dqc)
    if (present (z5)) call dqinp (iu, z5%dqc(3:4))
    if (present (z6)) call dqinp (iu, z6%dqc)
    if (present (z6)) call dqinp (iu, z6%dqc(3:4))
    if (present (z7)) call dqinp (iu, z7%dqc)
    if (present (z7)) call dqinp (iu, z7%dqc(3:4))
    if (present (z8)) call dqinp (iu, z8%dqc)
    if (present (z8)) call dqinp (iu, z8%dqc(3:4))
    if (present (z9)) call dqinp (iu, z9%dqc)
    if (present (z9)) call dqinp (iu, z9%dqc(3:4))
    return
  end subroutine

  function dq_ztoq (za)
    implicit none
    type (dq_real):: dq_ztoq
    type (dq_complex), intent (in):: za
    call dqeq (za%dqc, dq_ztoq%dqr)
    return
  end function

  function dq_dtoq (da)
    implicit none
    type (dq_real):: dq_dtoq
    real (dqknd), intent (in):: da
    dq_dtoq%dqr(1) = da
    dq_dtoq%dqr(2) = 0.q0
    return
  end function

  function dq_xtoq (xa)
    implicit none
    type (dq_real):: dq_xtoq
    complex (dqknd), intent (in):: xa
    dq_xtoq%dqr(1) = xa
    dq_xtoq%dqr(2) = 0.q0
    return
  end function

  function dq_atoq (aa)
    implicit none
    character(*), intent (in):: aa
    type (dq_real):: dq_atoq
    character(80) t
    t = aa
    call dqinpc (t, dq_atoq%dqr)
    return
  end function

  function dq_itoq (ia)
    implicit none
    type (dq_real):: dq_itoq
    integer, intent (in):: ia
    dq_itoq%dqr(1) = ia
    dq_itoq%dqr(2) = 0.q0
    return
  end function

  subroutine dq_outq (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    implicit none
    integer, intent (in):: iu
    type (dq_real), intent (in):: q1, q2, q3, q4, q5, q6, q7, q8, q9
    optional:: q2, q3, q4, q5, q6, q7, q8, q9
    call dqout (iu, q1%dqr)
    if (present (q2)) call dqout (iu, q2%dqr)
    if (present (q3)) call dqout (iu, q3%dqr)
    if (present (q4)) call dqout (iu, q4%dqr)
    if (present (q5)) call dqout (iu, q5%dqr)
    if (present (q6)) call dqout (iu, q6%dqr)
    if (present (q7)) call dqout (iu, q7%dqr)
    if (present (q8)) call dqout (iu, q8%dqr)
    if (present (q9)) call dqout (iu, q9%dqr)
     return
  end subroutine

  subroutine dq_outz (iu, z1, z2, z3, z4, z5, z6, z7, z8, z9)
    implicit none
    integer, intent (in):: iu
    type (dq_complex), intent (in):: z1, z2, z3, z4, z5, z6, z7, z8, z9
    optional:: z2, z3, z4, z5, z6, z7, z8, z9
    call dqout (iu, z1%dqc)
    call dqout (iu, z1%dqc(3:4))
    if (present (z2)) call dqout (iu, z2%dqc)
    if (present (z2)) call dqout (iu, z2%dqc(3:4))
    if (present (z3)) call dqout (iu, z3%dqc)
    if (present (z3)) call dqout (iu, z3%dqc(3:4))
    if (present (z4)) call dqout (iu, z4%dqc)
    if (present (z4)) call dqout (iu, z4%dqc(3:4))
    if (present (z5)) call dqout (iu, z5%dqc)
    if (present (z5)) call dqout (iu, z5%dqc(3:4))
    if (present (z6)) call dqout (iu, z6%dqc)
    if (present (z6)) call dqout (iu, z6%dqc(3:4))
    if (present (z7)) call dqout (iu, z7%dqc)
    if (present (z7)) call dqout (iu, z7%dqc(3:4))
    if (present (z8)) call dqout (iu, z8%dqc)
    if (present (z8)) call dqout (iu, z8%dqc(3:4))
    if (present (z9)) call dqout (iu, z9%dqc)
    if (present (z9)) call dqout (iu, z9%dqc(3:4))
     return
  end subroutine

  function dq_signq (qa, qb)
    implicit none
    type (dq_real):: dq_signq
    type (dq_real), intent (in):: qa, qb
    call dqeq (qa%dqr, dq_signq%dqr)
    dq_signq%dqr(1) = sign (dq_signq%dqr(1), qb%dqr(1))
    return
  end function

  function dq_sin (qa)
    implicit none
    type (dq_real):: dq_sin
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4)
    call dqcssnf (qa%dqr, qt1, dq_sin%dqr)
    return
  end function

  function dq_sinz (za)
    implicit none
    type (dq_complex):: dq_sinz
    type (dq_complex), intent (in):: za
    real (dqknd) qt1(4), qt2(4), qt3(4), qt4(4), qt5(4), qt6(4)
    call dqeq (za%dqc(3:4), qt2)
    qt2(1) = - qt2(1)
    qt2(2) = - qt2(2)
    call dqexp (qt2, qt1)
    qt3(1) = 1.q0
    qt3(2) = 0.q0
    call dqdiv (qt3, qt1, qt2)
    call dqcssnf (za%dqc, qt3, qt4)
    call dqadd (qt1, qt2, qt5)
    call dqmuld (qt5, 0.5q0, qt6)
    call dqmul (qt6, qt4, dq_sinz%dqc)
    call dqsub (qt1, qt2, qt5)
    call dqmuld (qt5, -0.5q0, qt6)
    call dqmul (qt6, qt3, dq_sinz%dqc(3:4))
    return
  end function

  function dq_sinh (qa)
    implicit none
    type (dq_real):: dq_sinh
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4)
    call dqcsshf (qa%dqr, qt1, dq_sinh%dqr)
    return
  end function

  function dq_sqrtq (qa)
    implicit none
    type (dq_real):: dq_sqrtq
    type (dq_real), intent (in):: qa
    call dqsqrt (qa%dqr, dq_sqrtq%dqr)
    return
  end function

  function dq_sqrtz (za)
    implicit none
    type (dq_complex):: dq_sqrtz
    type (dq_complex), intent (in):: za
    call dqcsqrt (za%dqc, dq_sqrtz%dqc)
    return
  end function

  function dq_tan (qa)
    implicit none
    type (dq_real):: dq_tan
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4), qt2(4)
    call dqcssnf (qa%dqr, qt1, qt2)
    call dqdiv (qt2, qt1, dq_tan%dqr)
    return
  end function

  function dq_tanh (qa)
    implicit none
    type (dq_real):: dq_tanh
    type (dq_real), intent (in):: qa
    real (dqknd) qt1(4), qt2(4)
    call dqcsshf (qa%dqr, qt1, qt2)
    call dqdiv (qt2, qt1, dq_tanh%dqr)
    return
  end function

end module
