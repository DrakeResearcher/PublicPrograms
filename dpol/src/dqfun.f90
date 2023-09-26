module dqfunmod

!  DQFUN: A thread-safe double-quad precision computation package

!  This software requires IEEE 128-bit floating-point arithmetic, in hardware or
!  software (it is currently provided, for instance, by the gfortran compiler).

!  Computational routine module (DQFUNMOD).

!  Revision date:  15 Jun 2017

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

!  DESCRIPTION OF THIS MODULE (DQFUNMOD):
!    This module contains most lower-level computational routines.

!  The following notational scheme is used to designate datatypes below:

!  A   Alphabetic [i.e. ASCII]
!  D   Quad precision [i.e. REAL (KIND (0.Q0))]
!  Q   Quad-double real
!  X   Quad complex  [i.e. COMPLEX (KIND (0.Q0))]
!  Z   Quad-double complex


integer, public:: dqknd
parameter (dqknd = kind (0.q0))

contains

subroutine dqabrt
implicit none

!   This permits one to insert a call to a vendor-specific traceback routine.

stop
end subroutine

subroutine dqadd (dqa, dqb, dqc)

!   This subroutine computes dqc = dqa + dqb.

implicit none
real (dqknd) dqa(2), dqb(2), dqc(2)
real (dqknd) e, t1, t2

!   Compute dqa + dqb using Knuth's trick.

t1 = dqa(1) + dqb(1)
e = t1 - dqa(1)
t2 = ((dqb(1) - e) + (dqa(1) - (t1 - e))) + dqa(2) + dqb(2)

!   The result is t1 + t2, after normalization.

dqc(1) = t1 + t2
dqc(2) = t2 - (dqc(1) - t1)
return
end subroutine

subroutine dqang (x, y, a)

!   This computes the dq angle A subtended by the dq pair (X, Y) considered as
!   a point in the x-y plane.  This is more useful than an arctan or arcsin
!   routine, since it places the result correctly in the full circle, i.e.
!   -Pi < A <= Pi.

!   The Taylor series for Sin converges much more slowly than that of Arcsin.
!   Thus this routine does not employ Taylor series, but instead computes
!   Arccos or Arcsin by solving Cos (a) = x or Sin (a) = y using one of the
!   following Newton iterations, both of which converge to a:

!           z_{k+1} = z_k - [x - Cos (z_k)] / Sin (z_k)
!           z_{k+1} = z_k + [y - Sin (z_k)] / Cos (z_k)

!   The first is selected if Abs (x) <= Abs (y); otherwise the second is used.

implicit none
integer i, ix, iy, k, kk, nx, ny
real (dqknd) t1, t2, t3
real (dqknd) a(2), pi(2), x(2), y(2), s0(2), s1(2), s2(2), s3(2), s4(2)
save pi

data pi / &
   3.1415926535897932384626433832795027974791q+00, &
   8.6718101301237810247970440260433494722008q-35/

!   Check if both X and Y are zero.

if (x(1) .eq. 0.q0 .and. y(1) .eq. 0.q0) then
  write (6, 1)
1 format ('*** dqANG: Both arguments are zero.')
  call dqabrt
  return
endif

!   Check if one of X or Y is zero.

if (x(1) .eq. 0.q0) then
  if (y(1) .gt. 0.q0) then
    call dqmuld (pi, 0.5q0, a)
  else
    call dqmuld (pi, -0.5q0, a)
  endif
  goto 120
elseif (y(1) .eq. 0.q0) then
  if (x(1) .gt. 0.q0) then
      a(1) = 0.q0
      a(2) = 0.q0
  else
    a(1) = pi(1)
    a(2) = pi(2)
  endif
  goto 120
endif

!   Normalize x and y so that x^2 + y^2 = 1.

call dqmul (x, x, s0)
call dqmul (y, y, s1)
call dqadd (s0, s1, s2)
call dqsqrt (s2, s3)
call dqdiv (x, s3, s1)
call dqdiv (y, s3, s2)

!   Compute initial approximation of the angle.

call dqdqdpc (s1, t1)
call dqdqdpc (s2, t2)
t3 = atan2 (t2, t1)
a(1) = t3
a(2) = 0.q0

!   The smaller of x or y will be used from now on to measure convergence.
!   This selects the Newton iteration (of the two listed above) that has the
!   largest denominator.

if (abs (t1) .le. abs (t2)) then
  kk = 1
  s0(1) = s1(1)
  s0(2) = s1(2)
else
  kk = 2
  s0(1) = s2(1)
  s0(2) = s2(2)
endif

!   Perform the Newton-Raphson iteration described.

do k = 1, 3
  call dqcssnf (a, s1, s2)
  if (kk .eq. 1) then
    call dqsub (s0, s1, s3)
    call dqdiv (s3, s2, s4)
    call dqsub (a, s4, s1)
  else
    call dqsub (s0, s2, s3)
    call dqdiv (s3, s1, s4)
    call dqadd (a, s4, s1)
  endif
  a(1) = s1(1)
  a(2) = s1(2)
enddo

 120  continue

return
end subroutine

subroutine dqcadd (a, b, c)

!   This computes the sum of the dqC numbers A and B and returns the dqC
!   result in C.

implicit none
real (dqknd) a(4), b(4), c(4)

call dqadd (a, b, c)
call dqadd (a(3), b(3), c(3))

return
end subroutine

subroutine dqcdiv (a, b, c)

!   This routine divides the dq complex numbers A and B to yield the dqC
!   quotient C.

!   This routine employs the formula described in dqCMUL to save multiprecision
!   multiplications.

implicit none
real (dqknd) a(4), b(4), c(4), f(2), s0(2), s1(2), s2(2), s3(2), s4(2)

if (b(1) .eq. 0.q0 .and. b(3) .eq. 0.q0) then
  write (6, 1)
1 format ('*** dqCDIV: Divisor is zero.')
  call dqabrt
  return
endif

f(1) = 1.q0
f(2) = 0.q0
call dqmul (a, b, s0)
call dqmul (a(3), b(3), s1)
call dqadd (s0, s1, s2)
call dqsub (s0, s1, s3)
call dqadd (a, a(3), s0)
call dqsub (b, b(3), s1)
call dqmul (s0, s1, s4)
call dqsub (s4, s3, s1)
call dqmul (b, b, s0)
call dqmul (b(3), b(3), s3)
call dqadd (s0, s3, s4)
call dqdiv (f, s4, s0)
call dqmul (s2, s0, c)
call dqmul (s1, s0, c(3))

return
end subroutine

subroutine dqceq (a, b)

!   This sets the dqC number B equal to the dqC number A.

implicit none
real (dqknd) a(4), b(4)

b(1) = a(1)
b(2) = a(2)
b(3) = a(3)
b(4) = a(4)

return
end subroutine

subroutine dqcmul (a, b, c)

!   This routine multiplies the dq complex numbers A and B to yield the dqC
!   product C.

implicit none
real (dqknd) a(4), b(4), c(4), s0(2), s1(2), s2(2), s3(2)

call dqmul (a, b, s0)
call dqmul (a(3), b(3), s1)
call dqmul (a, b(3), s2)
call dqmul (a(3), b, s3)
call dqsub (s0, s1, c)
call dqadd (s2, s3, c(3))

return
end subroutine

subroutine dqcpr (a, b, ic)

!   This routine compares the dq numbers A and B and returns in IC the value
!   -1, 0, or 1 depending on whether A < B, A = B, or A > B.  It is faster
!   than merely subtracting A and B and looking at the sign of the result.

implicit none
integer ic
real (dqknd) a(2), b(2)

if (a(1) .lt. b(1)) then
  ic = -1
elseif (a(1) .eq. b(1)) then
  if (a(2) .lt. b(2)) then
    ic = -1
  elseif (a(2) .eq. b(2)) then
    ic = 0
  else
    ic = 1
  endif
else
  ic = 1
endif

return
end subroutine

subroutine dqcpwr (a, n, b)

!   This computes the N-th power of the dqC number A and returns the dqC
!   result C in B.  When N is zero, 1 is returned.  When N is negative, the
!   reciprocal of A ^ |N| is returned.

!   This routine employs the binary method for exponentiation.

implicit none
integer j, kk, kn, l1, mn, n, na1, na2, nn
real (dqknd) cl2, t1
parameter (cl2 = 1.4426950408889633q0)
real (dqknd) a(4), b(4), s0(4), s1(4), s2(4), s3(4)

if (a(1) .eq. 0.q0 .and. a(3) .eq. 0.q0) then
  if (n .ge. 0) then
    b(1) = 0.q0
    b(2) = 0.q0
    b(3) = 0.q0
    b(4) = 0.q0
    goto 120
  else
    write (6, 1)
1   format ('*** dqCPWR: Argument is zero and N is negative or zero.')
    call dqabrt
    return
  endif
endif

nn = abs (n)
if (nn .eq. 0) then
  s2(1) = 1.q0
  s2(2) = 0.q0
  s2(3) = 0.q0
  s2(4) = 0.q0
  goto 120
elseif (nn .eq. 1) then
  s2(1) = a(1)
  s2(2) = a(2)
  s2(3) = a(3)
  s2(4) = a(4)
  goto 110
elseif (nn .eq. 2) then
  call dqcmul (a, a, s2)
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN .GT. NN.

t1 = nn
mn = cl2 * log (t1) + 1.q0 + 1.q-14

s0(1) = a(1)
s0(2) = a(2)
s0(3) = a(3)
s0(4) = a(4)
s2(1) = 1.q0
s2(2) = 0.q0
s2(3) = 0.q0
s2(4) = 0.q0
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn .ne. 2 * kk) then
    call dqcmul (s2, s0, s1)
    s2(1) = s1(1)
    s2(2) = s1(2)
    s2(3) = s1(3)
    s2(4) = s2(4)
  endif
  kn = kk
  if (j .lt. mn) then
    call dqcmul (s0, s0, s1)
    s0(1) = s1(1)
    s0(2) = s1(2)
    s0(3) = s1(3)
    s0(4) = s1(4)
  endif
enddo

!   Compute reciprocal if N is negative.

110  continue

if (n .lt. 0) then
  s1(1) = 1.q0
  s1(2) = 0.q0
  s1(3) = 0.q0
  s1(4) = 0.q0
  call dqcdiv (s1, s2, s0)
  s2(1) = s0(1)
  s2(2) = s0(2)
  s2(3) = s0(3)
  s2(4) = s0(4)
endif

b(1) = s2(1)
b(2) = s2(2)
b(3) = s2(3)
b(4) = s2(4)

120  continue
return
end subroutine

subroutine dqcsqrt (a, b)

!   This routine computes the complex square root of the dqC number C.
!   This routine uses the following formula, where A1 and A2 are the real and
!   imaginary parts of A, and where R = Sqrt [A1 ^ 2 + A2 ^2]:

!      B = Sqrt [(R + A1) / 2] + I Sqrt [(R - A1) / 2]

!   If the imaginary part of A is < 0, then the imaginary part of B is also
!   set to be < 0.

implicit none
real (dqknd) a(4), b(4), s0(2), s1(2), s2(2)

if (a(1) .eq. 0.q0 .and. a(3) .eq. 0.q0) then
  b(1) = 0.q0
  b(2) = 0.q0
  b(3) = 0.q0
  b(4) = 0.q0
  goto 100
endif

call dqmul (a, a, s0)
call dqmul (a(3), a(3), s1)
call dqadd (s0, s1, s2)
call dqsqrt (s2, s0)

s1(1) = a(1)
s1(2) = a(2)
if (s1(1) .lt. 0.q0) then
  s1(1) = - s1(1)
  s1(2) = - s1(2)
endif
call dqadd (s0, s1, s2)
call dqmuld (s2, 0.5q0, s1)
call dqsqrt (s1, s0)
call dqmuld (s0, 2.q0, s1)
if (a(1) .ge. 0.q0) then
  b(1) = s0(1)
  b(2) = s0(2)
  call dqdiv (a(3), s1, b(3))
else
  call dqdiv (a(3), s1, b)
  if (b(1) .lt. 0.q0) then
    b(1) = - b(1)
    b(2) = - b(2)
  endif
  b(3) = s0(1)
  b(4) = s0(2)
  if (a(3) .lt. 0.q0) then
    b(3) = - b(3)
    b(4) = - b(4)
  endif
endif

 100  continue
return
end subroutine

subroutine dqcsshf (a, x, y)

!   This computes the hyperbolic cosine and sine of the dq number A and
!   returns the two dq results in X and Y, respectively. 

implicit none
real (dqknd) a(2), f(2), x(2), y(2), s0(2), s1(2), s2(2)

f(1) = 1.q0
f(2) = 0.q0
call dqexp (a, s0)
call dqdiv (f, s0, s1)
call dqadd (s0, s1, s2)
call dqmuld (s2, 0.5q0, x)
call dqsub (s0, s1, s2)
call dqmuld (s2, 0.5q0, y)

return
end subroutine

subroutine dqcssnf (a, x, y)

!   This computes the cosine and sine of the dq number A and returns the two dq
!   results in X and Y, respectively.

!   This routine uses the conventional Taylor's series for Sin (s):

!   Sin (s) =  s - s^3 / 3! + s^5 / 5! - s^7 / 7! ...

!   where s = t - a * pi / 2 - b * pi / 16 and the integers a and b are chosen
!   to minimize the absolute value of s.  We can then compute

!   Sin (t) = Sin (s + a * pi / 2 + b * pi / 16)
!   Cos (t) = Cos (s + a * pi / 2 + b * pi / 16)

!   by applying elementary trig identities for sums.  The sine and cosine of
!   b * pi / 16 are of the form 1/2 * Sqrt {2 +- Sqrt [2 +- Sqrt(2)]}.
!   Reducing t in this manner insures that -Pi / 32 < s <= Pi / 32, which
!   accelerates convergence in the above series.

implicit none
integer ia, ka, kb, kc, l1
real (dqknd) t1, t2
real (dqknd) a(2), eps, f(2), pi(2), x(2), y(2), s0(2), s1(2), s2(2), s3(2), s4(2), &
  s5(2), s6(2)
real (dqknd) cs(2,2,4)
save cs, pi

data cs/ &
   9.8078528040323044912618223613423904736891q-01, &
  -1.0394975140095874525415508710255803756562q-35, &
   1.9509032201612826784828486847702224839222q-01, &
  -7.4645316836445610245134327435259958597567q-36, &
   9.2387953251128675612818318939678830976530q-01, &
  -2.2942887119059255531945272075494719330082q-35, &
   3.8268343236508977172845998403039888022151q-01, &
  -1.3460168589343374112996748555663487523819q-35, &
   8.3146961230254523707878837761790570978725q-01, &
   4.6951315383980835043512627079498950545099q-35, &
   5.5557023301960222474283081394853291620053q-01, &
  -4.1825596761694293165549520721665593806746q-35, &
   7.0710678118654752440084436210484899217362q-01, &
   4.7111212743109160328460583681970528405127q-35, &
   7.0710678118654752440084436210484899217362q-01, &
   4.7111212743109160328460583681970528405127q-35/

data pi / &
   3.1415926535897932384626433832795027974791q+00, &
   8.6718101301237810247970440260433494722008q-35/

if (a(1) .eq. 0.q0) then
  x(1) = 1.q0
  x(2) = 0.q0
  y(1) = 0.q0
  y(2) = 0.q0
  goto 120
endif

eps = 10.q0 ** (-70)
f(1) = 1.q0
f(2) = 0.q0

!   Reduce to between - Pi and Pi.

call dqmuld (pi, 2.q0, s0)
call dqdiv (a, s0, s1)
call dqnint (s1, s2)
call dqsub (s1, s2, s3)

!   Determine nearest multiple of Pi / 2, and within a quadrant, the nearest
!   multiple of Pi / 16.  Through most of the rest of this subroutine, KA and
!   KB are the integers a and b of the algorithm above.

t1 = s3(1)
t2 = 4.q0 * t1
ka = nint (t2)
kb = nint (8.q0 * (t2 - ka))
t1 = (8 * ka + kb) / 32.q0
s1(1) = t1
s1(2) = 0.q0
call dqsub (s3, s1, s2)
call dqmul (s0, s2, s1)

!   Compute cosine and sine of the reduced argument s using Taylor's series.

if (s1(1) .eq. 0.q0) then
  s0(1) = 0.q0
  s0(2) = 0.q0
  goto 110
endif
s0(1) = s1(1)
s0(2) = s1(2)
call dqmul (s0, s0, s2)
l1 = 0

100  l1 = l1 + 1
if (l1 .eq. 100) then
  write (6, 1)
1 format ('*** dqCSSN: Iteration limit exceeded.')
  call dqabrt
  return
endif

t2 = - (2.q0 * l1) * (2.q0 * l1 + 1.q0)
call dqmul (s2, s1, s3)
call dqdivd (s3, t2, s1)
call dqadd (s1, s0, s3)
s0(1) = s3(1)
s0(2) = s3(2)

!   Check for convergence of the series.

if (abs (s1(1)) .gt. eps * abs (s3(1))) goto 100

!   Compute Cos (s) = Sqrt [1 - Sin^2 (s)].

110  continue
s1(1) = s0(1)
s1(2) = s0(2)
call dqmul (s0, s0, s2)
call dqsub (f, s2, s3)
call dqsqrt (s3, s0)

!   Compute cosine and sine of b * Pi / 16.

kc = abs (kb)
f(1) = 2.
if (kc .eq. 0) then
  s2(1) = 1.q0
  s2(2) = 0.q0
  s3(1) = 0.q0
  s3(2) = 0.q0
else
  s2(1) = cs(1,1,kc)
  s2(2) = cs(2,1,kc)
  s3(1) = cs(1,2,kc)
  s3(2) = cs(2,2,kc)
endif
if (kb .lt. 0) then
  s3(1) = - s3(1)
  s3(2) = - s3(2)
endif

!   Apply the trigonometric summation identities to compute cosine and sine of
!   s + b * Pi / 16.

call dqmul (s0, s2, s4)
call dqmul (s1, s3, s5)
call dqsub (s4, s5, s6)
call dqmul (s1, s2, s4)
call dqmul (s0, s3, s5)
call dqadd (s4, s5, s1)
s0(1) = s6(1)
s0(2) = s6(2)

!   This code in effect applies the trigonometric summation identities for
!   (s + b * Pi / 16) + a * Pi / 2.

if (ka .eq. 0) then
  x(1) = s0(1)
  x(2) = s0(2)
  y(1) = s1(1)
  y(2) = s1(2)
elseif (ka .eq. 1) then
  x(1) = - s1(1)
  x(2) = - s1(2)
  y(1) = s0(1)
  y(2) = s0(2)
elseif (ka .eq. -1) then
  x(1) = s1(1)
  x(2) = s1(2)
  y(1) = - s0(1)
  y(2) = - s0(2)
elseif (ka .eq. 2 .or. ka .eq. -2) then
  x(1) = - s0(1)
  x(2) = - s0(2)
  y(1) = - s1(1)
  y(2) = - s1(2)
endif

120  continue
return
end subroutine

subroutine dqcsub (a, b, c)

!   This subracts the dqC numbers A and B and returns the dqC difference in
!   C.

implicit none
real (dqknd) a(4), b(4), c(4)

call dqsub (a, b, c)
call dqsub (a(3), b(3), c(3))

return
end subroutine

subroutine dqdiv (dqa, dqb, dqc)

!   This divides the dq number dqA by the dq number dqB to yield the dq
!   quotient dqC.

implicit none
real (dqknd) dqa(2), dqb(2), dqc(2)
real (dqknd) a1, a2, b1, b2, cona, conb, c11, c2, c21, e, split, s1, s2, &
  t1, t2, t11, t12, t21, t22
parameter (split = 144115188075855873.q0)

!   Compute a dq approximation to the quotient.

s1 = dqa(1) / dqb(1)

!   This splits s1 and dqb(1) into high-order and low-order words.

cona = s1 * split
conb = dqb(1) * split
a1 = cona - (cona - s1)
b1 = conb - (conb - dqb(1))
a2 = s1 - a1
b2 = dqb(1) - b1

!   Multiply s1 * dqb(1) using Dekker's method.

c11 = s1 * dqb(1)
c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2
!>
!   Compute s1 * dqb(2) (only high-order word is needed).

c2 = s1 * dqb(2)

!   Compute (c11, c21) + c2 using Knuth's trick.

t1 = c11 + c2
e = t1 - c11
t2 = ((c2 - e) + (c11 - (t1 - e))) + c21

!   The result is t1 + t2, after normalization.

t12 = t1 + t2
t22 = t2 - (t12 - t1)

!   Compute dqa - (t12, t22) using Knuth's trick.

t11 = dqa(1) - t12
e = t11 - dqa(1)
t21 = ((-t12 - e) + (dqa(1) - (t11 - e))) + dqa(2) - t22

!   Compute high-order word of (t11, t21) and divide by dqb(1).

s2 = (t11 + t21) / dqb(1)

!   The result is s1 + s2, after normalization.

dqc(1) = s1 + s2
dqc(2) = s2 - (dqc(1) - s1)

return
end subroutine

subroutine dqdivd (dqa, db, dqc)

!   This routine divides the dq number A by the DP number B to yield
!   the dq quotient C.  

implicit none
real (dqknd) dqa(2), db, dqc(2)
real (dqknd) a1, a2, b1, b2, cona, conb, e, split, t1, t2, t11, t12, t21, t22
parameter (split = 144115188075855873.q0)

!   Compute a DP approximation to the quotient.

t1 = dqa(1) / db
!>
!   On systems with a fused multiply add, such as IBM systems, it is faster to
!   uncomment the next two lines and comment out the following lines until !>.
!   On other systems, do the opposite.

! t12 = t1 * db
! t22 = t1 * db - t12

!   This splits t1 and db into high-order and low-order words.

cona = t1 * split
conb = db * split
a1 = cona - (cona - t1)
b1 = conb - (conb - db)
a2 = t1 - a1
b2 = db - b1

!   Multiply t1 * db using Dekker's method.

t12 = t1 * db
t22 = (((a1 * b1 - t12) + a1 * b2) + a2 * b1) + a2 * b2
!>
!   Compute dqa - (t12, t22) using Knuth's trick.

t11 = dqa(1) - t12
e = t11 - dqa(1)
t21 = ((-t12 - e) + (dqa(1) - (t11 - e))) + dqa(2) - t22

!   Compute high-order word of (t11, t21) and divide by db.

t2 = (t11 + t21) / db

!   The result is t1 + t2, after normalization.

dqc(1) = t1 + t2
dqc(2) = t2 - (dqc(1) - t1)
return
end subroutine

subroutine dqdpdqc (a, b)

!   This routine converts the DP number A to dq form in B.  All bits of
!   A are recovered in B.  However, note for example that if A = 0.1q0 and N
!   is 0, then B will NOT be the dq equivalent of 1/10.

implicit none
real (dqknd) a, b(2)

b(1) = a
b(2) = 0.q0
return
end subroutine

subroutine dqeq (a, b)

!   This routine sets the dq number B equal to the dq number A. 

implicit none
real (dqknd) a(2), b(2)

b(1) = a(1)
b(2) = a(2)
end subroutine

subroutine dqeform (a, n1, n2, b)

!   This routine converts the dq number A to E format, i.e. 1P,En1.n2.
!   B is the output array (type character*1 of size n1).  N1 must exceed
!   N2 by at least 8, and N1 must not exceed 80.  N2 must not exceed 30,
!   i.e., not more than 31 significant digits.

implicit none
integer i, ln, m1, n1, n2
parameter (ln = 80)
character*1 b(n1)
character*80 cs
real (dqknd) a(2)

if (n1 .lt. 0 .or. n2 .lt. 0 .or. n1 .gt. 80 .or. n2 .gt. 30 &
  .or. n2 .gt. n1 - 8) then
  write (6, 1) n1, n2
1 format ('*** dqEFORM: Improper n1, n2 =',2i6)
  call dqabrt
endif

call dqoutc (a, cs)
m1 = n1 - n2 - 8

do i = 1, m1
  b(i) = ' '
enddo

do i = 1, n2 + 3
  b(i+m1) = cs(i+2:i+2)
enddo

do i = 1, 5
  b(i+m1+n2+3) = cs(i+35:i+35)
enddo

return
end subroutine

subroutine dqexp (a, b)

!   This computes the exponential function of the dq number A and returns the
!   dq result in B.

!   This routine uses a modification of the Taylor's series for Exp (t):

!   Exp (t) =  (1 + r + r^2 / 2! + r^3 / 3! + r^4 / 4! ...) ^ q * 2 ^ n

!   where q = 64, r = t' / q, t' = t - n Log(2) and where n is chosen so
!   that -0.5 Log(2) < t' <= 0.5 Log(2).  Reducing t mod Log(2) and
!   dividing by 64 insures that -0.004 < r <= 0.004, which accelerates
!   convergence in the above series.

implicit none
integer i, ia, l1, na, nq, nz, n1
real (dqknd) t1, t2
parameter (nq = 6)
real (dqknd) a(2), b(2), al2(2), eps, f(2), s0(2), s1(2), s2(2), s3(2), tl
save al2
data al2/ &
   6.9314718055994530941723212145817657508364q-01, &
  -7.0081394745495851634126620087716238869662q-36/

!   Check for overflows and underflows.
eps = 10.q0 ** (-70)
if (abs (a(1)) .ge. 709.q0) then
  if (a(1) .gt. 0.q0) then
    write (6, 1) a(1)
1   format ('*** dqEXP: Argument is too large',f12.6)
    call dqabrt
    return
  else
    call dqdpdqc (0.q0, b)
    goto 130
  endif
endif

f(1) = 1.q0
f(2) = 0.q0

!   Compute the reduced argument A' = A - Log(2) * Nint [A / Log(2)].  Save
!   NZ = Nint [A / Log(2)] for correcting the exponent of the final result.

call dqdiv (a, al2, s0)
call dqnint (s0, s1)
t1 = s1(1)
nz = t1 + sign (1.q-14, t1)
call dqmul (al2, s1, s2)
call dqsub (a, s2, s0)

!   Check if the reduced argument is zero.

if (s0(1) .eq. 0.q0) then
  s0(1) = 1.q0
  s0(2) = 0.q0
  l1 = 0
  goto 120
endif

!   Divide the reduced argument by 2 ^ NQ.

call dqdivd (s0, 2.q0 ** nq, s1)

!   Compute Exp using the usual Taylor series.

s2(1) = 1.q0
s2(2) = 0.q0
s3(1) = 1.q0
s3(2) = 0.q0
l1 = 0

100  l1 = l1 + 1
if (l1 .eq. 100) then
  write (6, 2)
2 format ('*** dqEXP: Iteration limit exceeded.')
  call dqabrt
  return
endif

t2 = l1
call dqmul (s2, s1, s0)
call dqdivd (s0, t2, s2)
call dqadd (s3, s2, s0)
call dqeq (s0, s3)

!   Check for convergence of the series.

if (abs (s2(1)) .gt. eps * abs (s3(1))) goto 100

!   Raise to the (2 ^ NQ)-th power.

do i = 1, nq
  call dqmul (s0, s0, s1)
  s0(1) = s1(1)
  s0(2) = s1(2)
enddo

!  Multiply by 2 ^ NZ.

120  call dqmuld (s0, 2.q0 ** nz, b)

!   Restore original precision level.

 130  continue
return
end subroutine

subroutine dqfform (a, n1, n2, b)

!   This routine converts the dq number A to F format, i.e. Fn1.n2.
!   B is the output array (type character*1 of size n1).  N1 must exceed
!   N2 by at least 3, and N1 must not exceed 80.  N2 must not exceed 30.

implicit none
integer i, ix, kx, ln, ls, lz, mx, nx, n1, n2
parameter (ln = 80)
real (dqknd) a(2)
character*1 b(n1)
character*80 c
character*80 chr80

if (n1 .lt. 0 .or. n2 .lt. 0 .or. n1 .gt. 80 .or. n2 .gt. 30 &
  .or. n1 - n2 .lt. 3) then
  write (6, 1) n1, n2
1 format ('*** dqFFORM: Improper n1, n2 =',2i6)
  call dqabrt
endif

!   Call dqoutc and extract exponent.

call dqoutc (a, c)
ix = dqdigin (c(ln-3:ln), 4)

do i = 1, n1
  b(i) = ' '
enddo

if (a(1) .ge. 0.q0) then
  ls = 0
else
  ls = 1
endif
if (ix .ge. 0 .and. a(1) .ne. 0.q0) then
  lz = 0
else
  lz = 1
endif
mx = max (ix, 0)

!   Check for overflow of field length.

if (ls + lz + mx + n2 + 2 .gt. n1) then
  do i = 1, n1
    b(i) = '*'
  enddo

  goto 200
endif

!   Check if a zero should be output.

if (a(1) .eq. 0. .or. -ix .gt. n2) then
  do i = 1, n1 - n2 - 2
    b(i) = ' '
  enddo

  b(n1-n2-1) = '0'
  b(n1-n2) = '.'

  do i = 1, n2
    b(i+n1-n2) = '0'
  enddo

  goto 200
endif

!   Process other cases.

if (a(1) .lt. 0.) b(n1-n2-mx-2) = '-'
if (ix .ge. 0) then
  b(n1-n2-ix-1) = c(4:4)
  kx = min (ln - 9, ix)

  do i = 1, kx
    b(i+n1-n2-ix-1) = c(i+5:i+5)
  enddo

  do i = kx + 1, ix
    b(i+n1-n2-ix-1) = '0'
  enddo

  b(n1-n2) = '.'
  kx = max (min (ln - 9 - ix, n2), 0)

  do i = 1, kx
    b(i+n1-n2) = c(i+ix+5:i+ix+5)
  enddo

  do i = kx + 1, n2
    b(i+n1-n2) = '0'
  enddo
else
  nx = - ix
  b(n1-n2-1) = '0'
  b(n1-n2) = '.'

  do i = 1, nx - 1
    b(i+n1-n2) = '0'
  enddo

  b(n1-n2+nx) = c(4:4)
  kx = min (ln - 8, n2 - nx)

  do i = 1, kx
    b(i+n1-n2+nx) = c(i+5:i+5)
  enddo

  do i = kx + 1, n2 - nx
    b(i+n1-n2+nx) = '0'
  enddo
endif

200   continue

return
end subroutine
      
subroutine dqinfr (a, b, c)

!   Sets B to the integer part of the dq number A and sets C equal to the
!   fractional part of A.  Note that if A = -3.3, then B = -3 and C = -0.3.

implicit none
integer ic
real (dqknd) a(2), b(2), c(2), con(2), f(2), s0(2), s1(2), t225, t112
parameter (t225 = 2.q0 ** 225, t112 = 2.q0 ** 112)
save con
data con / t225, t112/

!   Check if  A  is zero.

if (a(1) .eq. 0.q0)  then
  b(1) = 0.q0
  b(2) = 0.q0
  c(1) = 0.q0
  c(2) = 0.q0
  goto 120
endif

if (a(1) .ge. t225) then
  write (6, 1)
1 format ('*** dqINFR: Argument is too large.')
  call dqabrt
  return
endif

f(1) = 1.q0
f(2) = 0.q0
if (a(1) .gt. 0.q0) then
  call dqadd (a, con, s0)
  call dqsub (s0, con, b)
  call dqcpr (a, b, ic)
  if (ic .ge. 0) then
    call dqsub (a, b, c)
  else
    call dqsub (b, f, s1)
    b(1) = s1(1)
    b(2) = s1(2)
    call dqsub (a, b, c)
  endif
else
  call dqsub (a, con, s0)
  call dqadd (s0, con, b)
  call dqcpr (a, b, ic)
  if (ic .le. 0) then
    call dqsub (a, b, c)
  else
    call dqadd (b, f, s1)
    b(1) = s1(1)
    b(2) = s1(2)
    call dqsub (a, b, c)
  endif
endif

120  continue
return
end subroutine

subroutine dqinp (iu, a)

!   This routine reads the dq number A from logical unit IU.  The input
!   value must be placed on a single line of not more than 80 characters.

implicit none
integer iu, ln
parameter (ln = 80)
character*80 cs
real (dqknd) a(2)

read (iu, '(a)', end = 100) cs
call dqinpc (cs, a)
goto 110

100 write (6, 1)
1  format ('*** dqINP: End-of-file encountered.')
call dqabrt

110 return
end subroutine

subroutine dqinpc (a, b)

!   Converts the CHARACTER*80 array A into the dq number B.

implicit none
integer i, ib, id, ie, inz, ip, is, ix, k, ln, lnn
parameter (ln = 80)
real (dqknd) bi
character*80 a
character*1 ai
character*10 ca, dig
parameter (dig = '0123456789')
real (dqknd) b(2), f(2), s0(2), s1(2), s2(2)

id = 0
ip = -1
is = 0
inz = 0
s1(1) = 0.q0
s1(2) = 0.q0

do i = 80, 1, -1
  if (a(i:i) /= ' ') goto 90
enddo

90 continue

lnn = i

!   Scan for digits, looking for the period also.

do i = 1, lnn
  ai = a(i:i)
  if (ai .eq. ' ' .and. id == 0) then
  elseif (ai .eq. '.') then
    if (ip >= 0) goto 210
    ip = id
    inz = 1
  elseif (ai .eq. '+') then
    if (id .ne. 0 .or. ip >= 0 .or. is .ne. 0) goto 210
    is = 1
  elseif (ai .eq. '-') then
    if (id .ne. 0 .or. ip >= 0 .or. is .ne. 0) goto 210
    is = -1
  elseif (ai .eq. 'e' .or. ai .eq. 'E' .or. ai .eq. 'd' .or. ai .eq. 'D') then
    goto 100
  elseif (index (dig, ai) .eq. 0) then
    goto 210
  else
!    read (ai, '(f1.0)') bi
    bi = index (dig, ai) - 1
    if (inz > 0 .or. bi > 0.q0) then
      inz = 1
      id = id + 1
      call dqmuld (s1, 10.q0, s0)
      f(1) = bi
      f(2) = 0.q0
      call dqdpdqc (bi, f)
      call dqadd (s0, f, s1)
    endif
  endif
enddo

100   continue
if (is .eq. -1) then
  s1(1) = - s1(1)
  s1(2) = - s1(2)
endif
k = i
if (ip == -1) ip = id
ie = 0
is = 0
ca = ' '

do i = k + 1, lnn
  ai = a(i:i)
  if (ai .eq. ' ') then
  elseif (ai .eq. '+') then
    if (ie .ne. 0 .or. is .ne. 0) goto 210
    is = 1
  elseif (ai .eq. '-') then
    if (ie .ne. 0 .or. is .ne. 0) goto 210
    is = -1
  elseif (index (dig, ai) .eq. 0) then
    goto 210
  else
    ie = ie + 1
    if (ie .gt. 3) goto 210
    ca(ie:ie) = ai
  endif
enddo

! read (ca, '(i4)') ie
ie = dqdigin (ca, 4)
if (is .eq. -1) ie = - ie
ie = ie + ip - id
s0(1) = 10.q0
s0(2) = 0.q0
call dqnpwr (s0, ie, s2)
call dqmul (s1, s2, b)
goto 220

210  write (6, 1)
1 format ('*** dqINPC: Syntax error in literal string.')
call dqabrt

220  continue

return
end subroutine

subroutine dqlog (a, b)

!   This computes the natural logarithm of the dq number A and returns the dq
!   result in B.

!   The Taylor series for Log converges much more slowly than that of Exp.
!   Thus this routine does not employ Taylor series, but instead computes
!   logarithms by solving Exp (b) = a using the following Newton iteration,
!   which converges to b:

!           x_{k+1} = x_k + [a - Exp (x_k)] / Exp (x_k)

!   These iterations are performed with a maximum precision level NW that
!   is dynamically changed, approximately doubling with each iteration.

implicit none
integer k
real (dqknd) t1, t2
real (dqknd) a(2), al2(2), b(2), s0(2), s1(2), s2(2)
save al2
data al2/ &
   6.9314718055994530941723212145817657508364q-01, &
  -7.0081394745495851634126620087716238869662q-36/

if (a(1) .le. 0.q0) then
  write (6, 1)
1 format ('*** dqLOG: Argument is less than or equal to zero.')
  call dqabrt
  return
endif

!   Compute initial approximation of Log (A).

t1 = a(1)
t2 = log (t1)
b(1) = t2
b(2) = 0.q0

!   Perform the Newton-Raphson iteration described above.

do k = 1, 3
  call dqexp (b, s0)
  call dqsub (a, s0, s1)
  call dqdiv (s1, s0, s2)
  call dqadd (b, s2, s1)
  b(1) = s1(1)
  b(2) = s1(2)
enddo

120  continue

return
end subroutine

subroutine dqdqdpc (a, b)

!   This converts the dq number A to DP.

implicit none
real (dqknd) a(2), b

b = a(1)
return
end subroutine

subroutine dqlog2c (alog2d)

!   This returns Pi to quad precision.

implicit none
real (dqknd) alog2d(2), alog2c(2)
save alog2c
data alog2c/ 6.9314718055994530941723212145817657508364q-01, &
  -7.0081394745495851634126620087716238869662q-36/

alog2d(1) = alog2c(1)
alog2d(2) = alog2c(2)
return
end subroutine

subroutine dqqqc (a, b, c)

!   This converts dq numbers A and B to dqC form in C, i.e. C = A + B i.

implicit none
real (dqknd) a(2), b(2), c(4)

c(1) = a(1)
c(2) = a(2)
c(3) = b(1)
c(4) = b(2)
return
end subroutine

subroutine dqmul (dqa, dqb, dqc)

!   This routine multiplies dq numbers dqA and dqB to yield the dq product dqC.

implicit none
real (dqknd) dqa(2), dqb(2), dqc(2)
real (dqknd) a1, a2, b1, b2, cona, conb, c11, c21, c2, e, split, t1, t2
parameter (split = 144115188075855873.q0)

!   This splits dqa(1) and dqb(1) into high-order and low-order words.

cona = dqa(1) * split
conb = dqb(1) * split
a1 = cona - (cona - dqa(1))
b1 = conb - (conb - dqb(1))
a2 = dqa(1) - a1
b2 = dqb(1) - b1

!   Multilply dqa(1) * dqb(1) using Dekker's method.

c11 = dqa(1) * dqb(1)
c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2
!>
!   Compute dqa(1) * dqb(2) + dqa(2) * dqb(1) (only high-order word is needed).

c2 = dqa(1) * dqb(2) + dqa(2) * dqb(1)

!   Compute (c11, c21) + c2 using Knuth's trick, also adding low-order product.

t1 = c11 + c2
e = t1 - c11
t2 = ((c2 - e) + (c11 - (t1 - e))) + c21 + dqa(2) * dqb(2)

!   The result is t1 + t2, after normalization.

dqc(1) = t1 + t2
dqc(2) = t2 - (dqc(1) - t1)

return
end subroutine

subroutine dqmuld (dqa, db, dqc)

!   This routine multiplies the dq number dqA by the DP number DB to yield
!   the dq product dqC.

implicit none
real (dqknd) dqa(2), db, dqc(2)
real (dqknd) a1, a2, b1, b2, cona, conb, c11, c21, c2, e, split, t1, t2
parameter (split = 144115188075855873.q0)

!   This splits dqa(1) and db into high-order and low-order words.

cona = dqa(1) * split
conb = db * split
a1 = cona - (cona - dqa(1))
b1 = conb - (conb - db)
a2 = dqa(1) - a1
b2 = db - b1

!   Multilply dqa(1) * db using Dekker's method.

c11 = dqa(1) * db
c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2
!>
!   Compute dqa(2) * db (only high-order word is needed).

c2 = dqa(2) * db

!   Compute (c11, c21) + c2 using Knuth's trick.

t1 = c11 + c2
e = t1 - c11
t2 = ((c2 - e) + (c11 - (t1 - e))) + c21

!   The result is t1 + t2, after normalization.

dqc(1) = t1 + t2
dqc(2) = t2 - (dqc(1) - t1)
return
end subroutine

subroutine dqmuldd (da, db, dqc)

!   This subroutine computes dqc = da x db.

implicit none
real (dqknd) a1, a2, b1, b2, cona, conb, da, db, dqc(2), split, s1, s2
parameter (split = 144115188075855873.q0)

!>
!   On systems with a fused multiply add, such as IBM systems, it is faster to
!   uncomment the next two lines and comment out the following lines until !>.
!   On other systems, do the opposite.

! s1 = da * db
! s2 = da * db - s1

!   This splits da and db into high-order and low-order words.

cona = da * split
conb = db * split
a1 = cona - (cona - da)
b1 = conb - (conb - db)
a2 = da - a1
b2 = db - b1

!   Multiply da * db using Dekker's method.

s1 = da * db
s2 = (((a1 * b1 - s1) + a1 * b2) + a2 * b1) + a2 * b2
!>
dqc(1) = s1
dqc(2) = s2

return
end subroutine

subroutine dqnint (a, b)

!   This sets B equal to the integer nearest to the dq number A.

implicit none
real (dqknd) a(2), b(2), con(2), s0(2), t225, t112
parameter (t225 = 2.q0 ** 225, t112 = 2.q0 ** 112)
save con
data con / t225, t112/


!   Check if  A  is zero.

if (a(1) .eq. 0.q0)  then
  b(1) = 0.q0
  b(2) = 0.q0
  goto 120
endif

if (a(1) .ge. t225) then
  write (6, 1)
1 format ('*** dqINFR: Argument is too large.')
  call dqabrt
  return
endif

if (a(1) .gt. 0.q0) then
  call dqadd (a, con, s0)
  call dqsub (s0, con, b)
else
  call dqsub (a, con, s0)
  call dqadd (s0, con, b)
endif

120  continue
return
end subroutine

subroutine dqnpwr (a, n, b)

!   This computes the N-th power of the dq number A and returns the dq result
!   in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
!   of A ^ |N| is returned. 

!   This routine employs the binary method for exponentiation.

implicit none
integer j, kk, kn, l1, mn, n, na1, na2, nn
real (dqknd) cl2, t1
parameter (cl2 = 1.4426950408889633q0)
real (dqknd) a(2), b(2), s0(2), s1(2), s2(2), s3(2)

if (a(1) .eq. 0.q0) then
  if (n .ge. 0) then
    s2(1) = 0.q0
    s2(2) = 0.q0
    goto 120
  else
    write (6, 1)
1   format ('*** dqCPWR: Argument is zero and N is negative or zero.')
    call dqabrt
    return
  endif
endif

nn = abs (n)
if (nn .eq. 0) then
  s2(1) = 1.q0
  s2(2) = 0.q0
  goto 120
elseif (nn .eq. 1) then
  s2(1) = a(1)
  s2(2) = a(2)
  goto 110
elseif (nn .eq. 2) then
  call dqmul (a, a, s2)
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN .GT. NN.

t1 = nn
mn = cl2 * log (t1) + 1.q0 + 1.d-14
s0(1) = a(1)
s0(2) = a(2)
s2(1) = 1.q0
s2(2) = 0.q0
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn .ne. 2 * kk) then
    call dqmul (s2, s0, s1)
    s2(1) = s1(1)
    s2(2) = s1(2)
  endif
  kn = kk
  if (j .lt. mn) then
    call dqmul (s0, s0, s1)
    s0(1) = s1(1)
    s0(2) = s1(2)
  endif
enddo

!   Compute reciprocal if N is negative.

110  continue

if (n .lt. 0) then
  s1(1) = 1.q0
  s1(2) = 0.q0
  call dqdiv (s1, s2, s0)
  s2(1) = s0(1)
  s2(2) = s0(2)
endif

120  continue

b(1) = s2(1)
b(2) = s2(2)
  
return
end subroutine

subroutine dqnrtf (a, n, b)

!   This computes the N-th root of the dq number A and returns the dq result
!   in B.  N must be at least one.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to A ^ (-1/N):

!    X_{k+1} = X_k + (X_k / N) * (1 - A * X_k^N)

!   The reciprocal of the final approximation to A ^ (-1/N) is the N-th root.

implicit none
integer i, k, n
real (dqknd) t1, t2, tn
real (dqknd) a(2), b(2), f1(2), s0(2), s1(2)

if (a(1) .eq. 0.q0) then
  b(1) = 0.q0
  b(2) = 0.q0
  goto 140
elseif (a(1) .lt. 0.q0) then
  write (6, 1)
1 format ('*** dqNRT: Argument is negative.')
  call dqabrt
  return
endif
if (n .le. 0) then
  write (6, 2) n
2 format ('*** dqNRT: Improper value of N',i10)
  call dqabrt
  return
endif

!   Handle cases N = 1 and 2.

if (n .eq. 1) then
  b(1) = a(1)
  b(2) = a(1)
  goto 140
elseif (n .eq. 2) then
  call dqsqrt (a, b)
  goto 140
endif

f1(1) = 1.q0
f1(2) = 0.q0

!   Compute the initial approximation of A ^ (-1/N).

tn = n
t1 = a(1)
t2 = exp (- log (t1) / tn)
b(1) = t2
b(2) = 0.q0

!   Perform the Newton-Raphson iteration described above.

do k = 1, 3
  call dqnpwr (b, n, s0)
  call dqmul (a, s0, s1)
  call dqsub (f1, s1, s0)
  call dqmul (b, s0, s1)
  call dqdivd (s1, tn, s0)
  call dqadd (b, s0, s1)
  b(1) = s1(1)
  b(2) = s1(2)
enddo

!   Take the reciprocal to give final result.

call dqdiv (f1, b, s1)
b(1) = s1(1)
b(2) = s1(2)

140  continue
return
end subroutine

subroutine dqout (iu, a)

!   This routine writes the dq number A on logical unit iu using a standard
!   E format, with lines 40 characters long.

implicit none
integer i, iu, ln
parameter (ln = 80)
character*80 cs
real (dqknd) a(2)

call dqoutc (a, cs)
write (iu, '(a)') cs

return
end subroutine

subroutine dqoutc (a, b)

!   Converts the dq number A into character form in the CHARACTER*80 array B.
!   The format is analogous to the Fortran E format.

!   This routine is called by dqOUT, but it may be directly called by the user
!   if desired for custom output.

implicit none
integer i, ii, ix, ln, nx
parameter (ln = 80)
integer ib(ln)
real (dqknd) t1
character*80 b
character*10 ca, digits
parameter (digits = '0123456789')
real (dqknd) a(2), f(2), s0(2), s1(2)

f(1) = 10.q0
f(2) = 0.q0

do i = 1, ln
  ib(i) = 0
enddo

!   Determine exact power of ten for exponent.

if (a(1) .ne. 0.q0) then
  t1 = log10 (abs (a(1)))
  if (t1 .ge. 0.q0) then
    nx = t1
  else
    nx = t1 - 1.q0
  endif
  call dqnpwr (f, nx, s0)
  call dqdiv (a, s0, s1)
  if (s1(1) .lt. 0.q0) then
    s1(1) = - s1(1)
    s1(2) = - s1(2)
  endif

!   If we didn't quite get it exactly right, multiply or divide by 10 to fix.

  i = 0

100 continue

  i = i + 1
  if (s1(1) .lt. 1.q0) then
    nx = nx - 1
    call dqmuld (s1, 10.q0, s0)
    s1(1) = s0(1)
    s1(2) = s0(2)
    if (i <= 3) goto 100
  elseif (s1(1) .ge. 10.q0) then
    nx = nx + 1
    call dqdivd (s1, 10.q0, s0)
    s1(1) = s0(1)
    s1(2) = s0(2)
    goto 100
  endif
else
  nx = 0
  s1(1) = 0.q0
  s1(2) = 0.q0
endif

!   Compute digits.

do i = 1, ln - 8
  ii = s1(1)
  ib(i) = ii
  f(1) = ii
  call dqsub (s1, f, s0)
  call dqmuld (s0, 10.q0, s1)
enddo

!   Fix negative digits.

do i = ln - 8, 2, -1
  if (ib(i) .lt. 0) then
    ib(i) = ib(i) + 10
    ib(i-1) = ib(i-1) - 1
  endif
enddo

if (ib(1) .lt. 0) then
  write (6, 1) 
1 format ('dqoutc: negative leading digit')
  call dqabrt
endif

!   Round.

if (ib(ln-8) .ge. 5) then
  ib(ln-9) = ib(ln-9) + 1

  do i = ln - 9, 2, -1
    if (ib(i) .eq. 10) then
      ib(i) = 0
      ib(i-1) = ib(i-1) + 1
    endif
  enddo

  if (ib(1) .eq. 10) then
    ib(1) = 1
    nx = nx + 1
  endif
endif

!   Insert digit characters in ib.

b(1:1) = ' '
b(2:2) = ' '
if (a(1) .ge. 0.q0) then
  b(3:3) = ' '
else
  b(3:3) = '-'
endif
ii = ib(1)
b(4:4) = digits(ii+1:ii+1)
b(5:5) = '.'
b(ln:ln) = ' '

do i = 2, ln - 9
  ii = ib(i)  
  b(i+4:i+4) = digits(ii+1:ii+1)
enddo

!   Insert exponent.

190  continue
! write (ca, '(i4)') nx
ca = dqdigout (real (nx, 16), 4)
b(ln-4:ln-4) = 'E'
ii = 0

do i = 1, 4
  if (ca(i:i) /= ' ') then
    ii = ii + 1
    b(ln-4+ii:ln-4+ii) = ca(i:i)
  endif
enddo

do i = ii + 1, 4
  b(ln-4+i:ln-4+i) = ' '
enddo

return
end subroutine

subroutine dqpic (pi)

!   This returns Pi to quad precision.

implicit none
real (dqknd) pi(2), pic(2)
save pic
data pic / &
   3.1415926535897932384626433832795027974791q+00, &
   8.6718101301237810247970440260433494722008q-35/

pi(1) = pic(1)
pi(2) = pic(2)

return
end subroutine

subroutine dqpoly (n, a, x0, x)

!   This finds the root x near x0 (input) for the nth-degree polynomial whose
!   coefficients are given in the n+1-long vector a.  It may be necessary to
!   adjust eps -- default value is 1.q-65.

implicit none
integer i, it, n
real (dqknd)  a(2,0:n), ad(2,0:n), t1(2), t2(2), t3(2), t4(2), t5(2), &
  x(2), x0(2), dt1, eps
parameter (eps = 1.q-65)

do i = 0, n - 1
  dt1 = i + 1
  call dqmuld (a(1,i+1), dt1, ad(1,i))
enddo

ad(1,n) = 0.q0
ad(2,n) = 0.q0
x(1) = x0(1)
x(2) = x0(2)

do it = 1, 20
  t1(1) = 0.q0
  t1(2) = 0.q0
  t2(1) = 0.q0
  t2(2) = 0.q0
  t3(1) = 1.q0
  t3(2) = 0.q0

  do i = 0, n
    call dqmul (a(1,i), t3, t4)
    call dqadd (t1, t4, t5)
    t1(1) = t5(1)
    t1(2) = t5(2)
    call dqmul (ad(1,i), t3, t4)
    call dqadd (t2, t4, t5)
    t2(1) = t5(1)
    t2(2) = t5(2)
    call dqmul (t3, x, t4)
    t3(1) = t4(1)
    t3(2) = t4(2)
  enddo

  call dqdiv (t1, t2, t3)
  call dqsub (x, t3, t4)
  x(1) = t4(1)
  x(2) = t4(2)
  if (abs (t3(1)) .le. eps) goto 110
enddo

write (6, 1)
1 format ('dqroot: failed to converge.')
call dqabrt

110 continue

return
end subroutine

subroutine dqsqrt (a, b)

!   This computes the square root of the dq number A and returns the dq result
!   in B.

!   This subroutine employs the following formula (due to Alan Karp):

!          Sqrt(A) = (A * X) + 0.5 * [A - (A * X)^2] * X  (approx.)

!   where X is a double precision approximation to the reciprocal square root,
!   and where the multiplications A * X and [] * X are performed with only
!   double precision.

implicit none
real (dqknd) t1, t2, t3
real (dqknd) a(2), b(2), f(2), s0(2), s1(2)

if (a(1) .eq. 0.q0) then
  b(1) = 0.q0
  b(2) = 0.q0
  goto 100
endif
t1 = 1.q0 / sqrt (a(1))
t2 = a(1) * t1
call dqmuldd (t2, t2, s0)
call dqsub (a, s0, s1)
t3 = 0.5q0 * s1(1) * t1
s0(1) = t2
s0(2) = 0.q0
s1(1) = t3
s1(2) = 0.q0
call dqadd (s0, s1, b)

100 continue

return
end subroutine

subroutine dqsub (dqa, dqb, dqc)

!   This subroutine computes dqc = dqa - dqb.

implicit none
real (dqknd) dqa(2), dqb(2), dqc(2)
real (dqknd) e, t1, t2

!   Compute dqa + dqb using Knuth's trick.

t1 = dqa(1) - dqb(1)
e = t1 - dqa(1)
t2 = ((-dqb(1) - e) + (dqa(1) - (t1 - e))) + dqa(2) - dqb(2)

!   The result is t1 + t2, after normalization.

dqc(1) = t1 + t2
dqc(2) = t2 - (dqc(1) - t1)
return
end subroutine

  real (dqknd) function dqdigin (ca, n)
    implicit none
    real (dqknd) d1
    character*(*), ca
    character*16 digits
    integer i, k, n
    parameter (digits = '0123456789')

    d1 = 0.q0

    do i = 1, n
      k = index (digits, ca(i:i)) - 1
      if (k < 0) then
        write (6, *) 'dqdigin: non-digit in character string'
      elseif (k <= 9) then
        d1 = 10.q0 * d1 + k
      endif
    enddo

    dqdigin = d1
  end function

  character*16 function dqdigout (a, n)
    implicit none
    real (dqknd) a, d1, d2
    character*16 ca, digits
    parameter (digits = '0123456789')
    integer i, is, k, n

    ca = ' '
    is = sign (1.q0, a)
    d1 = abs (a)

    do i = n, 1, -1
      d2 = aint (d1 / 10.q0)
      k = 1.q0 + (d1 - 10.q0 * d2)
      d1 = d2
      ca(i:i) = digits(k:k)
      if (d1 == 0.q0) goto 100
    enddo

    i = 0

100 continue

    if (is < 0 .and. i > 1) then
      ca(i-1:i-1) = '-'
    elseif (i == 0 .or. is < 0 .and. i == 1) then
      ca = '****************'
    endif

    dqdigout = ca
    return
  end function

end module
