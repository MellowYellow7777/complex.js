// glibc implementations
// this is included for reference only and is not meant to be a standalone file

#include <math.h>
#include <math_private.h>
#include <math-underflow.h>
#include <complex.h>
#include <float.h>
#include <fenv.h>



/* Return real part of complex float type. */

FLOAT
M_DECL_FUNC (__creal) (CFLOAT z)
{
  return __real__ z;
}

declare_mgen_alias (__creal, creal)



/* Return imaginary part of complex float type. */

FLOAT
M_DECL_FUNC (__cimag) (CFLOAT z)
{
  return __imag__ z;
}

declare_mgen_alias (__cimag, cimag)



/* Return the complex absolute value of complex float type. */

FLOAT
M_DECL_FUNC (__cabs) (CFLOAT z)
{
  return M_SUF (__hypot) (__real__ z, __imag__ z);
}

declare_mgen_alias (__cabs, cabs)



/* Compute argument of complex float type. */

FLOAT
M_DECL_FUNC (__carg) (CFLOAT x)
{
  return M_SUF (__atan2) (__imag__ x, __real__ x);
}

declare_mgen_alias (__carg, carg)



/* Return complex conjugate of complex float type. */

#include <complex.h>

CFLOAT
M_DECL_FUNC (__conj) (CFLOAT z)
{
  return ~z;
}

declare_mgen_alias (__conj, conj)



/* Compute projection of complex float type value to Riemann sphere. */

CFLOAT
M_DECL_FUNC (__cproj) (CFLOAT x)
{
  if (isinf (__real__ x) || isinf (__imag__ x))
    {
      CFLOAT res;

      __real__ res = INFINITY;
      __imag__ res = M_COPYSIGN (0, __imag__ x);

      return res;
    }

  return x;
}

declare_mgen_alias (__cproj, cproj)



/* Return value of complex exponential function for a float type. */

CFLOAT
M_DECL_FUNC (__cexp) (CFLOAT x)
{
  CFLOAT retval;
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  if (__glibc_likely (rcls >= FP_ZERO))
    {
      /* Real part is finite.  */
      if (__glibc_likely (icls >= FP_ZERO))
        {
          /* Imaginary part is finite.  */
          const int t = (int) ((M_MAX_EXP - 1) * M_MLIT (M_LN2));
          FLOAT sinix, cosix;

          if (__glibc_likely (M_FABS (__imag__ x) > M_MIN))
            {
              M_SINCOS (__imag__ x, &sinix, &cosix);
            }
          else
            {
              sinix = __imag__ x;
              cosix = 1;
            }

          if (__real__ x > t)
            {
              FLOAT exp_t = M_EXP (t);
              __real__ x -= t;
              sinix *= exp_t;
              cosix *= exp_t;
              if (__real__ x > t)
                {
                  __real__ x -= t;
                  sinix *= exp_t;
                  cosix *= exp_t;
                }
            }
          if (__real__ x > t)
            {
              /* Overflow (original real part of x > 3t).  */
              __real__ retval = M_MAX * cosix;
              __imag__ retval = M_MAX * sinix;
            }
          else
            {
              FLOAT exp_val = M_EXP (__real__ x);
              __real__ retval = exp_val * cosix;
              __imag__ retval = exp_val * sinix;
            }
          math_check_force_underflow_complex (retval);
        }
      else
        {
          /* If the imaginary part is +-inf or NaN and the real part
             is not +-inf the result is NaN + iNaN.  */
          __real__ retval = M_NAN;
          __imag__ retval = M_NAN;

          feraiseexcept (FE_INVALID);
        }
    }
  else if (__glibc_likely (rcls == FP_INFINITE))
    {
      /* Real part is infinite.  */
      if (__glibc_likely (icls >= FP_ZERO))
        {
          /* Imaginary part is finite.  */
          FLOAT value = signbit (__real__ x) ? 0 : M_HUGE_VAL;

          if (icls == FP_ZERO)
            {
              /* Imaginary part is 0.0.  */
              __real__ retval = value;
              __imag__ retval = __imag__ x;
            }
          else
            {
              FLOAT sinix, cosix;

              if (__glibc_likely (M_FABS (__imag__ x) > M_MIN))
                {
                  M_SINCOS (__imag__ x, &sinix, &cosix);
                }
              else
                {
                  sinix = __imag__ x;
                  cosix = 1;
                }

              __real__ retval = M_COPYSIGN (value, cosix);
              __imag__ retval = M_COPYSIGN (value, sinix);
            }
        }
      else if (signbit (__real__ x) == 0)
        {
          __real__ retval = M_HUGE_VAL;
          __imag__ retval = __imag__ x - __imag__ x;
        }
      else
        {
          __real__ retval = 0;
          __imag__ retval = M_COPYSIGN (0, __imag__ x);
        }
    }
  else
    {
      /* If the real part is NaN the result is NaN + iNaN unless the
         imaginary part is zero.  */
      __real__ retval = M_NAN;
      if (icls == FP_ZERO)
        __imag__ retval = __imag__ x;
      else
        {
          __imag__ retval = M_NAN;

          if (rcls != FP_NAN || icls != FP_NAN)
            feraiseexcept (FE_INVALID);
        }
    }

  return retval;
}
declare_mgen_alias (__cexp, cexp)



/* Compute complex natural logarithm. */

CFLOAT
M_DECL_FUNC (__clog) (CFLOAT x)
{
  CFLOAT result;
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  if (__glibc_unlikely (rcls == FP_ZERO && icls == FP_ZERO))
    {
      /* Real and imaginary part are 0.0.  */
      __imag__ result = signbit (__real__ x) ? M_MLIT (M_PI) : 0;
      __imag__ result = M_COPYSIGN (__imag__ result, __imag__ x);
      /* Yes, the following line raises an exception.  */
      __real__ result = -1 / M_FABS (__real__ x);
    }
  else if (__glibc_likely (rcls != FP_NAN && icls != FP_NAN))
    {
      /* Neither real nor imaginary part is NaN.  */
      FLOAT absx = M_FABS (__real__ x), absy = M_FABS (__imag__ x);
      int scale = 0;

      if (absx < absy)
        {
          FLOAT t = absx;
          absx = absy;
          absy = t;
        }

      if (absx > M_MAX / 2)
        {
          scale = -1;
          absx = M_SCALBN (absx, scale);
          absy = (absy >= M_MIN * 2 ? M_SCALBN (absy, scale) : 0);
        }
      else if (absx < M_MIN && absy < M_MIN)
        {
          scale = M_MANT_DIG;
          absx = M_SCALBN (absx, scale);
          absy = M_SCALBN (absy, scale);
        }

      if (absx == 1 && scale == 0)
        {
          __real__ result = M_LOG1P (absy * absy) / 2;
          math_check_force_underflow_nonneg (__real__ result);
        }
      else if (absx > 1 && absx < 2 && absy < 1 && scale == 0)
        {
          FLOAT d2m1 = (absx - 1) * (absx + 1);
          if (absy >= M_EPSILON)
            d2m1 += absy * absy;
          __real__ result = M_LOG1P (d2m1) / 2;
        }
      else if (absx < 1
               && absx >= M_LIT (0.5)
               && absy < M_EPSILON / 2
               && scale == 0)
        {
          FLOAT d2m1 = (absx - 1) * (absx + 1);
          __real__ result = M_LOG1P (d2m1) / 2;
        }
      else if (absx < 1
               && absx >= M_LIT (0.5)
               && scale == 0
               && absx * absx + absy * absy >= M_LIT (0.5))
        {
          FLOAT d2m1 = M_SUF (__x2y2m1) (absx, absy);
          __real__ result = M_LOG1P (d2m1) / 2;
        }
      else
        {
          FLOAT d = M_HYPOT (absx, absy);
          __real__ result = M_LOG (d) - scale * M_MLIT (M_LN2);
        }

      __imag__ result = M_ATAN2 (__imag__ x, __real__ x);
    }
  else
    {
      __imag__ result = M_NAN;
      if (rcls == FP_INFINITE || icls == FP_INFINITE)
        /* Real or imaginary part is infinite.  */
        __real__ result = M_HUGE_VAL;
      else
        __real__ result = M_NAN;
    }

  return result;
}

declare_mgen_alias (__clog, clog)



/* Compute complex base 10 logarithm. */

/* log_10 (2).  */
#define LOG10_2 M_LIT (0.3010299956639811952137388947244930267682)

/* pi * log10 (e).  */
#define PI_LOG10E M_LIT (1.364376353841841347485783625431355770210)

CFLOAT
M_DECL_FUNC (__clog10) (CFLOAT x)
{
  CFLOAT result;
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  if (__glibc_unlikely (rcls == FP_ZERO && icls == FP_ZERO))
    {
      /* Real and imaginary part are 0.0.  */
      __imag__ result = signbit (__real__ x) ? PI_LOG10E : 0;
      __imag__ result = M_COPYSIGN (__imag__ result, __imag__ x);
      /* Yes, the following line raises an exception.  */
      __real__ result = -1 / M_FABS (__real__ x);
    }
  else if (__glibc_likely (rcls != FP_NAN && icls != FP_NAN))
    {
      /* Neither real nor imaginary part is NaN.  */
      FLOAT absx = M_FABS (__real__ x), absy = M_FABS (__imag__ x);
      int scale = 0;

      if (absx < absy)
        {
          FLOAT t = absx;
          absx = absy;
          absy = t;
        }

      if (absx > M_MAX / 2)
        {
          scale = -1;
          absx = M_SCALBN (absx, scale);
          absy = (absy >= M_MIN * 2 ? M_SCALBN (absy, scale) : 0);
        }
      else if (absx < M_MIN && absy < M_MIN)
        {
          scale = M_MANT_DIG;
          absx = M_SCALBN (absx, scale);
          absy = M_SCALBN (absy, scale);
        }

      if (absx == 1 && scale == 0)
        {
          __real__ result = (M_LOG1P (absy * absy)
                             * (M_MLIT (M_LOG10E) / 2));
          math_check_force_underflow_nonneg (__real__ result);
        }
      else if (absx > 1 && absx < 2 && absy < 1 && scale == 0)
        {
          FLOAT d2m1 = (absx - 1) * (absx + 1);
          if (absy >= M_EPSILON)
            d2m1 += absy * absy;
          __real__ result = M_LOG1P (d2m1) * (M_MLIT (M_LOG10E) / 2);
        }
      else if (absx < 1
               && absx >= M_LIT (0.5)
               && absy < M_EPSILON / 2
               && scale == 0)
        {
          FLOAT d2m1 = (absx - 1) * (absx + 1);
          __real__ result = M_LOG1P (d2m1) * (M_MLIT (M_LOG10E) / 2);
        }
      else if (absx < 1
               && absx >= M_LIT (0.5)
               && scale == 0
               && absx * absx + absy * absy >= M_LIT (0.5))
        {
          FLOAT d2m1 = M_SUF (__x2y2m1) (absx, absy);
          __real__ result = M_LOG1P (d2m1) * (M_MLIT (M_LOG10E) / 2);
        }
      else
        {
          FLOAT d = M_HYPOT (absx, absy);
          __real__ result = M_SUF (__ieee754_log10) (d) - scale * LOG10_2;
        }

      __imag__ result = M_MLIT (M_LOG10E) * M_ATAN2 (__imag__ x, __real__ x);
    }
  else
    {
      __imag__ result = M_NAN;
      if (rcls == FP_INFINITE || icls == FP_INFINITE)
        /* Real or imaginary part is infinite.  */
        __real__ result = M_HUGE_VAL;
      else
        __real__ result = M_NAN;
    }

  return result;
}

declare_mgen_alias (__clog10, clog10)



/* Complex power of float type. */

#include <complex.h>
#include <math.h>

CFLOAT
M_DECL_FUNC (__cpow) (CFLOAT x, CFLOAT c)
{
  return M_SUF (__cexp) (c * M_SUF (__clog) (x));
}

declare_mgen_alias (__cpow, cpow)



/* Complex square root of a float type. */

CFLOAT
M_DECL_FUNC (__csqrt) (CFLOAT x)
{
  CFLOAT res;
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  if (__glibc_unlikely (rcls <= FP_INFINITE || icls <= FP_INFINITE))
    {
      if (icls == FP_INFINITE)
        {
          __real__ res = M_HUGE_VAL;
          __imag__ res = __imag__ x;
        }
      else if (rcls == FP_INFINITE)
        {
          if (__real__ x < 0)
            {
              __real__ res = icls == FP_NAN ? M_NAN : 0;
              __imag__ res = M_COPYSIGN (M_HUGE_VAL, __imag__ x);
            }
          else
            {
              __real__ res = __real__ x;
              __imag__ res = (icls == FP_NAN
                              ? M_NAN : M_COPYSIGN (0, __imag__ x));
            }
        }
      else
        {
          __real__ res = M_NAN;
          __imag__ res = M_NAN;
        }
    }
  else
    {
      if (__glibc_unlikely (icls == FP_ZERO))
        {
          if (__real__ x < 0)
            {
              __real__ res = 0;
              __imag__ res = M_COPYSIGN (M_SQRT (-__real__ x), __imag__ x);
            }
          else
            {
              __real__ res = M_FABS (M_SQRT (__real__ x));
              __imag__ res = M_COPYSIGN (0, __imag__ x);
            }
        }
      else if (__glibc_unlikely (rcls == FP_ZERO))
        {
          FLOAT r;
          if (M_FABS (__imag__ x) >= 2 * M_MIN)
            r = M_SQRT (M_LIT (0.5) * M_FABS (__imag__ x));
          else
            r = M_LIT (0.5) * M_SQRT (2 * M_FABS (__imag__ x));

          __real__ res = r;
          __imag__ res = M_COPYSIGN (r, __imag__ x);
        }
      else
        {
          FLOAT d, r, s;
          int scale = 0;

          if (M_FABS (__real__ x) > M_MAX / 4)
            {
              scale = 1;
              __real__ x = M_SCALBN (__real__ x, -2 * scale);
              __imag__ x = M_SCALBN (__imag__ x, -2 * scale);
            }
          else if (M_FABS (__imag__ x) > M_MAX / 4)
            {
              scale = 1;
              if (M_FABS (__real__ x) >= 4 * M_MIN)
                __real__ x = M_SCALBN (__real__ x, -2 * scale);
              else
                __real__ x = 0;
              __imag__ x = M_SCALBN (__imag__ x, -2 * scale);
            }
          else if (M_FABS (__real__ x) < 2 * M_MIN
                   && M_FABS (__imag__ x) < 2 * M_MIN)
            {
              scale = -((M_MANT_DIG + 1) / 2);
              __real__ x = M_SCALBN (__real__ x, -2 * scale);
              __imag__ x = M_SCALBN (__imag__ x, -2 * scale);
            }

          d = M_HYPOT (__real__ x, __imag__ x);
          /* Use the identity   2  Re res  Im res = Im x
             to avoid cancellation error in  d +/- Re x.  */
          if (__real__ x > 0)
            {
              r = M_SQRT (M_LIT (0.5) * (d + __real__ x));
              if (scale == 1 && M_FABS (__imag__ x) < 1)
                {
                  /* Avoid possible intermediate underflow.  */
                  s = __imag__ x / r;
                  r = M_SCALBN (r, scale);
                  scale = 0;
                }
              else
                s = M_LIT (0.5) * (__imag__ x / r);
            }
          else
            {
              s = M_SQRT (M_LIT (0.5) * (d - __real__ x));
              if (scale == 1 && M_FABS (__imag__ x) < 1)
                {
                  /* Avoid possible intermediate underflow.  */
                  r = M_FABS (__imag__ x / s);
                  s = M_SCALBN (s, scale);
                  scale = 0;
                }
              else
                r = M_FABS (M_LIT (0.5) * (__imag__ x / s));
            }

          if (scale)
            {
              r = M_SCALBN (r, scale);
              s = M_SCALBN (s, scale);
            }

          math_check_force_underflow (r);
          math_check_force_underflow (s);

          __real__ res = r;
          __imag__ res = M_COPYSIGN (s, __imag__ x);
        }
    }

  return res;
}
declare_mgen_alias (__csqrt, csqrt)




/* Complex sine function for float types. */

CFLOAT
M_DECL_FUNC (__csin) (CFLOAT x)
{
  CFLOAT retval;
  int negate = signbit (__real__ x);
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  __real__ x = M_FABS (__real__ x);

  if (__glibc_likely (icls >= FP_ZERO))
    {
      /* Imaginary part is finite.  */
      if (__glibc_likely (rcls >= FP_ZERO))
        {
          /* Real part is finite.  */
          const int t = (int) ((M_MAX_EXP - 1) * M_MLIT (M_LN2));
          FLOAT sinix, cosix;

          if (__glibc_likely (__real__ x > M_MIN))
            {
              M_SINCOS (__real__ x, &sinix, &cosix);
            }
          else
            {
              sinix = __real__ x;
              cosix = 1;
            }

          if (negate)
            sinix = -sinix;

          if (M_FABS (__imag__ x) > t)
            {
              FLOAT exp_t = M_EXP (t);
              FLOAT ix = M_FABS (__imag__ x);
              if (signbit (__imag__ x))
                cosix = -cosix;
              ix -= t;
              sinix *= exp_t / 2;
              cosix *= exp_t / 2;
              if (ix > t)
                {
                  ix -= t;
                  sinix *= exp_t;
                  cosix *= exp_t;
                }
              if (ix > t)
                {
                  /* Overflow (original imaginary part of x > 3t).  */
                  __real__ retval = M_MAX * sinix;
                  __imag__ retval = M_MAX * cosix;
                }
              else
                {
                  FLOAT exp_val = M_EXP (ix);
                  __real__ retval = exp_val * sinix;
                  __imag__ retval = exp_val * cosix;
                }
            }
          else
            {
              __real__ retval = M_COSH (__imag__ x) * sinix;
              __imag__ retval = M_SINH (__imag__ x) * cosix;
            }

          math_check_force_underflow_complex (retval);
        }
      else
        {
          if (icls == FP_ZERO)
            {
              /* Imaginary part is 0.0.  */
              __real__ retval = __real__ x - __real__ x;
              __imag__ retval = __imag__ x;
            }
          else
            {
              __real__ retval = M_NAN;
              __imag__ retval = M_NAN;

              feraiseexcept (FE_INVALID);
            }
        }
    }
  else if (icls == FP_INFINITE)
    {
      /* Imaginary part is infinite.  */
      if (rcls == FP_ZERO)
        {
          /* Real part is 0.0.  */
          __real__ retval = M_COPYSIGN (0, negate ? -1 : 1);
          __imag__ retval = __imag__ x;
        }
      else if (rcls > FP_ZERO)
        {
          /* Real part is finite.  */
          FLOAT sinix, cosix;

          if (__glibc_likely (__real__ x > M_MIN))
            {
              M_SINCOS (__real__ x, &sinix, &cosix);
            }
          else
            {
              sinix = __real__ x;
              cosix = 1;
            }

          __real__ retval = M_COPYSIGN (M_HUGE_VAL, sinix);
          __imag__ retval = M_COPYSIGN (M_HUGE_VAL, cosix);

          if (negate)
            __real__ retval = -__real__ retval;
          if (signbit (__imag__ x))
            __imag__ retval = -__imag__ retval;
        }
      else
        {
          __real__ retval = __real__ x - __real__ x;
          __imag__ retval = M_HUGE_VAL;
        }
    }
  else
    {
      if (rcls == FP_ZERO)
        __real__ retval = M_COPYSIGN (0, negate ? -1 : 1);
      else
        __real__ retval = M_NAN;
      __imag__ retval = M_NAN;
    }

  return retval;
}

declare_mgen_alias (__csin, csin)



/* Return cosine of complex float type. */

#include <complex.h>
#include <fenv.h>
#include <math.h>

CFLOAT
M_DECL_FUNC (__ccos) (CFLOAT x)
{
  CFLOAT y;

  __real__ y = -__imag__ x;
  __imag__ y = __real__ x;

  return M_SUF (__ccosh) (y);
}

declare_mgen_alias (__ccos, ccos);



/* Complex tangent function for a complex float type. */

CFLOAT
M_DECL_FUNC (__ctan) (CFLOAT x)
{
  CFLOAT res;

  if (__glibc_unlikely (!isfinite (__real__ x) || !isfinite (__imag__ x)))
    {
      if (isinf (__imag__ x))
        {
          if (isfinite (__real__ x) && M_FABS (__real__ x) > 1)
            {
              FLOAT sinrx, cosrx;
              M_SINCOS (__real__ x, &sinrx, &cosrx);
              __real__ res = M_COPYSIGN (0, sinrx * cosrx);
            }
          else
            __real__ res = M_COPYSIGN (0, __real__ x);
          __imag__ res = M_COPYSIGN (1, __imag__ x);
        }
      else if (__real__ x == 0)
        {
          res = x;
        }
      else
        {
          __real__ res = M_NAN;
          if (__imag__ x == 0)
            __imag__ res = __imag__ x;
          else
            __imag__ res = M_NAN;

          if (isinf (__real__ x))
            feraiseexcept (FE_INVALID);
        }
    }
  else
    {
      FLOAT sinrx, cosrx;
      FLOAT den;
      const int t = (int) ((M_MAX_EXP - 1) * M_MLIT (M_LN2) / 2);

      /* tan(x+iy) = (sin(2x) + i*sinh(2y))/(cos(2x) + cosh(2y))
         = (sin(x)*cos(x) + i*sinh(y)*cosh(y)/(cos(x)^2 + sinh(y)^2). */

      if (__glibc_likely (M_FABS (__real__ x) > M_MIN))
        {
          M_SINCOS (__real__ x, &sinrx, &cosrx);
        }
      else
        {
          sinrx = __real__ x;
          cosrx = 1;
        }

      if (M_FABS (__imag__ x) > t)
        {
          /* Avoid intermediate overflow when the real part of the
             result may be subnormal.  Ignoring negligible terms, the
             imaginary part is +/- 1, the real part is
             sin(x)*cos(x)/sinh(y)^2 = 4*sin(x)*cos(x)/exp(2y).  */
          FLOAT exp_2t = M_EXP (2 * t);

          __imag__ res = M_COPYSIGN (1, __imag__ x);
          __real__ res = 4 * sinrx * cosrx;
          __imag__ x = M_FABS (__imag__ x);
          __imag__ x -= t;
          __real__ res /= exp_2t;
          if (__imag__ x > t)
            {
              /* Underflow (original imaginary part of x has absolute
                 value > 2t).  */
              __real__ res /= exp_2t;
            }
          else
            __real__ res /= M_EXP (2 * __imag__ x);
        }
      else
        {
          FLOAT sinhix, coshix;
          if (M_FABS (__imag__ x) > M_MIN)
            {
              sinhix = M_SINH (__imag__ x);
              coshix = M_COSH (__imag__ x);
            }
          else
            {
              sinhix = __imag__ x;
              coshix = 1;
            }

          if (M_FABS (sinhix) > M_FABS (cosrx) * M_EPSILON)
            den = cosrx * cosrx + sinhix * sinhix;
          else
            den = cosrx * cosrx;
          __real__ res = sinrx * cosrx / den;
          __imag__ res = sinhix * coshix / den;
        }
      math_check_force_underflow_complex (res);
    }

  return res;
}

declare_mgen_alias (__ctan, ctan)



/* Complex sine hyperbole function for float types. */

CFLOAT
M_DECL_FUNC (__csinh) (CFLOAT x)
{
  CFLOAT retval;
  int negate = signbit (__real__ x);
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  __real__ x = M_FABS (__real__ x);

  if (__glibc_likely (rcls >= FP_ZERO))
    {
      /* Real part is finite.  */
      if (__glibc_likely (icls >= FP_ZERO))
        {
          /* Imaginary part is finite.  */
          const int t = (int) ((M_MAX_EXP - 1) * M_MLIT (M_LN2));
          FLOAT sinix, cosix;

          if (__glibc_likely (M_FABS (__imag__ x) > M_MIN))
            {
              M_SINCOS (__imag__ x, &sinix, &cosix);
            }
          else
            {
              sinix = __imag__ x;
              cosix = 1;
            }

          if (negate)
            cosix = -cosix;

          if (M_FABS (__real__ x) > t)
            {
              FLOAT exp_t = M_EXP (t);
              FLOAT rx = M_FABS (__real__ x);
              if (signbit (__real__ x))
                cosix = -cosix;
              rx -= t;
              sinix *= exp_t / 2;
              cosix *= exp_t / 2;
              if (rx > t)
                {
                  rx -= t;
                  sinix *= exp_t;
                  cosix *= exp_t;
                }
              if (rx > t)
                {
                  /* Overflow (original real part of x > 3t).  */
                  __real__ retval = M_MAX * cosix;
                  __imag__ retval = M_MAX * sinix;
                }
              else
                {
                  FLOAT exp_val = M_EXP (rx);
                  __real__ retval = exp_val * cosix;
                  __imag__ retval = exp_val * sinix;
                }
            }
          else
            {
              __real__ retval = M_SINH (__real__ x) * cosix;
              __imag__ retval = M_COSH (__real__ x) * sinix;
            }

          math_check_force_underflow_complex (retval);
        }
      else
        {
          if (rcls == FP_ZERO)
            {
              /* Real part is 0.0.  */
              __real__ retval = M_COPYSIGN (0, negate ? -1 : 1);
              __imag__ retval = __imag__ x - __imag__ x;
            }
          else
            {
              __real__ retval = M_NAN;
              __imag__ retval = M_NAN;

              feraiseexcept (FE_INVALID);
            }
        }
    }
  else if (rcls == FP_INFINITE)
    {
      /* Real part is infinite.  */
      if (__glibc_likely (icls > FP_ZERO))
        {
          /* Imaginary part is finite.  */
          FLOAT sinix, cosix;

          if (__glibc_likely (M_FABS (__imag__ x) > M_MIN))
            {
              M_SINCOS (__imag__ x, &sinix, &cosix);
            }
          else
            {
              sinix = __imag__ x;
              cosix = 1;
            }

          __real__ retval = M_COPYSIGN (M_HUGE_VAL, cosix);
          __imag__ retval = M_COPYSIGN (M_HUGE_VAL, sinix);

          if (negate)
            __real__ retval = -__real__ retval;
        }
      else if (icls == FP_ZERO)
        {
          /* Imaginary part is 0.0.  */
          __real__ retval = negate ? -M_HUGE_VAL : M_HUGE_VAL;
          __imag__ retval = __imag__ x;
        }
      else
        {
          __real__ retval = M_HUGE_VAL;
          __imag__ retval = __imag__ x - __imag__ x;
        }
    }
  else
    {
      __real__ retval = M_NAN;
      __imag__ retval = __imag__ x == 0 ? __imag__ x : M_NAN;
    }

  return retval;
}

declare_mgen_alias (__csinh, csinh)



/* Complex cosine hyperbolic function for float types. */

CFLOAT
M_DECL_FUNC (__ccosh) (CFLOAT x)
{
  CFLOAT retval;
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  if (__glibc_likely (rcls >= FP_ZERO))
    {
      /* Real part is finite.  */
      if (__glibc_likely (icls >= FP_ZERO))
        {
          /* Imaginary part is finite.  */
          const int t = (int) ((M_MAX_EXP - 1) * M_MLIT (M_LN2));
          FLOAT sinix, cosix;

          if (__glibc_likely (M_FABS (__imag__ x) > M_MIN))
            {
              M_SINCOS (__imag__ x, &sinix, &cosix);
            }
          else
            {
              sinix = __imag__ x;
              cosix = 1;
            }

          if (M_FABS (__real__ x) > t)
            {
              FLOAT exp_t = M_EXP (t);
              FLOAT rx = M_FABS (__real__ x);
              if (signbit (__real__ x))
                sinix = -sinix;
              rx -= t;
              sinix *= exp_t / 2;
              cosix *= exp_t / 2;
              if (rx > t)
                {
                  rx -= t;
                  sinix *= exp_t;
                  cosix *= exp_t;
                }
              if (rx > t)
                {
                  /* Overflow (original real part of x > 3t).  */
                  __real__ retval = M_MAX * cosix;
                  __imag__ retval = M_MAX * sinix;
                }
              else
                {
                  FLOAT exp_val = M_EXP (rx);
                  __real__ retval = exp_val * cosix;
                  __imag__ retval = exp_val * sinix;
                }
            }
          else
            {
              __real__ retval = M_COSH (__real__ x) * cosix;
              __imag__ retval = M_SINH (__real__ x) * sinix;
            }

          math_check_force_underflow_complex (retval);
        }
      else
        {
          __imag__ retval = __real__ x == 0 ? 0 : M_NAN;
          __real__ retval = __imag__ x - __imag__ x;
        }
    }
  else if (rcls == FP_INFINITE)
    {
      /* Real part is infinite.  */
      if (__glibc_likely (icls > FP_ZERO))
        {
          /* Imaginary part is finite.  */
          FLOAT sinix, cosix;

          if (__glibc_likely (M_FABS (__imag__ x) > M_MIN))
            {
              M_SINCOS (__imag__ x, &sinix, &cosix);
            }
          else
            {
              sinix = __imag__ x;
              cosix = 1;
            }

          __real__ retval = M_COPYSIGN (M_HUGE_VAL, cosix);
          __imag__ retval = (M_COPYSIGN (M_HUGE_VAL, sinix)
                             * M_COPYSIGN (1, __real__ x));
        }
      else if (icls == FP_ZERO)
        {
          /* Imaginary part is 0.0.  */
          __real__ retval = M_HUGE_VAL;
          __imag__ retval = __imag__ x * M_COPYSIGN (1, __real__ x);
        }
      else
        {
          __real__ retval = M_HUGE_VAL;
          __imag__ retval = __imag__ x - __imag__ x;
        }
    }
  else
    {
      __real__ retval = M_NAN;
      __imag__ retval = __imag__ x == 0 ? __imag__ x : M_NAN;
    }

  return retval;
}

declare_mgen_alias (__ccosh, ccosh);



/* Complex hyperbolic tangent for float types. */

CFLOAT
M_DECL_FUNC (__ctanh) (CFLOAT x)
{
  CFLOAT res;

  if (__glibc_unlikely (!isfinite (__real__ x) || !isfinite (__imag__ x)))
    {
      if (isinf (__real__ x))
        {
          __real__ res = M_COPYSIGN (1, __real__ x);
          if (isfinite (__imag__ x) && M_FABS (__imag__ x) > 1)
            {
              FLOAT sinix, cosix;
              M_SINCOS (__imag__ x, &sinix, &cosix);
              __imag__ res = M_COPYSIGN (0, sinix * cosix);
            }
          else
            __imag__ res = M_COPYSIGN (0, __imag__ x);
        }
      else if (__imag__ x == 0)
        {
          res = x;
        }
      else
        {
          if (__real__ x == 0)
            __real__ res = __real__ x;
          else
            __real__ res = M_NAN;
          __imag__ res = M_NAN;

          if (isinf (__imag__ x))
            feraiseexcept (FE_INVALID);
        }
    }
  else
    {
      FLOAT sinix, cosix;
      FLOAT den;
      const int t = (int) ((M_MAX_EXP - 1) * M_MLIT (M_LN2) / 2);

      /* tanh(x+iy) = (sinh(2x) + i*sin(2y))/(cosh(2x) + cos(2y))
         = (sinh(x)*cosh(x) + i*sin(y)*cos(y))/(sinh(x)^2 + cos(y)^2).  */

      if (__glibc_likely (M_FABS (__imag__ x) > M_MIN))
        {
          M_SINCOS (__imag__ x, &sinix, &cosix);
        }
      else
        {
          sinix = __imag__ x;
          cosix = 1;
        }

      if (M_FABS (__real__ x) > t)
        {
          /* Avoid intermediate overflow when the imaginary part of
             the result may be subnormal.  Ignoring negligible terms,
             the real part is +/- 1, the imaginary part is
             sin(y)*cos(y)/sinh(x)^2 = 4*sin(y)*cos(y)/exp(2x).  */
          FLOAT exp_2t = M_EXP (2 * t);

          __real__ res = M_COPYSIGN (1, __real__ x);
          __imag__ res = 4 * sinix * cosix;
          __real__ x = M_FABS (__real__ x);
          __real__ x -= t;
          __imag__ res /= exp_2t;
          if (__real__ x > t)
            {
              /* Underflow (original real part of x has absolute value
                 > 2t).  */
              __imag__ res /= exp_2t;
            }
          else
            __imag__ res /= M_EXP (2 * __real__ x);
        }
      else
        {
          FLOAT sinhrx, coshrx;
          if (M_FABS (__real__ x) > M_MIN)
            {
              sinhrx = M_SINH (__real__ x);
              coshrx = M_COSH (__real__ x);
            }
          else
            {
              sinhrx = __real__ x;
              coshrx = 1;
            }

          if (M_FABS (sinhrx) > M_FABS (cosix) * M_EPSILON)
            den = sinhrx * sinhrx + cosix * cosix;
          else
            den = cosix * cosix;
          __real__ res = sinhrx * coshrx / den;
          __imag__ res = sinix * cosix / den;
        }
      math_check_force_underflow_complex (res);
    }

  return res;
}

declare_mgen_alias (__ctanh, ctanh)



/* Return arc sine of a complex float type. */

CFLOAT
M_DECL_FUNC (__casin) (CFLOAT x)
{
  CFLOAT res;

  if (isnan (__real__ x) || isnan (__imag__ x))
    {
      if (__real__ x == 0)
        {
          res = x;
        }
      else if (isinf (__real__ x) || isinf (__imag__ x))
        {
          __real__ res = M_NAN;
          __imag__ res = M_COPYSIGN (M_HUGE_VAL, __imag__ x);
        }
      else
        {
          __real__ res = M_NAN;
          __imag__ res = M_NAN;
        }
    }
  else
    {
      CFLOAT y;

      __real__ y = -__imag__ x;
      __imag__ y = __real__ x;

      y = M_SUF (__casinh) (y);

      __real__ res = __imag__ y;
      __imag__ res = -__real__ y;
    }

  return res;
}

declare_mgen_alias (__casin, casin)



/* Return cosine of a complex type. */

CFLOAT
M_DECL_FUNC (__cacos) (CFLOAT x)
{
  CFLOAT y;
  CFLOAT res;
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE
      || (rcls == FP_ZERO && icls == FP_ZERO))
    {
      y = M_SUF (__casin) (x);

      __real__ res = M_MLIT (M_PI_2) - __real__ y;
      if (__real__ res == 0)
        __real__ res = 0;
      __imag__ res = -__imag__ y;
    }
  else
    {
      __real__ y = -__imag__ x;
      __imag__ y = __real__ x;

      y = M_SUF (__kernel_casinh) (y, 1);

      __real__ res = __imag__ y;
      __imag__ res = __real__ y;
    }

  return res;
}

declare_mgen_alias (__cacos, cacos);



/* Return arc tangent of complex float type. */

CFLOAT
M_DECL_FUNC (__catan) (CFLOAT x)
{
  CFLOAT res;
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  if (__glibc_unlikely (rcls <= FP_INFINITE || icls <= FP_INFINITE))
    {
      if (rcls == FP_INFINITE)
        {
          __real__ res = M_COPYSIGN (M_MLIT (M_PI_2), __real__ x);
          __imag__ res = M_COPYSIGN (0, __imag__ x);
        }
      else if (icls == FP_INFINITE)
        {
          if (rcls >= FP_ZERO)
            __real__ res = M_COPYSIGN (M_MLIT (M_PI_2), __real__ x);
          else
            __real__ res = M_NAN;
          __imag__ res = M_COPYSIGN (0, __imag__ x);
        }
      else if (icls == FP_ZERO || icls == FP_INFINITE)
        {
          __real__ res = M_NAN;
          __imag__ res = M_COPYSIGN (0, __imag__ x);
        }
      else
        {
          __real__ res = M_NAN;
          __imag__ res = M_NAN;
        }
    }
  else if (__glibc_unlikely (rcls == FP_ZERO && icls == FP_ZERO))
    {
      res = x;
    }
  else
    {
      if (M_FABS (__real__ x) >= 16 / M_EPSILON
          || M_FABS (__imag__ x) >= 16 / M_EPSILON)
        {
          __real__ res = M_COPYSIGN (M_MLIT (M_PI_2), __real__ x);
          if (M_FABS (__real__ x) <= 1)
            __imag__ res = 1 / __imag__ x;
          else if (M_FABS (__imag__ x) <= 1)
            __imag__ res = __imag__ x / __real__ x / __real__ x;
          else
            {
              FLOAT h = M_HYPOT (__real__ x / 2, __imag__ x / 2);
              __imag__ res = __imag__ x / h / h / 4;
            }
        }
      else
        {
          FLOAT den, absx, absy;

          absx = M_FABS (__real__ x);
          absy = M_FABS (__imag__ x);
          if (absx < absy)
            {
              FLOAT t = absx;
              absx = absy;
              absy = t;
            }

          if (absy < M_EPSILON / 2)
            {
              den = (1 - absx) * (1 + absx);
              if (den == 0)
                den = 0;
            }
          else if (absx >= 1)
            den = (1 - absx) * (1 + absx) - absy * absy;
          else if (absx >= M_LIT (0.75) || absy >= M_LIT (0.5))
            den = -M_SUF (__x2y2m1) (absx, absy);
          else
            den = (1 - absx) * (1 + absx) - absy * absy;

          __real__ res = M_LIT (0.5) * M_ATAN2 (2 * __real__ x, den);

          if (M_FABS (__imag__ x) == 1
              && M_FABS (__real__ x) < M_EPSILON * M_EPSILON)
            __imag__ res = (M_COPYSIGN (M_LIT (0.5), __imag__ x)
                            * (M_MLIT (M_LN2)
                               - M_LOG (M_FABS (__real__ x))));
          else
            {
              FLOAT r2 = 0, num, f;

              if (M_FABS (__real__ x) >= M_EPSILON * M_EPSILON)
                r2 = __real__ x * __real__ x;

              num = __imag__ x + 1;
              num = r2 + num * num;

              den = __imag__ x - 1;
              den = r2 + den * den;

              f = num / den;
              if (f < M_LIT (0.5))
                __imag__ res = M_LIT (0.25) * M_LOG (f);
              else
                {
                  num = 4 * __imag__ x;
                  __imag__ res = M_LIT (0.25) * M_LOG1P (num / den);
                }
            }
        }

      math_check_force_underflow_complex (res);
    }

  return res;
}

declare_mgen_alias (__catan, catan)



/* Return arc hyperbolic sine for a complex float type. */

CFLOAT
M_DECL_FUNC (__casinh) (CFLOAT x)
{
  CFLOAT res;
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE)
    {
      if (icls == FP_INFINITE)
        {
          __real__ res = M_COPYSIGN (M_HUGE_VAL, __real__ x);

          if (rcls == FP_NAN)
            __imag__ res = M_NAN;
          else
            __imag__ res = M_COPYSIGN ((rcls >= FP_ZERO
                                        ? M_MLIT (M_PI_2) : M_MLIT (M_PI_4)),
                                       __imag__ x);
        }
      else if (rcls <= FP_INFINITE)
        {
          __real__ res = __real__ x;
          if ((rcls == FP_INFINITE && icls >= FP_ZERO)
              || (rcls == FP_NAN && icls == FP_ZERO))
            __imag__ res = M_COPYSIGN (0, __imag__ x);
          else
            __imag__ res = M_NAN;
        }
      else
        {
          __real__ res = M_NAN;
          __imag__ res = M_NAN;
        }
    }
  else if (rcls == FP_ZERO && icls == FP_ZERO)
    {
      res = x;
    }
  else
    {
      res = M_SUF (__kernel_casinh) (x, 0);
    }

  return res;
}

declare_mgen_alias (__casinh, casinh)



/* Return arc hyperbolic cosine for a complex type. */

CFLOAT
M_DECL_FUNC (__cacosh) (CFLOAT x)
{
  CFLOAT res;
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE)
    {
      if (icls == FP_INFINITE)
        {
          __real__ res = M_HUGE_VAL;

          if (rcls == FP_NAN)
            __imag__ res = M_NAN;
          else
            __imag__ res = M_COPYSIGN ((rcls == FP_INFINITE
                                        ? (__real__ x < 0
                                           ? M_MLIT (M_PI) - M_MLIT (M_PI_4)
                                           : M_MLIT (M_PI_4))
                                        : M_MLIT (M_PI_2)), __imag__ x);
        }
      else if (rcls == FP_INFINITE)
        {
          __real__ res = M_HUGE_VAL;

          if (icls >= FP_ZERO)
            __imag__ res = M_COPYSIGN (signbit (__real__ x)
                                       ? M_MLIT (M_PI) : 0, __imag__ x);
          else
            __imag__ res = M_NAN;
        }
      else
        {
          __real__ res = M_NAN;
          if (rcls == FP_ZERO)
            __imag__ res = M_MLIT (M_PI_2);
          else
            __imag__ res = M_NAN;
        }
    }
  else if (rcls == FP_ZERO && icls == FP_ZERO)
    {
      __real__ res = 0;
      __imag__ res = M_COPYSIGN (M_MLIT (M_PI_2), __imag__ x);
    }
  else
    {
      CFLOAT y;

      __real__ y = -__imag__ x;
      __imag__ y = __real__ x;

      y = M_SUF (__kernel_casinh) (y, 1);

      if (signbit (__imag__ x))
        {
          __real__ res = __real__ y;
          __imag__ res = -__imag__ y;
        }
      else
        {
          __real__ res = -__real__ y;
          __imag__ res = __imag__ y;
        }
    }

  return res;
}

declare_mgen_alias (__cacosh, cacosh)



/* Return arc hyperbolic tangent for a complex float type. */

CFLOAT
M_DECL_FUNC (__catanh) (CFLOAT x)
{
  CFLOAT res;
  int rcls = fpclassify (__real__ x);
  int icls = fpclassify (__imag__ x);

  if (__glibc_unlikely (rcls <= FP_INFINITE || icls <= FP_INFINITE))
    {
      if (icls == FP_INFINITE)
        {
          __real__ res = M_COPYSIGN (0, __real__ x);
          __imag__ res = M_COPYSIGN (M_MLIT (M_PI_2), __imag__ x);
        }
      else if (rcls == FP_INFINITE || rcls == FP_ZERO)
        {
          __real__ res = M_COPYSIGN (0, __real__ x);
          if (icls >= FP_ZERO)
            __imag__ res = M_COPYSIGN (M_MLIT (M_PI_2), __imag__ x);
          else
            __imag__ res = M_NAN;
        }
      else
        {
          __real__ res = M_NAN;
          __imag__ res = M_NAN;
        }
    }
  else if (__glibc_unlikely (rcls == FP_ZERO && icls == FP_ZERO))
    {
      res = x;
    }
  else
    {
      if (M_FABS (__real__ x) >= 16 / M_EPSILON
          || M_FABS (__imag__ x) >= 16 / M_EPSILON)
        {
          __imag__ res = M_COPYSIGN (M_MLIT (M_PI_2), __imag__ x);
          if (M_FABS (__imag__ x) <= 1)
            __real__ res = 1 / __real__ x;
          else if (M_FABS (__real__ x) <= 1)
            __real__ res = __real__ x / __imag__ x / __imag__ x;
          else
            {
              FLOAT h = M_HYPOT (__real__ x / 2, __imag__ x / 2);
              __real__ res = __real__ x / h / h / 4;
            }
        }
      else
        {
          if (M_FABS (__real__ x) == 1
              && M_FABS (__imag__ x) < M_EPSILON * M_EPSILON)
            __real__ res = (M_COPYSIGN (M_LIT (0.5), __real__ x)
                            * (M_MLIT (M_LN2)
                               - M_LOG (M_FABS (__imag__ x))));
          else
            {
              FLOAT i2 = 0;
              if (M_FABS (__imag__ x) >= M_EPSILON * M_EPSILON)
                i2 = __imag__ x * __imag__ x;

              FLOAT num = 1 + __real__ x;
              num = i2 + num * num;

              FLOAT den = 1 - __real__ x;
              den = i2 + den * den;

              FLOAT f = num / den;
              if (f < M_LIT (0.5))
                __real__ res = M_LIT (0.25) * M_LOG (f);
              else
                {
                  num = 4 * __real__ x;
                  __real__ res = M_LIT (0.25) * M_LOG1P (num / den);
                }
            }

          FLOAT absx, absy, den;

          absx = M_FABS (__real__ x);
          absy = M_FABS (__imag__ x);
          if (absx < absy)
            {
              FLOAT t = absx;
              absx = absy;
              absy = t;
            }

          if (absy < M_EPSILON / 2)
            {
              den = (1 - absx) * (1 + absx);
              if (den == 0)
                den = 0;
            }
          else if (absx >= 1)
            den = (1 - absx) * (1 + absx) - absy * absy;
          else if (absx >= M_LIT (0.75) || absy >= M_LIT (0.5))
            den = -M_SUF (__x2y2m1) (absx, absy);
          else
            den = (1 - absx) * (1 + absx) - absy * absy;

          __imag__ res = M_LIT (0.5) * M_ATAN2 (2 * __imag__ x, den);
        }

      math_check_force_underflow_complex (res);
    }

  return res;
}

declare_mgen_alias (__catanh, catanh)
