(() => {

// BEGIN POLYFILL //////////////////////////////////////////////////////////////

var M_E          = 2.71828182845904523536; // e
var M_LOG2E      = 1.44269504088896340736; // log2(e)
var M_LOG10E     = 0.43429448190325182765; // log10(e)
var M_LN2        = 0.69314718055994530942; // ln(2)
var M_LN10       = 2.30258509299404568402; // ln(10)
var M_LNPI       = 1.1447298858494002    ; // ln(pi)
var M_PI         = 3.14159265358979323846; // pi
var M_2PI        = 6.28318530717958647692; // 2*pi
var M_PI_2       = 1.57079632679489661923; // pi/2
var M_3PI_2      = 4.71238898038468985779; // 3*pi/2
var M_PI_4       = 0.78539816339744830962; // pi/4
var M_1_PI       = 0.31830988618379067154; // 1/pi
var M_2_PI       = 0.63661977236758134308; // 2/pi
var M_2_SQRTPI   = 1.12837916709551257390; // 2/sqrt(pi)
var M_SQRT2PI    = 2.5066282746310002    ; // sqrt(2*pi)
var M_SQRT2      = 1.41421356237309504880; // sqrt(2)
var M_SQRT1_2    = 0.70710678118654752440; // 1/sqrt(2)

var M_MANT_DIG   = 53;
var M_MIN_EXP    = -1021;                   // 2**-1021 underflows
var M_MAX_EXP    = 1024;                    // 2**1024 overflows
var M_EPSILON    = 2.220446049250313e-16;
var M_MIN        = 2.2250738585072014e-308; // smallest normal
var M_MAX        = 1.7976931348623157e+308;
var M_NAN        = NaN;
var M_HUGE_VAL   = Infinity;

var INFINITY     = Infinity;

var FE_INVALID   = -1;
var FP_NAN       = 0;
var FP_INFINITE  = 1;
var FP_ZERO      = 2;
var FP_SUBNORMAL = 3;
var FP_NORMAL    = 4;

function fpclassify(x) {
  if (Number.isNaN(x)) return FP_NAN;
  if (x === Infinity) return FP_INFINITE;
  if (x === -Infinity) return FP_INFINITE;
  if (x === 0) return FP_ZERO;
  if (x > 0 && x < M_MIN) return FP_SUBNORMAL;
  if (x < 0 && x > -M_MIN) return FP_SUBNORMAL;
  return FP_NORMAL;
}

function feraiseexcept() {} // not relevant to javascript

function math_check_force_underflow_complex() {} // same here

function signbit(x) {
  if (x === 0) return 1/x > 0 ? 0 : 1;
  return x > 0 ? 0 : 1;
}

function isinf(x) {
  return x === Infinity || x === -Infinity;
}

var isnan = Number.isNaN;

function isfinite(x) {
  return x !== Infinity && x !== -Infinity;
}

// need to optimize
function kernel_casinh(z,i) {
  var I = c(0,1);
  var asinhz = clog(cadd(z,csqrt(cmul(cadd(z,I),csub(z,I)))));
  if (i) buffer[asinhz+1] = M_PI_2 - buffer[asinhz+1];
  return asinhz;
}

function M_HYPOT(x,y) {
  x = Math.abs(x);
  y = Math.abs(y);
  if (x < y) {var t = x; x = y; y = t;}
  if (x < 1e154) return Math.sqrt(x*x + y*y);
  y /= x;
  return x*Math.sqrt(1 + y*y);
}

var M_ATAN2 = Math.atan2;

var M_SIN = Math.sin;

var M_COS = Math.cos;

var M_SINH = Math.sinh;

var M_COSH = Math.cosh;

function M_COPYSIGN(x, y) { // looks messy but is very fast and correct
  if (x === 0) { if (1/x > 0) {
    if (y === 0) return 1/y > 0 ? x : -x; else return y > 0 ? x : -x;
  } else {
    if (y === 0) return 1/y > 0 ? -x : x; else return y > 0 ? -x : x;
  }} else { if (x > 0) {
    if (y === 0) return 1/y > 0 ? x : -x; else return y > 0 ? x : -x;
  } else {
    if (y === 0) return 1/y > 0 ? -x : x; else return y > 0 ? -x : x;
  }}
}

function M_FABS(x) { // beats Math.abs
  if (x > 0) return x;
  return -x;
}

var M_SQRT = Math.sqrt;

function M_X2Y2M1(x, y) {
  return (x - 1) * (x + 1) + y * y;
}

var M_LOG = Math.log;

var M_LOG1P = Math.log1p;

var M_EXP = Math.exp;

var M_SCALBN = (() => { // not sure of a better way to do this
  var POWLUT = (new Float64Array(2098)).map((_,i) => Math.pow(2, i - 1074));
  function M_SCALBN(x,n) {
    // if (x === 0 || !Number.isFinite(x)) return x;
    // if (n > 1023) return x*Infinity;
    // if (n < -1074) return x*0;
    return x*POWLUT[n + 1074];
  }
  return M_SCALBN;
})();

var M_LOG10 = Math.log10;

// END POLYFILL ////////////////////////////////////////////////////////////////

// BEGIN BUFFER INITILIZATION & CONFIGURATION //////////////////////////////////

var bufferSize = 10_000_000;
var buffer = new Float64Array(bufferSize);
var index = 0;

function resetBuffer() {
  index = 0;
  return index;
}

function setBufferSize(size) {
  bufferSize = size;
  buffer = new Float64Array(bufferSize);
  index = 0;
}

// END BUFFER INITILIZATION & CONFIGURATION ////////////////////////////////////

// BEGIN STANDARD OPERATION IMPLEMENTATIONS ////////////////////////////////////
// * NOTE: THE FUNCTIONS WITHIN THIS SECION ARE NOT A PART OF GLIBC BUT ARE
// * NONETHELESS INVALUALBLE FOR INSPECTION PURPOSES AND DUE TO THE LACK OF
// * OPERATOR OVERLOADING IN JAVASCRIPT. BE MINDFUL THAT THESE FUNCTIONS AS
// * WRITTEN MAY NOT FULLY MIMIC C BEHAVIOR IN ALL EDGE CASES.

/* Compute complex string representation. */
// string ctoString(complex z)
function ctoString(z) {
  var re = buffer[z];
  var im = buffer[z+1];
  if (Number.isNaN(re) || Number.isNaN(im)) return 'NaN';
  if (re === 0) {
    if (im === 0) {
      if (1/re < 0) return 1/im < 0 ? '-0 - 0i' : '-0';
      else return 1/im < 0 ? '0 - 0i' : '0';
    }
  } else if (im === 0) {
    if (1/im < 0) return re.toString() + ' - 0i';
    return re.toString();
  } else if (re === 0) {
    if (1/re < 0) {
      if (im === 1) return '-0 - i';
      else if (im === -1) return '-0 - i';
      else if (im === Infinity) return '-0 + Infinity*i';
      else if (im === -Infinity) return '-0 - Infinity*i';
      else if (im > 0) return '-0 + ' + im.toString() + 'i';
      else return '-0 - ' + Math.abs(im) + 'i';
    }
    if (im === 1) return 'i';
    else if (im === -1) return '-i';
    else if (!Number.isFinite(im)) return im.toString() + '*i';
    else return im.toString() + 'i';
  }
  if (im === 1) return re.toString() + ' + i';
  else if (im === -1) return re.toString() + ' + i';
  else if (im === Infinity) return re.toString() + ' + Infinity*i';
  else if (im === -Infinity) return re.toString() + ' - Infinity*i';
  else if (im > 0)
    return re.toString() + ' + ' + im.toString() + 'i';
  else if (im < 0)
    return re.toString() + ' - ' + (-im).toString() + 'i';
  return 'NaN';
}

/* Create complex number. */
// complex c(float x, float y)
function c(x,y) {
  var retval = index;
  index += 2;
  buffer[retval] = x;
  buffer[retval+1] = y;
  return retval;
}

/* Compute complex negation. */
// complex cneg(complex z)
function cneg(z) {
  var retval = index;
  index += 2;
  buffer[retval] = -buffer[z];
  buffer[retval+1] = -buffer[z+1];
  return retval;
}

/* Compute complex sign. */
// complex csgn(complex z)
function csgn(z) {
  var retval = index;
  index += 2;
  var absz = M_HYPOT(buffer[z], buffer[z+1]);
  if (absz === 0) {
    buffer[retval] = M_COPYSIGN(0,buffer[z]);
    buffer[retval+1] = M_COPYSIGN(0,buffer[z+1]);
    return retval;
  }
  var dinv = 1/absz;
  buffer[retval] = buffer[z]*dinv;
  buffer[retval+1] = buffer[z+1]*dinv;
  return retval;
}

/* Compute complex reciprocal. */
// complex cinv(complex z)
function cinv(z) {
  var retval = index;
  index += 2;
  var dinv = 1/(buffer[z]*buffer[z] + buffer[z+1]*buffer[z+1]);
  buffer[retval] = buffer[z]*dinv;
  buffer[retval+1] = -buffer[z+1]*dinv;
  return retval;
}

/* Compute complex addition. */
// complex cadd(complex z, complex w)
function cadd(z, w) {
  var retval = index;
  index += 2;
  buffer[retval] = buffer[z] + buffer[w];
  buffer[retval+1] = buffer[z+1] + buffer[w+1];
  return retval;
}

/* Compute complex subtraction. */
// complex csub(complex z, complex w)
function csub(z, w) {
  var retval = index;
  index += 2;
  buffer[retval] = buffer[z] - buffer[w];
  buffer[retval+1] = buffer[z+1] - buffer[w+1];
  return retval;
}

/* Compute complex multiplication. */
// complex cmul(complex z, complex w)
function cmul(z, w) {
  var retval = index;
  index += 2;
  buffer[retval] = buffer[z]*buffer[w] - buffer[z+1]*buffer[w+1];
  buffer[retval+1] = buffer[z]*buffer[w+1] + buffer[z+1]*buffer[w];
  return retval;
}

/* Compute complex division. */
// complex cdiv(complex z, complex w)
function cdiv(z, w) {
  var retval = index;
  index += 2;
  var dinv = 1/(buffer[w]*buffer[w] + buffer[w+1]*buffer[w+1]);
  buffer[retval] = (buffer[z]*buffer[w] - buffer[z+1]*buffer[w+1])*dinv;
  buffer[retval+1] = (buffer[z]*buffer[w+1] + buffer[z+1]*buffer[w])*dinv;
  return retval;
}

/* Compute complex scalar addition */
// complex caddScalar(complex z, float n)
function caddScalar(z, n) {
  var retval = index;
  index += 2;
  buffer[retval] = buffer[z] + n;
  buffer[retval+1] = buffer[z+1];
  return retval;
}

/* Compute complex scalar subtraction */
// complex csubScalar(complex z, float n)
function csubScalar(z, n) {
  var retval = index;
  index += 2;
  buffer[retval] = buffer[z] - n;
  buffer[retval+1] = buffer[z+1];
  return retval;
}

/* Compute complex scalar multiplication */
// complex cmulScalar(complex z, float n)
function cmulScalar(z, n) {
  var retval = index;
  index += 2;
  buffer[retval] = buffer[z] * n;
  buffer[retval+1] = buffer[z+1] * n;
  return retval;
}

/* Compute complex scalar division */
// complex cdivScalar(complex z, float n)
function cdivScalar(z, n) {
  var retval = index;
  index += 2;
  n = 1/n;
  buffer[retval] = buffer[z] * n;
  buffer[retval+1] = buffer[z+1] * n;
  return retval;
}

/* Compute complex imaginary multiplication */
// complex cmulScalar(complex z, float n)
function cmuli(z,n=1) {
  var retval = index;
  index += 2;
  buffer[retval] = -buffer[z+1] * n;
  buffer[retval+1] = buffer[z] * n;
  return retval;
}

// END STANDARD OPERATION IMPLEMENTATIONS //////////////////////////////////////

// BEGIN GLIBC PORT ////////////////////////////////////////////////////////////

/* Return real part of complex float type. */
// float creal(complex z)
function creal(z) {
  return buffer[z];
}

/* Return imaginary part of complex float type. */
// float cimag(complex z)
function cimag(z) {
  return buffer[z+1];
}

/* Return the complex absolute value of complex float type. */
// float cabs(complex z)
function cabs(z) {
  return M_HYPOT(buffer[z], buffer[z+1]);
}

/* Compute argument of complex float type. */
// float carg(complex x)
function carg(x) {
  return M_ATAN2(buffer[x+1], buffer[x]);
}

/* Return complex conjugate of complex float type. */
// complex conj(complex z)
function conj(z) {
  var retval = index;
  index += 2;
  buffer[retval] = buffer[z];
  buffer[retval+1] = -buffer[z+1];
  return retval;
}

/* Compute projection of complex float type value to Riemann sphere. */
// complex cproj(complex x)
function cproj(x) {
  if (isinf(buffer[x]) || isinf(buffer[x+1])) {
    var res = index;
    index += 2;

    buffer[res] = INFINITY;
    buffer[res+1] = M_COPYSIGN(0, buffer[x+1]);

    return res;
  }

  return x;
}

/* Return value of complex exponential function for a float type. */
// complex cexp(complex x)
function cexp(x) {
  var retval = index;
  index += 2;
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  if (rcls >= FP_ZERO) {
    /* Real part is finite.  */
    if (icls >= FP_ZERO) {
      /* Imaginary part is finite.  */
      var t = ((M_MAX_EXP - 1) * M_LN2) | 0;
      var sinix, cosix;

      if (M_FABS(buffer[x+1]) > M_MIN) {
        sinix = M_SIN(buffer[x+1]);
        cosix = M_COS(buffer[x+1]);
      } else {
        sinix = buffer[x+1];
        cosix = 1;
      }

      if (buffer[x] > t) {
        var exp_t = M_EXP(t);
        buffer[x] -= t;
        sinix *= exp_t;
        cosix *= exp_t;
        if (buffer[x] > t) {
          buffer[x] -= t;
          sinix *= exp_t;
          cosix *= exp_t;
        }
      }
      if (buffer[x] > t) {
        /* Overflow (original real part of x > 3t).  */
        buffer[retval] = M_MAX * cosix;
        buffer[retval+1] = M_MAX * sinix;
      } else {
        var exp_val = M_EXP(buffer[x]);
        buffer[retval] = exp_val * cosix;
        buffer[retval+1] = exp_val * sinix;
      }
      math_check_force_underflow_complex(retval);
    } else {
      /* If the imaginary part is +-inf or NaN and the real part
         is not +-inf the result is NaN + iNaN.  */
      buffer[retval] = M_NAN;
      buffer[retval+1] = M_NAN;

      feraiseexcept(FE_INVALID);
    }
  } else if (rcls == FP_INFINITE) {
    /* Real part is infinite.  */
    if (icls >= FP_ZERO) {
      /* Imaginary part is finite.  */
      var value = signbit(buffer[x]) ? 0 : M_HUGE_VAL;

      if (icls == FP_ZERO) {
        /* Imaginary part is 0.0.  */
        buffer[retval] = value;
        buffer[retval+1] = buffer[x+1];
      } else {
        var sinix, cosix;

        if (M_FABS(buffer[x+1]) > M_MIN) {
          sinix = M_SIN(buffer[x+1]);
          cosix = M_COS(buffer[x+1]);
        } else {
          sinix = buffer[x+1];
          cosix = 1;
        }

        buffer[retval] = M_COPYSIGN(value, cosix);
        buffer[retval+1] = M_COPYSIGN(value, sinix);
      }
    } else if (signbit(buffer[x]) == 0) {
      buffer[retval] = M_HUGE_VAL;
      buffer[retval+1] = buffer[x+1] - buffer[x+1];
    } else {
      buffer[retval] = 0;
      buffer[retval+1] = M_COPYSIGN(0, buffer[x+1]);
    }
  } else {
    /* If the real part is NaN the result is NaN + iNaN unless the
       imaginary part is zero.  */
    buffer[retval] = M_NAN;
    if (icls == FP_ZERO)
      buffer[retval+1] = buffer[x+1];
    else {
      buffer[retval+1] = M_NAN;

      if (rcls != FP_NAN || icls != FP_NAN)
        feraiseexcept(FE_INVALID);
    }
  }

  return retval;
}

/* Compute complex natural logarithm. */
// complex clog(complex x)
function clog(x) {
  var result = index;
  index += 2;
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  if (rcls == FP_ZERO && icls == FP_ZERO) {
    /* Real and imaginary part are 0.0.  */
    buffer[result+1] = signbit(buffer[x]) ? M_PI : 0;
    buffer[result+1] = M_COPYSIGN(buffer[result+1], buffer[x+1]);
    /* Yes, the following line raises an exception.  */
    buffer[result] = -1 / M_FABS(buffer[x]);
  } else if (rcls != FP_NAN && icls != FP_NAN) {
    /* Neither real nor imaginary part is NaN.  */
    var absx = M_FABS(buffer[x]), absy = M_FABS(buffer[x+1]);
    var scale = 0;

    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }

    if (absx > M_MAX / 2) {
      scale = -1;
      absx = M_SCALBN(absx, scale);
      absy = (absy >= M_MIN * 2 ? M_SCALBN(absy, scale) : 0);
    } else if (absx < M_MIN && absy < M_MIN) {
      scale = M_MANT_DIG;
      absx = M_SCALBN(absx, scale);
      absy = M_SCALBN(absy, scale);
    }

    if (absx == 1 && scale == 0) {
      buffer[result] = M_LOG1P(absy * absy) / 2;
      // math_check_force_underflow_nonneg (buffer[result]);
    } else if (absx > 1 && absx < 2 && absy < 1 && scale == 0) {
      var d2m1 = (absx - 1) * (absx + 1);
      if (absy >= M_EPSILON)
        d2m1 += absy * absy;
      buffer[result] = M_LOG1P(d2m1) / 2;
    } else if (absx < 1
            && absx >= 0.5
            && absy < M_EPSILON / 2
            && scale == 0) {
      var d2m1 = (absx - 1) * (absx + 1);
      buffer[result] = M_LOG1P(d2m1) / 2;
    } else if (absx < 1
            && absx >= 0.5
            && scale == 0
            && absx * absx + absy * absy >= 0.5) {
      var d2m1 = M_X2Y2M1(absx, absy);
      buffer[result] = M_LOG1P(d2m1) / 2;
    } else {
      var d = M_HYPOT(absx, absy);
      buffer[result] = M_LOG(d) - scale * M_LN2;
    }

    buffer[result+1] = M_ATAN2(buffer[x+1], buffer[x]);
  } else {
    buffer[result+1] = M_NAN;
    if (rcls == FP_INFINITE || icls == FP_INFINITE)
      /* Real or imaginary part is infinite.  */
      buffer[result] = M_HUGE_VAL;
    else
      buffer[result] = M_NAN;
  }

  return result;
}

/* Compute complex base 10 logarithm. */
/* log_10 (2).  */
var LOG10_2 = 0.3010299956639811952137388947244930267682;
/* pi * log10 (e).  */
var PI_LOG10E = 1.364376353841841347485783625431355770210;
// complex clog10(complex x)
function clog10(x) {
  var result = index;
  index += 2;
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  if (rcls == FP_ZERO && icls == FP_ZERO) {
    /* Real and imaginary part are 0.0.  */
    buffer[result+1] = signbit(buffer[x]) ? PI_LOG10E : 0;
    buffer[result+1] = M_COPYSIGN(buffer[result+1], buffer[x+1]);
    /* Yes, the following line raises an exception.  */
    buffer[result] = -1 / M_FABS(buffer[x]);
  } else if (rcls != FP_NAN && icls != FP_NAN) {
    /* Neither real nor imaginary part is NaN.  */
    var absx = M_FABS(buffer[x]), absy = M_FABS(buffer[x+1]);
    var scale = 0;

    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }

    if (absx > M_MAX / 2) {
      scale = -1;
      absx = M_SCALBN(absx, scale);
      absy = (absy >= M_MIN * 2 ? M_SCALBN(absy, scale) : 0);
    } else if (absx < M_MIN && absy < M_MIN) {
      scale = M_MANT_DIG;
      absx = M_SCALBN(absx, scale);
      absy = M_SCALBN(absy, scale);
    }

    if (absx == 1 && scale == 0) {
      buffer[result] = (M_LOG1P(absy * absy)
                     * (M_LOG10E / 2));
      // math_check_force_underflow_nonneg (buffer[result]);
    } else if (absx > 1 && absx < 2 && absy < 1 && scale == 0) {
      var d2m1 = (absx - 1) * (absx + 1);
      if (absy >= M_EPSILON)
        d2m1 += absy * absy;
      buffer[result] = M_LOG1P(d2m1) * (M_LOG10E / 2);
    } else if (absx < 1
            && absx >= 0.5
            && absy < M_EPSILON / 2
            && scale == 0) {
      var d2m1 = (absx - 1) * (absx + 1);
      buffer[result] = M_LOG1P(d2m1) * (M_LOG10E / 2);
    } else if (absx < 1
            && absx >= 0.5
            && scale == 0
            && absx * absx + absy * absy >= 0.5) {
      var d2m1 = M_X2Y2M1(absx, absy);
      buffer[result] = M_LOG1P(d2m1) * (M_LOG10E / 2);
    } else {
      var d = M_HYPOT(absx, absy);
      buffer[result] = M_LOG10(d) - scale * LOG10_2;
    }

    buffer[result+1] = M_LOG10E * M_ATAN2(buffer[x+1], buffer[x]);
  } else {
    buffer[result+1] = M_NAN;
    if (rcls == FP_INFINITE || icls == FP_INFINITE)
      /* Real or imaginary part is infinite.  */
      buffer[result] = M_HUGE_VAL;
    else
      buffer[result] = M_NAN;
  }

  return result;
}

/* Complex power of float type. */
// complex cpow(complex x, complex c)
function cpow(x, c) {
  return cexp(cmul(c, clog(x)));
}

/* Complex square root of a float type. */
// complex csqrt(complex x)
function csqrt(x) {
  var res = index;
  index += 2;
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      buffer[res] = M_HUGE_VAL;
      buffer[res+1] = buffer[x+1];
    } else if (rcls == FP_INFINITE) {
      if (buffer[x] < 0) {
        buffer[res] = icls == FP_NAN ? M_NAN : 0;
        buffer[res+1] = M_COPYSIGN(M_HUGE_VAL, buffer[x+1]);
      } else {
        buffer[res] = buffer[x];
        buffer[res+1] = (icls == FP_NAN
                      ? M_NAN : M_COPYSIGN(0, buffer[x+1]));
      }
    } else {
      buffer[res] = M_NAN;
      buffer[res+1] = M_NAN;
    }
  } else {
    if (icls == FP_ZERO) {
      if (buffer[x] < 0) {
        buffer[res] = 0;
        buffer[res+1] = M_COPYSIGN(M_SQRT(-buffer[x]), buffer[x+1]);
      } else {
        buffer[res] = M_FABS(M_SQRT(buffer[x]));
        buffer[res+1] = M_COPYSIGN(0, buffer[x+1]);
      }
    } else if (rcls == FP_ZERO) {
      var r;
      if (M_FABS(buffer[x+1]) >= 2 * M_MIN)
        r = M_SQRT(0.5 * M_FABS(buffer[x+1]));
      else
        r = 0.5 * M_SQRT(2 * M_FABS(buffer[x+1]));

      buffer[res] = r;
      buffer[res+1] = M_COPYSIGN(r, buffer[x+1]);
    } else {
      var d, r, s;
      var scale = 0;

      if (M_FABS(buffer[x]) > M_MAX / 4) {
        scale = 1;
        buffer[x] = M_SCALBN(buffer[x], -2 * scale);
        buffer[x+1] = M_SCALBN(buffer[x+1], -2 * scale);
      } else if (M_FABS(buffer[x+1]) > M_MAX / 4) {
        scale = 1;
        if (M_FABS(buffer[x]) >= 4 * M_MIN)
          buffer[x] = M_SCALBN(buffer[x], -2 * scale);
        else
          buffer[x] = 0;
        buffer[x+1] = M_SCALBN(buffer[x+1], -2 * scale);
      } else if (M_FABS(buffer[x]) < 2 * M_MIN
              && M_FABS(buffer[x+1]) < 2 * M_MIN) {
        scale = -((M_MANT_DIG + 1) / 2);
        buffer[x] = M_SCALBN(buffer[x], -2 * scale);
        buffer[x+1] = M_SCALBN(buffer[x+1], -2 * scale);
      }

      d = M_HYPOT(buffer[x], buffer[x+1]);
      /* Use the identity   2  Re res  Im res = Im x
         to avoid cancellation error in  d +/- Re x.  */
      if (buffer[x] > 0) {
        r = M_SQRT(0.5 * (d + buffer[x]));
        if (scale == 1 && M_FABS(buffer[x+1]) < 1) {
          /* Avoid possible intermediate underflow.  */
          s = buffer[x+1] / r;
          r = M_SCALBN(r, scale);
          scale = 0;
        } else
        s = 0.5 * (buffer[x+1] / r);
      } else {
        s = M_SQRT(0.5 * (d - buffer[x]));
        if (scale == 1 && M_FABS(buffer[x+1]) < 1) {
          /* Avoid possible intermediate underflow.  */
          r = M_FABS(buffer[x+1] / s);
          s = M_SCALBN(s, scale);
          scale = 0;
        } else
        r = M_FABS(0.5 * (buffer[x+1] / s));
      }

      if (scale) {
        r = M_SCALBN(r, scale);
        s = M_SCALBN(s, scale);
      }

      math_check_force_underflow(r);
      math_check_force_underflow(s);

      buffer[res] = r;
      buffer[res+1] = M_COPYSIGN(s, buffer[x+1]);
    }
  }

  return res;
}

/* Complex sine function for float types. */
// complex csin(complex x)
function csin(x) {
  var retval = index;
  index += 2;
  var negate = signbit(buffer[x]);
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  buffer[x] = M_FABS(buffer[x]);

  if (icls >= FP_ZERO) {
    /* Imaginary part is finite.  */
    if (rcls >= FP_ZERO) {
      /* Real part is finite.  */
      var t = ((M_MAX_EXP - 1) * M_LN2) | 0;
      var sinix, cosix;

      if (buffer[x] > M_MIN) {
        sinix = M_SIN(buffer[x]);
        cosix = M_COS(buffer[x]);
      } else {
        sinix = buffer[x];
        cosix = 1;
      }

      if (negate)
        sinix = -sinix;

      if (M_FABS(buffer[x+1]) > t) {
        var exp_t = M_EXP(t);
        var ix = M_FABS(buffer[x+1]);
        if (signbit(buffer[x+1]))
          cosix = -cosix;
        ix -= t;
        sinix *= exp_t / 2;
        cosix *= exp_t / 2;
        if (ix > t) {
          ix -= t;
          sinix *= exp_t;
          cosix *= exp_t;
        }
        if (ix > t) {
          /* Overflow (original imaginary part of x > 3t).  */
          buffer[retval] = M_MAX * sinix;
          buffer[retval+1] = M_MAX * cosix;
        } else {
          var exp_val = M_EXP(ix);
          buffer[retval] = exp_val * sinix;
          buffer[retval+1] = exp_val * cosix;
        }
      } else {
        buffer[retval] = M_COSH(buffer[x+1]) * sinix;
        buffer[retval+1] = M_SINH(buffer[x+1]) * cosix;
      }

      math_check_force_underflow_complex(retval);
    } else {
      if (icls == FP_ZERO) {
        /* Imaginary part is 0.0.  */
        buffer[retval] = buffer[x] - buffer[x];
        buffer[retval+1] = buffer[x+1];
      } else {
        buffer[retval] = M_NAN;
        buffer[retval+1] = M_NAN;

        feraiseexcept(FE_INVALID);
      }
    }
  } else if (icls == FP_INFINITE) {
    /* Imaginary part is infinite.  */
    if (rcls == FP_ZERO) {
      /* Real part is 0.0.  */
      buffer[retval] = M_COPYSIGN(0, negate ? -1 : 1);
      buffer[retval+1] = buffer[x+1];
    } else if (rcls > FP_ZERO) {
      /* Real part is finite.  */
      var sinix, cosix;

      if (buffer[x] > M_MIN) {
        sinix = M_SIN(buffer[x]);
        cosix = M_COS(buffer[x]);
      } else {
        sinix = buffer[x];
        cosix = 1;
      }

      buffer[retval] = M_COPYSIGN(M_HUGE_VAL, sinix);
      buffer[retval+1] = M_COPYSIGN(M_HUGE_VAL, cosix);

      if (negate)
        buffer[retval] = -buffer[retval];
      if (signbit(buffer[x+1]))
        buffer[retval+1] = -buffer[retval+1];
    } else {
      buffer[retval] = buffer[x] - buffer[x];
      buffer[retval+1] = M_HUGE_VAL;
    }
  } else {
    if (rcls == FP_ZERO)
      buffer[retval] = M_COPYSIGN(0, negate ? -1 : 1);
    else
      buffer[retval] = M_NAN;
    buffer[retval+1] = M_NAN;
  }

  return retval;
}

/* Return cosine of complex float type. */
// complex ccos(complex x)
function ccos(x) {
  var y = index;
  index += 2;

  buffer[y] = -buffer[x+1];
  buffer[y+1] = buffer[x];

  return ccosh(y);
}

/* Complex tangent function for a complex float type. */
// complex ctan(complex x)
function ctan(x) {
  var res = index;
  index += 2;

  if (!isfinite(buffer[x]) || !isfinite(buffer[x+1])) {
    if (isinf(buffer[x+1])) {
      if (isfinite(buffer[x]) && M_FABS(buffer[x]) > 1) {
        var sinrx, cosrx;
        sinrx = M_SIN(buffer[x]);
        cosrx = M_COS(buffer[x]);
        buffer[res] = M_COPYSIGN(0, sinrx * cosrx);
      } else
        buffer[res] = M_COPYSIGN(0, buffer[x]);
      buffer[res+1] = M_COPYSIGN(1, buffer[x+1]);
    } else if (buffer[x] == 0) {
      res = x;
    } else {
      buffer[res] = M_NAN;
      if (buffer[x+1] == 0)
        buffer[res+1] = buffer[x+1];
      else
        buffer[res+1] = M_NAN;

      if (isinf(buffer[x]))
        feraiseexcept(FE_INVALID);
    }
  } else {
    var sinrx, cosrx;
    var den;
    var t = ((M_MAX_EXP - 1) * M_LN2 / 2) | 0;

    /* tan(x+iy) = (sin(2x) + i*sinh(2y))/(cos(2x) + cosh(2y))
       = (sin(x)*cos(x) + i*sinh(y)*cosh(y)/(cos(x)^2 + sinh(y)^2). */

    if (M_FABS(buffer[x]) > M_MIN) {
      sinrx = M_SIN(buffer[x]);
      cosrx = M_COS(buffer[x]);
    } else {
      sinrx = buffer[x];
      cosrx = 1;
    }

    if (M_FABS(buffer[x+1]) > t) {
      /* Avoid intermediate overflow when the real part of the
         result may be subnormal.  Ignoring negligible terms, the
         imaginary part is +/- 1, the real part is
         sin(x)*cos(x)/sinh(y)^2 = 4*sin(x)*cos(x)/exp(2y).  */
      var exp_2t = M_EXP(2 * t);

      buffer[res+1] = M_COPYSIGN(1, buffer[x+1]);
      buffer[res] = 4 * sinrx * cosrx;
      buffer[x+1] = M_FABS(buffer[x+1]);
      buffer[x+1] -= t;
      buffer[res] /= exp_2t;
      if (buffer[x+1] > t) {
        /* Underflow (original imaginary part of x has absolute
           value > 2t).  */
        buffer[res] /= exp_2t;
      } else
        buffer[res] /= M_EXP(2 * buffer[x+1]);
    } else {
      var sinhix, coshix;
      if (M_FABS(buffer[x+1]) > M_MIN) {
        sinhix = M_SINH(buffer[x+1]);
        coshix = M_COSH(buffer[x+1]);
      } else {
        sinhix = buffer[x+1];
        coshix = 1;
      }

      if (M_FABS(sinhix) > M_FABS(cosrx) * M_EPSILON)
        den = cosrx * cosrx + sinhix * sinhix;
      else
        den = cosrx * cosrx;
      buffer[res] = sinrx * cosrx / den;
      buffer[res+1] = sinhix * coshix / den;
    }
    math_check_force_underflow_complex(res);
  }

  return res;
}

/* Complex sine hyperbole function for float types. */
// complex csinh(complex x)
function csinh(x) {
  var retval = index;
  index += 2;
  var negate = signbit(buffer[x]);
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  buffer[x] = M_FABS(buffer[x]);

  if (rcls >= FP_ZERO) {
    /* Real part is finite.  */
    if (icls >= FP_ZERO) {
      /* Imaginary part is finite.  */
      var t = ((M_MAX_EXP - 1) * M_LN2) | 0;
      var sinix, cosix;

      if (M_FABS(buffer[x+1]) > M_MIN) {
        sinix = M_SIN(buffer[x+1]);
        cosix = M_COS(buffer[x+1]);
      } else {
        sinix = buffer[x+1];
        cosix = 1;
      }

      if (negate)
        cosix = -cosix;

      if (M_FABS(buffer[x]) > t) {
        var exp_t = M_EXP(t);
        var rx = M_FABS(buffer[x]);
        if (signbit(buffer[x]))
          cosix = -cosix;
        rx -= t;
        sinix *= exp_t / 2;
        cosix *= exp_t / 2;
        if (rx > t) {
          rx -= t;
          sinix *= exp_t;
          cosix *= exp_t;
        }
        if (rx > t) {
          /* Overflow (original real part of x > 3t).  */
          buffer[retval] = M_MAX * cosix;
          buffer[retval+1] = M_MAX * sinix;
        } else {
          var exp_val = M_EXP(rx);
          buffer[retval] = exp_val * cosix;
          buffer[retval+1] = exp_val * sinix;
        }
      } else {
        buffer[retval] = M_SINH(buffer[x]) * cosix;
        buffer[retval+1] = M_COSH(buffer[x]) * sinix;
      }

      math_check_force_underflow_complex(retval);
    } else {
      if (rcls == FP_ZERO) {
        /* Real part is 0.0.  */
        buffer[retval] = M_COPYSIGN(0, negate ? -1 : 1);
        buffer[retval+1] = buffer[x+1] - buffer[x+1];
      } else {
        buffer[retval] = M_NAN;
        buffer[retval+1] = M_NAN;

        feraiseexcept(FE_INVALID);
      }
    }
  } else if (rcls == FP_INFINITE) {
    /* Real part is infinite.  */
    if (icls > FP_ZERO) {
      /* Imaginary part is finite.  */
      var sinix, cosix;

      if (M_FABS(buffer[x+1]) > M_MIN) {
        sinix = M_SIN(buffer[x+1]);
        cosix = M_COS(buffer[x+1]);
      } else {
        sinix = buffer[x+1];
        cosix = 1;
      }

      buffer[retval] = M_COPYSIGN(M_HUGE_VAL, cosix);
      buffer[retval+1] = M_COPYSIGN(M_HUGE_VAL, sinix);

      if (negate)
        buffer[retval] = -buffer[retval];
    } else if (icls == FP_ZERO) {
      /* Imaginary part is 0.0.  */
      buffer[retval] = negate ? -M_HUGE_VAL : M_HUGE_VAL;
      buffer[retval+1] = buffer[x+1];
    } else {
      buffer[retval] = M_HUGE_VAL;
      buffer[retval+1] = buffer[x+1] - buffer[x+1];
    }
  } else {
    buffer[retval] = M_NAN;
    buffer[retval+1] = buffer[x+1] == 0 ? buffer[x+1] : M_NAN;
  }

  return retval;
}

/* Complex cosine hyperbolic function for float types. */
// complex ccosh(complex x)
function ccosh(x) {
  var retval = index;
  index += 2;
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  if (rcls >= FP_ZERO) {
    /* Real part is finite.  */
    if (icls >= FP_ZERO) {
      /* Imaginary part is finite.  */
      var t = ((M_MAX_EXP - 1) * M_LN2) | 0;
      var sinix, cosix;

      if (M_FABS(buffer[x+1]) > M_MIN) {
        sinix = M_SIN(buffer[x+1]);
        cosix = M_COS(buffer[x+1]);
      } else {
        sinix = buffer[x+1];
        cosix = 1;
      }

      if (M_FABS(buffer[x]) > t) {
        var exp_t = M_EXP(t);
        var rx = M_FABS(buffer[x]);
        if (signbit(buffer[x]))
          sinix = -sinix;
        rx -= t;
        sinix *= exp_t / 2;
        cosix *= exp_t / 2;
        if (rx > t) {
          rx -= t;
          sinix *= exp_t;
          cosix *= exp_t;
        }
        if (rx > t) {
          /* Overflow (original real part of x > 3t).  */
          buffer[retval] = M_MAX * cosix;
          buffer[retval+1] = M_MAX * sinix;
        } else {
          var exp_val = M_EXP(rx);
          buffer[retval] = exp_val * cosix;
          buffer[retval+1] = exp_val * sinix;
        }
      } else {
        buffer[retval] = M_COSH(buffer[x]) * cosix;
        buffer[retval+1] = M_SINH(buffer[x]) * sinix;
      }

      math_check_force_underflow_complex(retval);
    } else {
      buffer[retval+1] = buffer[x] == 0 ? 0 : M_NAN;
      buffer[retval] = buffer[x+1] - buffer[x+1];
    }
  } else if (rcls == FP_INFINITE) {
    /* Real part is infinite.  */
    if (icls > FP_ZERO) {
      /* Imaginary part is finite.  */
      var sinix, cosix;

      if (M_FABS(buffer[x+1]) > M_MIN) {
        sinix = M_SIN(buffer[x+1]);
        cosix = M_COS(buffer[x+1]);
      } else {
        sinix = buffer[x+1];
        cosix = 1;
      }

      buffer[retval] = M_COPYSIGN(M_HUGE_VAL, cosix);
      buffer[retval+1] = (M_COPYSIGN(M_HUGE_VAL, sinix)
                         * M_COPYSIGN(1, buffer[x]));
    } else if (icls == FP_ZERO) {
      /* Imaginary part is 0.0.  */
      buffer[retval] = M_HUGE_VAL;
      buffer[retval+1] = buffer[x+1] * M_COPYSIGN(1, buffer[x]);
    } else {
      buffer[retval] = M_HUGE_VAL;
      buffer[retval+1] = buffer[x+1] - buffer[x+1];
    }
  } else {
    buffer[retval] = M_NAN;
    buffer[retval+1] = buffer[x+1] == 0 ? buffer[x+1] : M_NAN;
  }

  return retval;
}

/* Complex hyperbolic tangent for float types. */
// complex ctanh(complex x)
function ctanh(x) {
  var res = index;
  index += 2;

  if (!isfinite(buffer[x]) || !isfinite(buffer[x+1])) {
    if (isinf(buffer[x])) {
      buffer[res] = M_COPYSIGN(1, buffer[x]);
      if (isfinite(buffer[x+1]) && M_FABS(buffer[x+1]) > 1) {
        var sinix, cosix;
        sinix = M_SIN(buffer[x+1]);
        cosix = M_COS(buffer[x+1]);
        buffer[res+1] = M_COPYSIGN(0, sinix * cosix);
      } else
        buffer[res+1] = M_COPYSIGN(0, buffer[x+1]);
    } else if (buffer[x+1] == 0) {
      res = x;
    } else {
      if (buffer[x] == 0)
        buffer[res] = buffer[x];
      else
        buffer[res] = M_NAN;
      buffer[res+1] = M_NAN;

      if (isinf(buffer[x+1]))
        feraiseexcept(FE_INVALID);
    }
  } else {
    var sinix, cosix;
    var den;
    var t = ((M_MAX_EXP - 1) * M_LN2 / 2) | 0;

    /* tanh(x+iy) = (sinh(2x) + i*sin(2y))/(cosh(2x) + cos(2y))
       = (sinh(x)*cosh(x) + i*sin(y)*cos(y))/(sinh(x)^2 + cos(y)^2).  */

    if (M_FABS(buffer[x+1]) > M_MIN) {
      sinix = M_SIN(buffer[x+1]);
      cosix = M_COS(buffer[x+1]);
    } else {
      sinix = buffer[x+1];
      cosix = 1;
    }

    if (M_FABS(buffer[x]) > t) {
      /* Avoid intermediate overflow when the imaginary part of
         the result may be subnormal.  Ignoring negligible terms,
         the real part is +/- 1, the imaginary part is
         sin(y)*cos(y)/sinh(x)^2 = 4*sin(y)*cos(y)/exp(2x).  */
      var exp_2t = M_EXP(2 * t);

      buffer[res] = M_COPYSIGN(1, buffer[x]);
      buffer[res+1] = 4 * sinix * cosix;
      buffer[x] = M_FABS(buffer[x]);
      buffer[x] -= t;
      buffer[res+1] /= exp_2t;
      if (buffer[x] > t) {
        /* Underflow (original real part of x has absolute value
           > 2t).  */
        buffer[res+1] /= exp_2t;
      } else
        buffer[res+1] /= M_EXP(2 * buffer[x]);
    } else {
      var sinhrx, coshrx;
      if (M_FABS(buffer[x]) > M_MIN) {
        sinhrx = M_SINH(buffer[x]);
        coshrx = M_COSH(buffer[x]);
      } else {
        sinhrx = buffer[x];
        coshrx = 1;
      }

      if (M_FABS(sinhrx) > M_FABS(cosix) * M_EPSILON)
        den = sinhrx * sinhrx + cosix * cosix;
      else
        den = cosix * cosix;
      buffer[res] = sinhrx * coshrx / den;
      buffer[res+1] = sinix * cosix / den;
    }
    math_check_force_underflow_complex(res);
  }

  return res;
}

/* Return arc sine of a complex float type. */
// complex casin(complex x)
function casin(x) {
  var res = index;
  index += 2;

  if (isnan(buffer[x]) || isnan(buffer[x+1])) {
    if (buffer[x] == 0) {
      res = x;
    } else if (isinf(buffer[x]) || isinf(buffer[x+1])) {
      buffer[res] = M_NAN;
      buffer[res+1] = M_COPYSIGN(M_HUGE_VAL, buffer[x+1]);
    } else {
      buffer[res] = M_NAN;
      buffer[res+1] = M_NAN;
    }
  } else {
    var y = index;
    index += 2;

    buffer[y] = -buffer[x+1];
    buffer[y+1] = buffer[x];

    y = casinh(y);

    buffer[res] = buffer[y+1];
    buffer[res+1] = -buffer[y];
  }

  return res;
}

/* Return cosine of a complex type. */
// complex cacos(complex x)
function cacos(x) {
  var y = index;
  index += 2;
  var res = index;
  index += 2;
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE
    || (rcls == FP_ZERO && icls == FP_ZERO)) {
    y = casin(x);

    buffer[res] = M_PI_2 - buffer[y];
    if (buffer[res] == 0)
      buffer[res] = 0;
    buffer[res+1] = -buffer[y+1];
  } else {
    buffer[y] = -buffer[x+1];
    buffer[y+1] = buffer[x];

    y = kernel_casinh(y, 1);

    buffer[res] = buffer[y+1];
    buffer[res+1] = buffer[y];
  }

  return res;
}

/* Return arc tangent of complex float type. */
// complex catan(complex x)
function catan(x) {
  var res = index;
  index += 2;
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (rcls == FP_INFINITE) {
        buffer[res] = M_COPYSIGN(M_PI_2, buffer[x]);
        buffer[res+1] = M_COPYSIGN(0, buffer[x+1]);
      } else if (icls == FP_INFINITE) {
        if (rcls >= FP_ZERO)
          buffer[res] = M_COPYSIGN(M_PI_2, buffer[x]);
        else
          buffer[res] = M_NAN;
        buffer[res+1] = M_COPYSIGN(0, buffer[x+1]);
      } else if (icls == FP_ZERO || icls == FP_INFINITE) {
        buffer[res] = M_NAN;
        buffer[res+1] = M_COPYSIGN(0, buffer[x+1]);
      } else {
        buffer[res] = M_NAN;
        buffer[res+1] = M_NAN;
      }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    res = x;
  } else {
    if (M_FABS(buffer[x]) >= 16 / M_EPSILON
      || M_FABS(buffer[x+1]) >= 16 / M_EPSILON) {
      buffer[res] = M_COPYSIGN(M_PI_2, buffer[x]);
      if (M_FABS(buffer[x]) <= 1)
        buffer[res+1] = 1 / buffer[x+1];
      else if (M_FABS(buffer[x+1]) <= 1)
        buffer[res+1] = buffer[x+1] / buffer[x] / buffer[x];
      else {
          var h = M_HYPOT(buffer[x] / 2, buffer[x+1] / 2);
          buffer[res+1] = buffer[x+1] / h / h / 4;
        }
    } else {
      var den, absx, absy;

      absx = M_FABS(buffer[x]);
      absy = M_FABS(buffer[x+1]);
      if (absx < absy) {
          var t = absx;
          absx = absy;
          absy = t;
        }

      if (absy < M_EPSILON / 2) {
          den = (1 - absx) * (1 + absx);
          if (den == 0)
            den = 0;
        } else if (absx >= 1)
        den = (1 - absx) * (1 + absx) - absy * absy;
      else if (absx >= 0.75 || absy >= 0.5)
        den = -M_X2Y2M1(absx, absy);
      else
        den = (1 - absx) * (1 + absx) - absy * absy;

      buffer[res] = 0.5 * M_ATAN2(2 * buffer[x], den);

      if (M_FABS(buffer[x+1]) == 1
          && M_FABS(buffer[x]) < M_EPSILON * M_EPSILON)
        buffer[res+1] = (M_COPYSIGN(0.5, buffer[x+1])
                      * (M_LN2 - M_LOG(M_FABS(buffer[x]))));
      else {
        var r2 = 0, num, f;

        if (M_FABS(buffer[x]) >= M_EPSILON * M_EPSILON)
          r2 = buffer[x] * buffer[x];

        num = buffer[x+1] + 1;
        num = r2 + num * num;

        den = buffer[x+1] - 1;
        den = r2 + den * den;

        f = num / den;
        if (f < 0.5)
          buffer[res+1] = 0.25 * M_LOG(f);
        else {
          num = 4 * buffer[x+1];
          buffer[res+1] = 0.25 * M_LOG1P(num / den);
        }
      }
    }

    math_check_force_underflow_complex(res);
  }

  return res;
}

/* Return arc hyperbolic sine for a complex float type. */
// complex casinh(complex x)
function casinh(x) {
  var res = index;
  index += 2;
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      buffer[res] = M_COPYSIGN(M_HUGE_VAL, buffer[x]);

      if (rcls == FP_NAN)
        buffer[res+1] = M_NAN;
      else
        buffer[res+1] = M_COPYSIGN((rcls >= FP_ZERO
                                    ? M_PI_2 : M_PI_4),
                                    buffer[x+1]);
    } else if (rcls <= FP_INFINITE) {
      buffer[res] = buffer[x];
      if ((rcls == FP_INFINITE && icls >= FP_ZERO)
          || (rcls == FP_NAN && icls == FP_ZERO))
        buffer[res+1] = M_COPYSIGN(0, buffer[x+1]);
      else
        buffer[res+1] = M_NAN;
    } else {
      buffer[res] = M_NAN;
      buffer[res+1] = M_NAN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    res = x;
  } else {
    res = kernel_casinh(x, 0);
  }

  return res;
}

/* Return arc hyperbolic cosine for a complex type. */
// complex cacosh(complex x)
function cacosh(x) {
  var res = index;
  index += 2;
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      buffer[res] = M_HUGE_VAL;

      if (rcls == FP_NAN)
        buffer[res+1] = M_NAN;
      else
        buffer[res+1] = M_COPYSIGN((rcls == FP_INFINITE
                                    ? (buffer[x] < 0
                                       ? M_PI - M_PI_4
                                       : M_PI_4)
                                    : M_PI_2), buffer[x+1]);
    } else if (rcls == FP_INFINITE) {
      buffer[res] = M_HUGE_VAL;

      if (icls >= FP_ZERO)
        buffer[res+1] = M_COPYSIGN(signbit(buffer[x])
                                   ? M_PI : 0, buffer[x+1]);
      else
        buffer[res+1] = M_NAN;
    } else {
      buffer[res] = M_NAN;
      if (rcls == FP_ZERO)
        buffer[res+1] = M_PI_2;
      else
        buffer[res+1] = M_NAN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    buffer[res] = 0;
    buffer[res+1] = M_COPYSIGN(M_PI_2, buffer[x+1]);
  } else {
    var y = index;
    index += 2;

    buffer[y] = -buffer[x+1];
    buffer[y+1] = buffer[x];

    y = kernel_casinh(y, 1);

    if (signbit(buffer[x+1])) {
      buffer[res] = buffer[y];
      buffer[res+1] = -buffer[y+1];
    } else {
      buffer[res] = -buffer[y];
      buffer[res+1] = buffer[y+1];
    }
  }

  return res;
}

/* Return arc hyperbolic tangent for a complex float type. */
// complex catanh(complex x)
function catanh(x) {
  var res = index;
  index += 2;
  var rcls = fpclassify(buffer[x]);
  var icls = fpclassify(buffer[x+1]);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      buffer[res] = M_COPYSIGN(0, buffer[x]);
      buffer[res+1] = M_COPYSIGN(M_PI_2, buffer[x+1]);
    } else if (rcls == FP_INFINITE || rcls == FP_ZERO) {
      buffer[res] = M_COPYSIGN(0, buffer[x]);
      if (icls >= FP_ZERO)
        buffer[res+1] = M_COPYSIGN(M_PI_2, buffer[x+1]);
      else
        buffer[res+1] = M_NAN;
    } else {
      buffer[res] = M_NAN;
      buffer[res+1] = M_NAN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    res = x;
  } else {
    if (M_FABS(buffer[x]) >= 16 / M_EPSILON
        || M_FABS(buffer[x+1]) >= 16 / M_EPSILON) {
      buffer[res+1] = M_COPYSIGN(M_PI_2, buffer[x+1]);
      if (M_FABS(buffer[x+1]) <= 1)
        buffer[res] = 1 / buffer[x];
      else if (M_FABS(buffer[x]) <= 1)
        buffer[res] = buffer[x] / buffer[x+1] / buffer[x+1];
      else {
        var h = M_HYPOT(buffer[x] / 2, buffer[x+1] / 2);
        buffer[res] = buffer[x] / h / h / 4;
      }
    } else {
      if (M_FABS(buffer[x]) == 1
          && M_FABS(buffer[x+1]) < M_EPSILON * M_EPSILON)
        buffer[res] = (M_COPYSIGN(0.5, buffer[x])
                      * (M_LN2 - M_LOG(M_FABS(buffer[x+1]))));
      else {
        var i2 = 0;
        if (M_FABS(buffer[x+1]) >= M_EPSILON * M_EPSILON)
          i2 = buffer[x+1] * buffer[x+1];

        var num = 1 + buffer[x];
        num = i2 + num * num;

        var den = 1 - buffer[x];
        den = i2 + den * den;

        var f = num / den;
        if (f < 0.5)
          buffer[res] = 0.25 * M_LOG(f);
        else {
          num = 4 * buffer[x];
          buffer[res] = 0.25 * M_LOG1P(num / den);
        }
      }

      var absx, absy, den;

      absx = M_FABS(buffer[x]);
      absy = M_FABS(buffer[x+1]);
      if (absx < absy) {
        var t = absx;
        absx = absy;
        absy = t;
      }

      if (absy < M_EPSILON / 2) {
        den = (1 - absx) * (1 + absx);
        if (den == 0)
          den = 0;
      } else if (absx >= 1)
      den = (1 - absx) * (1 + absx) - absy * absy;
      else if (absx >= 0.75 || absy >= 0.5)
        den = -M_X2Y2M1(absx, absy);
      else
        den = (1 - absx) * (1 + absx) - absy * absy;

      buffer[res+1] = 0.5 * M_ATAN2(2 * buffer[x+1], den);
    }

    math_check_force_underflow_complex(res);
  }

  return res;
}

// END GLIBC PORT //////////////////////////////////////////////////////////////

// BEGIN MISCELLANEOUS FUNCTION IMPLEMENTATIONS ////////////////////////////////
// * NOTE: THE FUNCTIONS WITHIN THIS SECION ARE NOT A PART OF GLIBC. THEY HAVE
// * NOT BEEN OPTIMIZED NOR HAVE THEY BEEN TESTED OVER ALL EDGE CASES.

// TYPE CONVERSION
// * NOTE: THESE WERE IMPORTED FROM AN EARLIER LIBRARY I WROTE AND THUS ARE
// * SLIGHTLY INCONSISTENT WITH THE STYLISTIC CONVENTIONS USED THROUGHOUT THE
// * REST OF THE FILE.

function fromScalar(n) {
  buffer[index++] = n;
  buffer[index++] = 0;
  return index - 2;
}

function fromPolar(r,phi) {
  if (Number.isNaN(r) || !Number.isFinite(phi)) {
    buffer[index++] = NaN;
    buffer[index++] = NaN;
    return index - 2;
  }
  if (!Number.isFinite(r)) {
    var sin = Math.sin(phi);
    var cos = Math.cos(phi);
    if (Math.abs(sin) < Number.EPSILON) {
      buffer[index++] = r*cos;
      buffer[index++] = 0;
      return index - 2;
    }
    if (Math.abs(cos) < Number.EPSILON) {
      buffer[index++] = 0;
      buffer[index++] = r*sin;
      return index - 2;
    }
    buffer[index++] = r*cos;
    buffer[index++] = r*sin;
    return index - 2;
  }
}

function fromString(str) {
  str = str.replaceAll(/\s/g,'');
  str = str.replaceAll('*',''); // support 'Infinity*i'
  if (!str.length) {
    buffer[index++] = 0;
    buffer[index++] = 0;
    return index - 2;
  }
  if (str.includes('NaN')) {
    buffer[index++] = NaN;
    buffer[index++] = NaN;
    return index - 2;
  }
  var imagRegex = /^([+-]?(?:\d*\.\d+|\d+|Infinity)?i)$/;
  if (imagRegex.test(str)) {
    if (str === 'i') {
      buffer[index++] = 0;
      buffer[index++] = 1;
      return index - 2;
    }
    if (str === '-i') {
      buffer[index++] = 0;
      buffer[index++] = -1;
      return index - 2;
    }
    var im = +str.slice(0,-1);
    if (Number.isNaN(im)) {
      buffer[index++] = NaN;
      buffer[index++] = NaN;
      return index - 2;
    }
    buffer[index++] = 0;
    buffer[index++] = im;
    return index - 2;
  }
  var regex =
    /^([+-]?(?:\d*\.\d+|\d+|Infinity))?([+-]?(?:\d*\.\d+|\d+|Infinity)?i)?$/;
  var match = regex.exec(str);
  var re = +match[1];
  if (Number.isNaN(re)) {
    buffer[index++] = NaN;
    buffer[index++] = NaN;
    return index - 2;
  }
  var im = match[2] ?? '0';
  if (im === '+i' || im === 'i') im = 1;
  else if (im === '-i') im = -1;
  else im = +im.slice(0,-1);
  if (Number.isNaN(im)) {
    buffer[index++] = NaN;
    buffer[index++] = NaN;
    return index - 2;
  }
  buffer[index++] = re;
  buffer[index++] = im;
  return index - 2;
}

function fromArray(array,offset=0) {
  buffer[index++] = array[offset];
  buffer[index++] = array[offset+1];
  return index - 2;
}

function fromArrayPolar(array,offset=0) {
  return fromPolar(
    array[offset],
    array[offset+1],
  );
}

// slow, keep away from hot loops
function from() {
  switch (arguments.length) {
    case 0:
      buffer[index++] = NaN;
      buffer[index++] = NaN;
      return index - 2;
    case 1:
      var arg = arguments[0];
      if (Array.isArray(arg)) {
        switch (arg.length) {
          case 0:
            buffer[index++] = NaN;
            buffer[index++] = NaN;
            return index - 2;
          case 1:
            buffer[index++] = arg[0];
            buffer[index++] = 0;
            return index - 2;
          default:
            buffer[index++] = arg[0];
            buffer[index++] = arg[1];
            return index - 2;
        }
      }
      if (typeof arg === 'function') arg = arg();
      switch (typeof arg) {
        case 'number':
          buffer[index++] = arg;
          buffer[index++] = 0;
          return index - 2;
        case 'bigint':
          buffer[index++] = Number(arg);
          buffer[index++] = 0;
          return index - 2;
        case 'string': return fromString(arg);
        case 'object':
          var keys = Object.getOwnPropertyNames(arg);
          var re = arg.re ?? arg.real ?? arg.x ?? arg[2] ?? null; 
          var im = arg.im ?? arg.imag ?? arg.y ?? arg[1] ?? null;
          var r = arg.r ?? arg.rho ?? null;
          var phi = arg.phi ?? arg.theta ?? null;
          if (re === null) {
            if (im === null) {
              if (r === null) {
                if (phi === null) {
                  buffer[index++] = NaN;
                  buffer[index++] = NaN;
                  return index - 2;
                }
                return fromPolar(1,phi);
              }
              if (phi === null) {
                buffer[index++] = r;
                buffer[index++] = 0;
                return index - 2;
              }
              return fromPolar(r,phi);
            }
            if (r === null || phi === null) {
              buffer[index++] = 0;
              buffer[index++] = im;
              return index - 2;
            }
            return fromPolar(r,phi);
          }
          if (im === null) {
            if (r === null || phi === null) {
              buffer[index++] = re;
              buffer[index++] = 0;
              return index - 2;
            }
            return fromPolar(r,phi);
          }
          buffer[index++] = re;
          buffer[index++] = im;
          return index - 2;
        default:
          arg = +arg;
          if (Number.isNaN(arg)) {
            buffer[index++] = NaN;
            buffer[index++] = NaN;
            return index - 2;
          }
          buffer[index++] = arg;
          buffer[index++] = 0;
          return index - 2;
      }
    default:
      var re = +arguments[0];
      var im = +arguments[1];
      if (Number.isNaN(re) || Number.isNaN(im)) {
        buffer[index++] = NaN;
        buffer[index++] = NaN;
        return index - 2;
      }
      buffer[index++] = re;
      buffer[index++] = im;
      return index - 2;
  }
}

function toArray(z,array=[],offset=0) {
  array[offset] = buffer[z];
  array[offset+1] = buffer[z+1];
  return array;
}

function toArrayPolar(z,array=[],offset=0) {
  var re = buffer[z];
  var im = buffer[z+1];
  array[offset] = M_HYPOT(re,im);
  array[offset+1] = Math.atan2(im,re);
  return array;
}

// BOOLEAN OPERATIONS

/* Complex equality test */
// boolean ceq(compelx z, complex w, float t)
function ceq(z,w,t=M_EPSILON) {
  if (z === w) return true;
  var t2 = t*t;
  return M_FABS(buffer[z] - buffer[w]) < t
      && M_FABS(buffer[z+1] - buffer[w+1]) < t;
}

/* Complex NaN test */
// boolean cisNaN(compelx z)
function cisNaN(z) {
  var rcls = fpclassify(buffer[z]);
  var icls = fpclassify(buffer[z+1]);
  return rcls === FP_NAN || icls === FP_NAN;
}

/* Complex finite test */
// boolean cisFinite(compelx z)
function cisFinite(z) {
  var rcls = fpclassify(buffer[z]);
  var icls = fpclassify(buffer[z+1]);
  return rcls > FP_INFINITE && icls > FP_INFINITE;
}

/* Complex zero test */
// boolean cisZero(compelx z)
function cisZero(z) {
  var rcls = fpclassify(buffer[z]);
  var icls = fpclassify(buffer[z+1]);
  return rcls === FP_ZERO && icls === FP_ZERO;
}

/* Complex real test */
// boolean cisReal(compelx z)
function cisReal(z) {
  var icls = fpclassify(buffer[z+1]);
  return icls === FP_ZERO;
}

// COMPONENT OPERATIONS

/* Compute floor of complex components */
// complex floor(complex z)
function cfloor(z) {
  var retval = index;
  index += 2;
  buffer[retval] = Math.floor(buffer[z]);
  buffer[retval+1] = Math.floor(buffer[z+1]);
  return retval;
}

/* Compute ceiling of complex components */
// complex ceil(complex z)
function cceil(z) {
  var retval = index;
  index += 2;
  buffer[retval] = Math.ceil(buffer[z]);
  buffer[retval+1] = Math.ceil(buffer[z+1]);
  return retval;
}

/* Compute rounded complex components */
// complex round(complex z)
function cround(z) {
  var retval = index;
  index += 2;
  buffer[retval] = Math.round(buffer[z]);
  buffer[retval+1] = Math.round(buffer[z+1]);
  return retval;
}

/* Compute complex components rounded to the nearest 32-bit float */
// complex fround(complex z)
function cfround(z) {
  var retval = index;
  index += 2;
  buffer[retval] = Math.fround(buffer[z]);
  buffer[retval+1] = Math.fround(buffer[z+1]);
  return retval;
}

/* Compute complex components truncated towards zero */
// complex fround(complex z)
function ctrunc(z) {
  var retval = index;
  index += 2;
  buffer[retval] = Math.trunc(buffer[z]);
  buffer[retval+1] = Math.trunc(buffer[z+1]);
  return retval;
}

/* Compute random complex components */
// complex crandom()
function crandom() {
  var retval = index;
  index += 2;
  buffer[retval] = Math.random();
  buffer[retval+1] = Math.random();
  return retval;
}

/* Compute complex component minimum */
// complex cmin(...complex z)
function cmin(...z) {
  var x = Infinity;
  var y = Infinity;
  var l = z.length;
  var i = 0;
  for (;i < l;) {
    if (buffer[z[i]] < x) x = buffer[z[i]];
    if (buffer[z[i]+1] < y) y = buffer[z[i]+1];
    i++;
  }
  return c(x, y);
}

/* Compute complex component maximum */
// complex cmax(...complex z)
function cmax(...z) {
  var x = -Infinity;
  var y = -Infinity;
  var l = z.length;
  var i = 0;
  for (;i < l;) {
    if (buffer[z[i]] > x) x = buffer[z[i]];
    if (buffer[z[i]+1] > y) y = buffer[z[i]+1];
    i++;
  }
  return c(x, y);
}

// RECIPROCAL TRIG FUNCTIONS

/* Compute complex secant */
// sec(z) = 1/cos(z)
// complex csec(complex z)
function csec(z) {
  return cinv(ccos(z));
}

/* Compute complex cosecant */
// csc(z) = 1/sin(z)
// complex ccsc(complex z)
function ccsc(z) {
  return cinv(csin(z));
}

/* Compute complex cotangent */
// cot(z) = 1/tan(z)
// complex ccot(complex z)
function ccot(z) {
  return cinv(ctan(z));
}

/* Compute complex hyperbolic secant */
// sech(z) = 1/cosh(z)
// complex csech(complex z)
function csech(z) {
  return cinv(ccosh(z));
}

/* Compute complex hyperbolic cosecant */
// csch(z) = 1/sinh(z)
// complex ccsch(complex z)
function ccsch(z) {
  return cinv(csinh(z));
}

/* Compute complex hyperbolic cotangent */
// coth(z) = 1/tanh(z)
// complex ccoth(complex z)
function ccoth(z) {
  return cinv(ctanh(z));
}

/* Compute complex arc secant */
// asec(z) = acos(1/z)
// complex casec(complex z)
function casec(z) {
  return cacos(cinv(z));
}

/* Compute complex arc cosecant */
// acsc(z) = asin(1/z)
// complex cacsc(complex z)
function cacsc(z) {
  return casin(cinv(z));
}

/* Compute complex arc cotangent */
// acot(z) = atan(1/z)
// complex cacot(complex z)
function cacot(z) {
  return catan(cinv(z));
}

/* Compute complex hyperbolic arc secant */
// asech(z) = acosh(1/z)
// complex casech(complex z)
function casech(z) {
  return cacosh(cinv(z));
}

/* Compute complex hyperbolic arc cosecant */
// acsch(z) = asinh(1/z)
// complex cacsch(complex z)
function cacsch(z) {
  return casinh(cinv(z));
}

/* Compute complex hyperbolic arc cotangent */
// acoth(z) = atanh(1/z)
// complex cacoth(complex z)
function cacoth(z) {
  return catanh(cinv(z));
}

// OTHER MISCELLANEOUS FUNCTIONS

/* Compute complex square */
// complex csquare(complex z)
function csquare(z) {
  var retval = index;
  index += 2;
  var x = buffer[z];
  var y = buffer[z+1];
  buffer[retval] = x*x - y*y;
  buffer[retval+1] = 2*x*y;
  return retval;
}

/* Compute complex cube */
// complex ccube(complex z)
function ccube(z) {
  var retval = index;
  index += 2;
  var x = buffer[z];
  var y = buffer[z+1];
  var x2 = x*x;
  var y2 = y*y;
  buffer[retval] = x*x2 - 3*x*y2;
  buffer[retval+1] = 3*y*x2 - y*y2;
  return retval;
}

/* Compute complex sum */
// complex csum(...complex z)
function csum(...z) {
  return z.reduce(cadd);
}

/* Compute complex product */
// complex cprod(...complex z)
function cprod(...z) {
  return z.reduce(cmul);
}

/* Compute complex unit given angle */
// cis(x) = cos(x) + i*sin(x)
// complex cis(float x)
function cis(x) {
  var retval = index;
  index += 2;
  buffer[retval] = Math.cos(x);
  buffer[retval+1] = Math.sin(x);
  return retval;
}

/* Compute complex sine cardinal */
// sinc(0) = 1
// sinc(z) = sin(z)/z, z != 0
// complex csinc(complex z)
function csinc(z) {
  if (buffer[z] === 0 && buffer[z+1] === 0) {
    var retval = index;
    index += 2;
    buffer[retval] = 1;
    buffer[retval+1] = 0;
    return retval;
  }
  return cdiv(csin(z),z);
}

/* Compute complex mod */
// mod(z,w) = z - w*floor(z/w)
// complex cmod(complex z, complex w)
function cmod(z, w) {
  var flrzw = cdiv(z,w);
  buffer[flrzw] = Math.floor(buffer[flrzw]);
  buffer[flrzw+1] = Math.floor(buffer[flrzw+1]);
  return csub(z,cmul(w,flrzw));
}

/* Compute complex log2 */
// lg(z) = log(z) / log(2)
// complex clg(complex z)
function clg(z) {
  return cdivScalar(clog(z),M_LN2);
}

/* Compute complex gamma */
// complex cgamma(complex z)
function cgamma(z) {
  if (buffer[z] < 0.5) {
    return gamma_left(z);
  } else {
    return gamma_right(z);
  }
}

function gamma_left(z) {
  // reflection formula
  // Γ(z) = π / (sin(πz) * Γ(1-z))
  var sinpiz = csin(cmulScalar(z,M_PI));
  var gamma1z = gamma_right(caddScalar(cneg(z),1));
  var dinv = cinv(cmul(sinpiz, gamma1z));
  return cmulScalar(dinv,M_PI);
}

function gamma_right(z) {
  // Lanczos approximation
  var w = caddScalar(z, -1);
  var t = caddScalar(w, 7.5);
  var x = c(0.99999999999980993, 0);
  x = cadd(x, cmulScalar(cinv(caddScalar(w, 1)), 676.5203681218851));
  x = csub(x, cmulScalar(cinv(caddScalar(w, 2)), 1259.1392167224028));
  x = cadd(x, cmulScalar(cinv(caddScalar(w, 3)), 771.32342877765313));
  x = csub(x, cmulScalar(cinv(caddScalar(w, 4)), 176.61502916214059));
  x = cadd(x, cmulScalar(cinv(caddScalar(w, 5)), 12.507343278686905));
  x = csub(x, cmulScalar(cinv(caddScalar(w, 6)), .13857109526572012));
  x = cadd(x, cmulScalar(cinv(caddScalar(w, 7)), 9.9843695780195716e-6));
  x = cadd(x, cmulScalar(cinv(caddScalar(w, 8)), 1.5056327351493116e-7));
  // sqrt(2π) * x * t^(w+0.5) * e^(-t)
  var xroottau = cmulScalar(x,M_SQRT2PI);
  var twp1_2 = cpow(t,caddScalar(w,0.5));
  var expnt = cexp(cneg(t));
  return cmul(cmul(xroottau,twp1_2),expnt);
}

/* Compute complex factorial */
// z! = Γ(z + 1)
// complex cfact(complex z)
function cfact(z) {
  return cgamma(caddScalar(z,1));
}

/* Compute complex beta */
// B(a,b) = Γ(a) Γ(b) / Γ(a+b)
// complex cbeta(complex z, complex w)
function cbeta(z, w) {
  return cdiv(cmul(cgamma(z),cgamma(w)),cgamma(cadd(z,w)));
}

/* Compute complex binomial coefficient */
// C(z,w) = 1 / ((z+1) * B(w+1, z-w+1))
// complex cbinom(complex z, complex w)
function cbinom(z, w) {
  var zmw = csub(z,w);
  return cdiv(z, cmul(cmul(w, zmw), cbeta(zmw, w)));
}

/* Compute complex derivative */
// f'(z) = (f(z+h) - f(z)) / h
// (comlex -> complex) cprime((complex -> complex) f, float h)
function cprime(f, h=1e-7) {
  return z => {
    var zh = caddScalar(z, h);
    var fzh = f(zh);
    var fz = f(z);
    return cmulScalar(csub(fzh,fz), 1/h);
  };
}

// END MISCELLANEOUS FUNCTION IMPLEMENTATIONS //////////////////////////////////

var Complex = {resetBuffer: resetBuffer, setBufferSize: setBufferSize, toString: ctoString, create: c, neg: cneg, sgn: csgn, inv: cinv, add: cadd, sub: csub, mul: cmul, div: cdiv, addScalar: caddScalar, subScalar: csubScalar, mulScalar: cmulScalar, divScalar: cdivScalar, mulImaginary: cmuli, real: creal, imag: cimag, abs: cabs, arg: carg, conj: conj, proj: cproj, exp: cexp, log: clog, log10: clog10, pow: cpow, sqrt: csqrt, sin: csin, cos: ccos, tan: ctan, sinh: csinh, cosh: ccosh, tanh: ctanh, asin: casin, acos: cacos, atan: catan, asinh: casinh, acosh: cacosh, atanh: catanh, sec: csec, csc: ccsc, cot: ccot, sech: csech, csch: ccsch, coth: ccoth, asec: casec, acsc: cacsc, acot: cacot, asech: casech, acsch: cacsch, acoth: cacoth, sinc: csinc, cis: cis, square: csquare, cube: ccube, floor: cfloor, ceil: cceil, round: cround, fround: cfround, trunc: ctrunc, random: crandom, min: cmin, max: cmax, mod: cmod, log2: clg, gamma: cgamma, fact: cfact, beta: cbeta, binom: cbinom, derivative: cprime, equal: ceq, isNaN: cisNaN, isFinite: cisFinite, isZero: cisZero, isReal: cisReal, polar: fromPolar, fromString: fromString, fromArray: fromArray, fromArrayPolar: fromArrayPolar, from: from, toArray: toArray, toArrayPolar: toArrayPolar};

Object.defineProperties(Complex, {
  BUFFER: { get() { return buffer; } },
  BUFFER_SIZE: {
    get() { return bufferSize; },
    set: setBufferSize,
  },
  INDEX: {
    get() { return index; },
    set(v) { index = v; },
  },
  STRIDE: {
    get() { return 2; },
  }
});

if (typeof module !== 'undefined') module.exports = Complex;
else globalThis.Complex = Complex;

})();