(() => {// ──── BEGIN IIFE ─────────────────────────────────────────────────────

// ──── CONSTANTS ──────────────────────────────────────────────────────────────
var M_E             = 2.71828182845904523536028747135266; // e
var M_1_E           = 0.36787944117144232159552377016146; // 1/e
var M_LNPI          = 1.14472988584940017414342735135306; // ln(pi)
var M_53_LN2        = 36.7368005696771013991133024372834; // 53*ln(2)
var M_LOG2          = 0.30102999566398119521373889472449; // log10(2)
var M_53_LOG2       = 15.9545897701910033463281614203981; // 53*log10(2)

var M_LN10          = 2.30258509299404568401799145468436; // ln(10)
var M_1_LN10        = 0.43429448190325182765112891891661; // 1/ln(10)
var M_1_2_LN10      = 0.21714724095162591382556445945830; // 1/2/ln(10)
var M_PI_LN10       = 1.36437635384184134748578362543136; // pi/ln(10)
var M_3PI4_LN10     = 1.02328226538138101061433771907352; // 3*pi/4/ln(10)
var M_PI_2_LN10     = 0.68218817692092067374289181271568; // pi/2/ln(10)
var M_PI_4_LN10     = 0.34109408846046033687144590635784; // pi/4/ln(10)

var M_LN2           = 0.69314718055994530941723212145818; // ln(2)
var M_1_LN2         = 1.44269504088896340735992468100189; // 1/ln(2)
var M_1_2_LN2       = 0.72134752044448170367996234050095; // 1/2/ln(2)
var M_PI_LN2        = 4.53236014182719380962768294571667; // pi/ln(2)
var M_3PI4_LN2      = 3.39927010637039535722076220928750; // 3*pi/4/ln(2)
var M_PI_2_LN2      = 2.26618007091359690481384147285833; // pi/2/ln(2)
var M_PI_4_LN2      = 1.13309003545679845240692073642917; // pi/4/ln(2)

var M_2PI           = 6.28318530717958647692528676655901; // 2*pi
var M_PI            = 3.14159265358979323846264338327950; // pi
var M_3PI_2         = 4.71238898038468985769396507491925; // 3*pi/2
var M_3PI_4         = 2.35619449019234492884698253745963; // 3*pi/4
var M_PI_2          = 1.57079632679489661923132169163975; // pi/2
var M_PI_4          = 0.78539816339744830961566084581988; // pi/4
var M_1_PI          = 0.31830988618379067153776752674503; // 1/pi
var M_2_PI          = 0.63661977236758134307553505349006; // 2/pi
var M_1_SQRTPI      = 0.56418958354775628694807945156077; // 1/sqrt(pi)
var M_2_SQRTPI      = 0.63661977236758134307553505349006; // 2/sqrt(pi)
var M_SQRT2PI       = 2.50662827463100050241576528481105; // sqrt(2*pi)
var M_SQRT2         = 1.41421356237309504880168872420970; // sqrt(2)
var M_SQRT1_2       = 0.70710678118654752440084436210485; // 1/sqrt(2)

var M_MANT_DIG      = 53;
var M_MIN_EXP       = -1021;                              // 2**-1021 underflows
var M_MAX_EXP       = 1024;                               // 2**1024 overflows
var M_MIN           = 2.225073858507201383090232717e-308; // 2**-1022
var M_2_MIN         = 4.450147717014402766180465435e-308; // 2**-1021
var M_4_MIN         = 8.900295434028805532360930869e-308; // 2**-1020
var M_EPSILON_SQ    = 4.9303806576313237838233035330e-32; // 2**-104
var M_EPSILON_8     = 2.7755575615628913510590791702e-17; // 2**-55
var M_EPSILON_2     = 1.1102230246251565404236316681e-16; // 2**-53
var M_EPSILON       = 2.2204460492503130808472633362e-16; // 2**-52
var M_SQRTEPSI      = 1.49011611938476562500000000000e-8; // sqrt(2**-52)
var M_2_27P1        = 1.34217729000000000000000000000e+8; // 2**27+1
var M_1_EPSILON     = 4.5035996273704960000000000000e+15; // 2**52
var M_2_EPSILON     = 9.0071992547409920000000000000e+15; // 2**53
var M_16_EPSILON    = 7.2057594037927936000000000000e+16; // 2**56
var M_MAX_4         = 4.494232837155789769323262977e+307; // 2**1022
var M_MAX_2         = 8.988465674311579538646525954e+307; // 2**1023
var M_MAX           = 1.797693134862315907729305191e+308; // 2**1024
var M_EXP709        = 8.218407461554972189241372387e+307; // e**709

var FP_NAN          = 0;
var FP_INFINITE     = 1;
var FP_ZERO         = 2;
var FP_SUBNORMAL    = 3;
var FP_NORMAL       = 4;


// ──── UTILITY ────────────────────────────────────────────────────────────────
function fpclassify(x) {
  if (x !== x) return FP_NAN;
  if (x === Infinity) return FP_INFINITE;
  if (x === -Infinity) return FP_INFINITE;
  if (x === 0) return FP_ZERO;
  if (x > 0 && x < M_MIN) return FP_SUBNORMAL;
  if (x < 0 && x > -M_MIN) return FP_SUBNORMAL;
  return FP_NORMAL;
}

function signbit(x) {
  if (x === 0) return 1/x > 0 ? 0 : 1;
  return x > 0 ? 0 : 1;
}

function isinf(x) {
  return x === Infinity || x === -Infinity;
}

function isfinite(x) {
  return x !== Infinity && x !== -Infinity;
}

var hypot = (() => {
  var SCALE     = 2.409919865102884e-181;
  var LARGE_VAL = 6.703903964971299e+153;
  var TINY_VAL  = 6.717876107567089e-139;
  var EPS       = 5.551115123125783e-17;

  function hypot(x,y) {
    x = fabs(x);
    y = fabs(y);
    if (x < y) {var t = x; x = y; y = t;}
    if (x === Infinity) return Infinity;
    if (x !== x) return NaN;
    if (y <= x * EPS) return x + y;
    if (x > LARGE_VAL) {
      x *= SCALE;
      y *= SCALE;
      var h = Math.sqrt(x*x + y*y);
      var t1, t2;
      if (h <= 2 * y) {
        var delta = h - y;
        t1 = x * (2*delta - x);
        t2 = (delta - 2*(x - y)) * delta;
      } else {
        var delta = h - x;
        t1 = 2*delta * (x - 2*y);
        t2 = (4*delta - y)*y + delta*delta;
      }
      return (h - (t1 + t2) / (2 * h)) / SCALE;
    } else if (y < TINY_VAL) {
      x /= SCALE;
      y /= SCALE;
      var h = Math.sqrt(x*x + y*y);
      var t1, t2;
      if (h <= 2 * y) {
        var delta = h - y;
        t1 = x * (2*delta - x);
        t2 = (delta - 2*(x - y)) * delta;
      } else {
        var delta = h - x;
        t1 = 2*delta * (x - 2*y);
        t2 = (4*delta - y)*y + delta*delta;
      }
      return (h - (t1 + t2) / (2 * h)) * SCALE;
    } else {
      var h = Math.sqrt(x*x + y*y);
      var t1, t2;
      if (h <= 2 * y) {
        var delta = h - y;
        t1 = x * (2*delta - x);
        t2 = (delta - 2*(x - y)) * delta;
      } else {
        var delta = h - x;
        t1 = 2*delta * (x - 2*y);
        t2 = (4*delta - y)*y + delta*delta;
      }
      return h - (t1 + t2) / (2 * h);
    }
  };

  return hypot;
})();

function copysign(x, y) {
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

function fabs(x) { // non-zero input or for magnitude comparisons
  if (x > 0) return x;
  return -x;
}

function abs(x) {
  if (x > 0) return x;
  if (x < 0) return -x;
  return 0;
}

var x2y2m1 = (() => {
  var C = 134217729;
  return function(x,y) {
    var p = C * x;
    var xh = p - (p - x), xl = x - xh;
    var xxHi = x * x;
    var xxLo = (xh * xh - xxHi) + 2 * xh * xl + xl * xl;
    p = C * y;
    var yh = p - (p - y), yl = y - yh;
    var yyHi = y * y;
    var yyLo = (yh * yh - yyHi) + 2 * yh * yl + yl * yl;

    var v0 = xxLo, v1 = xxHi, v2 = yyLo, v3 = yyHi, v4 = -1.0;
    var t, h, l;
    var a0, a1, a2, a3, a4;

    // insertion sort by abs ascending
    a0=v0>0?v0:-v0; a1=v1>0?v1:-v1; a2=v2>0?v2:-v2; a3=v3>0?v3:-v3; a4=v4>0?v4:-v4;
    if(a0>a1){t=v0;v0=v1;v1=t;t=a0;a0=a1;a1=t;}
    if(a1>a2){t=v1;v1=v2;v2=t;t=a1;a1=a2;a2=t;
      if(a0>a1){t=v0;v0=v1;v1=t;t=a0;a0=a1;a1=t;}}
    if(a2>a3){t=v2;v2=v3;v3=t;t=a2;a2=a3;a3=t;
      if(a1>a2){t=v1;v1=v2;v2=t;t=a1;a1=a2;a2=t;
        if(a0>a1){t=v0;v0=v1;v1=t;t=a0;a0=a1;a1=t;}}}
    if(a3>a4){t=v3;v3=v4;v4=t;t=a3;a3=a4;a4=t;
      if(a2>a3){t=v2;v2=v3;v3=t;t=a2;a2=a3;a3=t;
        if(a1>a2){t=v1;v1=v2;v2=t;t=a1;a1=a2;a2=t;
          if(a0>a1){t=v0;v0=v1;v1=t;t=a0;a0=a1;a1=t;}}}}

    h=v1+v0; l=(v1-h)+v0; v0=l; v1=h;
    a1=v1>0?v1:-v1; a2=v2>0?v2:-v2; a3=v3>0?v3:-v3; a4=v4>0?v4:-v4;
    if(a1>a2){t=v1;v1=v2;v2=t;t=a1;a1=a2;a2=t;}
    if(a2>a3){t=v2;v2=v3;v3=t;t=a2;a2=a3;a3=t;
      if(a1>a2){t=v1;v1=v2;v2=t;t=a1;a1=a2;a2=t;}}
    if(a3>a4){t=v3;v3=v4;v4=t;t=a3;a3=a4;a4=t;
      if(a2>a3){t=v2;v2=v3;v3=t;t=a2;a2=a3;a3=t;
        if(a1>a2){t=v1;v1=v2;v2=t;t=a1;a1=a2;a2=t;}}}

    h=v2+v1; l=(v2-h)+v1; v1=l; v2=h;
    a2=v2>0?v2:-v2; a3=v3>0?v3:-v3; a4=v4>0?v4:-v4;
    if(a2>a3){t=v2;v2=v3;v3=t;t=a2;a2=a3;a3=t;}
    if(a3>a4){t=v3;v3=v4;v4=t;t=a3;a3=a4;a4=t;
      if(a2>a3){t=v2;v2=v3;v3=t;t=a2;a2=a3;a3=t;}}

    h=v3+v2; l=(v3-h)+v2; v2=l; v3=h;
    a3=v3>0?v3:-v3; a4=v4>0?v4:-v4;
    if(a3>a4){t=v3;v3=v4;v4=t;}

    h=v4+v3; l=(v4-h)+v3; v3=l; v4=h;
    return v4+v3+v2+v1+v0;
  }
})();

var scalbn = (() => {
  var POW2 = (new Float64Array(2098)).map((_,i) => Math.pow(2, i - 1074));
  function scalbn(x,n) {
    // if (x === 0 || !Number.isFinite(x)) return x;
    // if (n > 1023) return x*Infinity;
    // if (n < -1074) return x*0;
    return x*POW2[n + 1074];
  }
  return scalbn;
})();

/* Return arc hyperbolic sine for a complex float type, with the
   imaginary part of the result possibly adjusted for use in
   computing other functions. */
function kernel_casinh(x, adj) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  var y = new Complex(0,0);

  var rx = xre > 0 ? xre : -xre;
  var ix = xim > 0 ? xim : -xim;

  if (rx >= M_1_EPSILON || ix >= M_1_EPSILON) {
    y.re = rx;
    y.im = ix;

    if (adj) {
      var t = y.re;
      y.re = copysign(y.im, xim);
      y.im = t;
    }

    var res = y.log();
    re = res.re + M_LN2;
    im = res.im;
  } else if (rx >= 0.5 && ix < M_EPSILON_8) {
    var s = hypot(1, rx);

    re = Math.log(rx + s);
    if (adj)
      im = Math.atan2(s, xim);
    else
      im = Math.atan2(ix, s);
  } else if (rx < M_EPSILON_8 && ix >= 1.5) {
    var s = Math.sqrt((ix + 1) * (ix - 1));

    re = Math.log(ix + s);
    if (adj)
      im = Math.atan2(rx, copysign(s, xim));
    else
      im = Math.atan2(s, rx);
  } else if (ix > 1 && ix < 1.5 && rx < 0.5) {
    if (rx < M_EPSILON_SQ) {
      var ix2m1 = (ix + 1) * (ix - 1);
      var s = Math.sqrt(ix2m1);

      re = Math.log1p(2 * (ix2m1 + ix * s)) / 2;
      if (adj)
        im = Math.atan2(rx, copysign(s, xim));
      else
        im = Math.atan2(s, rx);
    } else {
      var ix2m1 = (ix + 1) * (ix - 1);
      var rx2 = rx * rx;
      var f = rx2 * (2 + rx2 + 2 * ix * ix);
      var d = Math.sqrt(ix2m1 * ix2m1 + f);
      var dp = d + ix2m1;
      var dm = f / dp;
      var r1 = Math.sqrt((dm + rx2) / 2);
      var r2 = rx * ix / r1;

      re = Math.log1p(rx2 + dp + 2 * (rx * r1 + ix * r2)) / 2;
      if (adj)
        im = Math.atan2(rx + r1, copysign(ix + r2, xim));
      else
        im = Math.atan2(ix + r2, rx + r1);
    }
  } else if (ix == 1 && rx < 0.5) {
    if (rx < M_EPSILON_8) {
      re = Math.log1p(2 * (rx + Math.sqrt(rx))) / 2;
      if (adj)
        im = Math.atan2(Math.sqrt (rx), copysign(1, xim));
      else
        im = Math.atan2(1, Math.sqrt (rx));
    } else {
      var d = rx * Math.sqrt(4 + rx * rx);
      var s1 = Math.sqrt((d + rx * rx) / 2);
      var s2 = Math.sqrt((d - rx * rx) / 2);

      re = Math.log1p(rx * rx + d + 2 * (rx * s1 + s2)) / 2;
      if (adj)
        im = Math.atan2(rx + s1, copysign(1 + s2, xim));
      else
        im = Math.atan2(1 + s2, rx + s1);
    }
  } else if (ix < 1 && rx < 0.5) {
    if (ix >= M_EPSILON) {
      if (rx < M_EPSILON_SQ) {
        var onemix2 = (1 + ix) * (1 - ix);
        var s = Math.sqrt(onemix2);

        re = Math.log1p(2 * rx / s) / 2;
        if (adj)
          im = Math.atan2(s, xim);
        else
          im = Math.atan2(ix, s);
      }
      else {
        var onemix2 = (1 + ix) * (1 - ix);
        var rx2 = rx * rx;
        var f = rx2 * (2 + rx2 + 2 * ix * ix);
        var d = Math.sqrt(onemix2 * onemix2 + f);
        var dp = d + onemix2;
        var dm = f / dp;
        var r1 = Math.sqrt((dp + rx2) / 2);
        var r2 = rx * ix / r1;

        re = Math.log1p(rx2 + dm + 2 * (rx * r1 + ix * r2)) / 2;
        if (adj)
          im = Math.atan2(rx + r1, copysign(ix + r2, xim));
        else
          im = Math.atan2(ix + r2, rx + r1);
      }
    } else {
      var s = hypot(1, rx);

      re = Math.log1p(2 * rx * (rx + s)) / 2;
      if (adj)
        im = Math.atan2(s, xim);
      else
        im = Math.atan2(ix, s);
    }
  } else {
    y.re = (rx - ix) * (rx + ix) + 1;
    y.im = 2 * rx * ix;

    y.sqrtEq();

    y.re += rx;
    y.im += ix;

    if (adj) {
      var t = y.re;
      y.re = copysign(y.im, xim);
      y.im = t;
    }

    var res = y.log();
    re = res.re;
    im = res.im;
  }

  /* Give results the correct sign for the original argument.  */
  re = copysign(re, xre);
  im = copysign(im, (adj ? 1 : xim));

  return new Complex(re,im);
}


// ──── CONSTRUCTORS ───────────────────────────────────────────────────────────
// create a new complex number given euclidean components
// constructor Complex(x,y)
// static create(x=0,y=0)
// static fromScalar(n=0)
// static fromArray(array,offset=0)

function Complex(x,y) {
  this.re = +x;
  this.im = +y;
}

Complex.create = function(x=0, y=0) {
  return new Complex(x,y);
}

Complex.fromScalar = function(n=0) {
  return new Complex(n,0);
}

Complex.fromArray = function(array, offset=0) {
  return new Complex(array[offset], array[offset+1]);
}

// create a new complex number given polar components
// static polar(r=0,phi=0)
// static cis(phi=0)
// static fromArrayPolar(array,offset=0)

Complex.polar = function(r=0, phi=0) {
  if (phi === Infinity || phi !== phi) return new Complex(NaN, NaN);
  if (phi === -Infinity) return new Complex(NaN, -NaN);
  if (phi === 0) return new Complex(r, phi);
  return new Complex(
    r * Math.cos(phi),
    r * Math.sin(phi),
  );
}

Complex.cis = function(phi=0) {
  if (phi === Infinity || phi !== phi) return new Complex(NaN, NaN);
  if (phi === -Infinity) return new Complex(NaN, -NaN);
  if (phi === 0) return new Complex(1, phi);
  return new Complex(
    Math.cos(phi),
    Math.sin(phi),
  );
}

Complex.fromArrayPolar = function(array, offset=0) {
  var r = array[offset];
  var phi = array[offset+1];
  if (phi === Infinity || phi !== phi) return new Complex(NaN, NaN);
  if (phi === -Infinity) return new Complex(NaN, -NaN);
  if (phi === 0) return new Complex(r, phi);
  return new Complex(
    r * Math.cos(phi),
    r * Math.sin(phi),
  );
}

// ──── COPY / CLONE ───────────────────────────────────────────────────────────
// create a clone of a complex number with the same components
// static clone(z)
// clone()

Complex.clone = function(z) {
  return new Complex(z.re, z.im);
}

Complex.prototype.clone = function() {
  return new Complex(this.re, this.im);
}

// copy the components of a complex number to another
// static clone(z)    : new Complex = z
// static copy(z,w)   : z = w
// static copyTo(z,w) : w = z
// clone()            : new Complex = z
// copy(z)            : this = z
// copyTo(z)          : z = this

Complex.copy = function(z,w) {
  z.re = w.re;
  z.im = w.im;
  return z;
}

Complex.prototype.copy = function(z) {
  this.re = z.re;
  this.im = z.im;
  return this;
}

Complex.copyTo = function(z,w) {
  w.re = z.re;
  w.im = z.im;
  return w;
}

Complex.prototype.copyTo = function(z) {
  z.re = this.re;
  z.im = this.im;
  return z;
}


// ──── TYPE CONVERSION ────────────────────────────────────────────────────────
// write a complex numbers components into a new or existing array
// static toArray(z,array=[],offset=0)
// static toArrayPolar(z,array=[],offset=0)
// toArray(array=[],offset=0)
// toArrayPolar(array=[],offset=0)

Complex.toArray = function(z, array=[], offset=0) {
  array[offset] = z.re;
  array[offset+1] = z.im;
  return array;
}

Complex.toArrayPolar = function(z, array=[], offset=0) {
  var x = z.re;
  var y = z.im;
  array[offset] = hypot(x,y);
  array[offset+1] = Math.atan2(y,x);
  return array;
}

Complex.prototype.toArray = function(z, array=[], offset=0) {
  array[offset] = this.re;
  array[offset+1] = this.im;
  return array;
}

Complex.prototype.toArrayPolar = function(z, array=[], offset=0) {
  var x = this.re;
  var y = this.im;
  array[offset] = hypot(x,y);
  array[offset+1] = Math.atan2(y,x);
  return array;
}

// get a string representation of a complex number
// toString()

Complex.prototype.toString = function() {
  if (arguments.length) {
    var x = this.re.toString(...arguments);
    var y = this.im.toString(...arguments);
    return `${x} + ${y}*i`;
  }
  var x = this.re;
  var y = this.im;
  var str;
  if (x === 0) {
    if (1/x > 0) {
      if (y !== y || y === Infinity || y === -Infinity) return `${y}*i`;
      else if (y === 0) return 1/y > 0 ? "0" : "-0i";
      else if (y === 1) return "i";
      else if (y === -1) return "-i";
      else return `${y}i`;
    } else str = "-0"
  } else str = `${x}`;
  if (y !== y || y === Infinity) return str + ` + ${y}*i`;
  else if (y === -Infinity) return str + " - Infinity*i";
  else if (y === 0 && 1/y > 0) return str;
  else if (y === 1) return str + " + i";
  else if (y === -1) return str + " - i";
  else if (y <= 0) return str + ` - ${-y}i`;
  else return str + ` + ${y}i`;
}

// ──── SETTERS ────────────────────────────────────────────────────────────────
// set the euclidean real/imaginary components of a complex number
// static setReal(z,x)
// static setImag(z,y)
// static set(z,x,y)
// static setScalar(z,n)
// static setComponent(z,i,n)
// static setFromArray(z,array,offset=0)
// setReal(x)
// setImag(y)
// set(x,y)
// setScalar(n)
// setComponent(i,n)
// setFromArray(array,offset=0)

Complex.setReal = function(z, x) {
  z.re = x;
  return z;
}

Complex.setImag = function(z, y) {
  z.im = y;
  return z;
}

Complex.set = function(z, x, y) {
  z.re = x;
  z.im = y;
  return z;
}

Complex.setScalar = function(z, n) {
  z.re = n;
  z.im = 0;
  return z;
}

Complex.setComponent = function(z, i, n) {
  if (i === 0 || i === -2) z.re = n;
  else if (i === 1 || i === -1) z.im = n;
  return z;
}

Complex.setfromArray = function(z, array, offset=0) {
  z.re = array[offset];
  z.im = array[offset+1] ?? 0;
  return z;
}

Complex.prototype.setReal = function(x) {
  this.re = x;
  return z;
}

Complex.prototype.setImag = function(y) {
  this.im = y;
  return z;
}

Complex.prototype.set = function(x, y) {
  this.re = x;
  this.im = y;
  return z;
}

Complex.prototype.setScalar = function(n) {
  this.re = n;
  this.im = 0;
  return z;
}

Complex.prototype.setComponent = function(i, n) {
  if (i === 0 || i === -2) this.re = n;
  else if (i === 1 || i === -1) this.im = n;
  return z;
}

Complex.prototype.setfromArray = function(array, offset=0) {
  this.re = array[offset];
  this.im = array[offset+1] ?? 0;
  return z;
}

// set the polar magnitude/argument components of a complex number
// static setAbs(z,r)
// static setArg(z,phi)
// static setPolar(z,x,y)
// static setFromArrayPolar(z,array,offset=0)
// setAbs(r)
// setArg(phi)
// setPolar(r,phi)
// setFromArrayPolar(array,offset=0)

Complex.setAbs = function(z, r) {/*
  var x = z.re;
  var y = z.im;
  if (x === Infinity) {
    if (y === Infinity) {
      z.re = r * M_SQRT1_2;
      z.im = r * M_SQRT1_2;
    } else if (y === -Infinity) {
      z.re = r * M_SQRT1_2;
      z.im = -r * M_SQRT1_2;
    } else if (y === 0) {
      z.re = r;
    } else {
      z.re = r;
      z.im *= 0;
    }
  } else if (x === -Infinity) {
    if (y === Infinity) {
      z.re = -r * M_SQRT1_2;
      z.im = -r * M_SQRT1_2;
    } else if (y === -Infinity) {
      z.re = -r * M_SQRT1_2;
      z.im = r * M_SQRT1_2;
    } else if (y === 0) {
      z.re = -r;
    } else {
      z.re = -r;
      z.im *= 0;
    }
  }
  r /= hypot(z.re, z.im);
  z.re *= r;
  z.im *= r;
  return z;*/
  var phi = Math.atan2(z.im, z.re);
  z.re = r * Math.cos(phi);
  z.im = r * Math.sin(phi);
  return z;
}

Complex.setArg = function(z, phi) {
  if (phi === Infinity || phi !== phi) {
    z.re = NaN;
    z.im = NaN;
  } else if (phi === -Infinity) {
    z.re = NaN;
    z.im = -NaN;
  } else if (phi === 0) {
    var r = hypot(z.re, z.im);
    z.re = r;
    z.im = phi;
  } else {
    var r = hypot(z.re, z.im);
    z.re = r * Math.cos(phi);
    z.im = r * Math.sin(phi);
  }
  return z;
}

Complex.setPolar = function(z, r, phi) {
  if (phi === Infinity || phi !== phi) {
    z.re = NaN;
    z.im = NaN;
  } else if (phi === -Infinity) {
    z.re = NaN;
    z.im = -NaN;
  } else if (phi === 0) {
    z.re = r;
    z.im = phi;
  } else {
    z.re = r * Math.cos(phi);
    z.im = r * Math.sin(phi);
  }
  return z;
}

Complex.setFromArrayPolar = function(z, array, offset=0) {
  var r = array[offset];
  var phi = array[offset+1] ?? 0;
  if (phi === Infinity || phi !== phi) {
    z.re = NaN;
    z.im = NaN;
  } else if (phi === -Infinity) {
    z.re = NaN;
    z.im = -NaN;
  } else if (phi === 0) {
    z.re = r;
    z.im = phi;
  } else {
    z.re = r * Math.cos(phi);
    z.im = r * Math.sin(phi);
  }
  return z;
}

Complex.prototype.setAbs = function(r) {
  var phi = Math.atan2(this.im, this.re);
  this.re = r * Math.cos(phi);
  this.im = r * Math.sin(phi);
  return this;
}

Complex.prototype.setArg = function(phi) {
  if (phi === Infinity || phi !== phi) {
    this.re = NaN;
    this.im = NaN;
  } else if (phi === -Infinity) {
    this.re = NaN;
    this.im = -NaN;
  } else if (phi === 0) {
    var r = hypot(this.re, this.im);
    this.re = r;
    this.im = phi;
  } else {
    var r = hypot(this.re, this.im);
    this.re = r * Math.cos(phi);
    this.im = r * Math.sin(phi);
  }
  return this;
}

Complex.prototype.setPolar = function(r, phi) {
  if (phi === Infinity || phi !== phi) {
    this.re = NaN;
    this.im = NaN;
  } else if (phi === -Infinity) {
    this.re = NaN;
    this.im = -NaN;
  } else if (phi === 0) {
    this.re = r;
    this.im = phi;
  } else {
    this.re = r * Math.cos(phi);
    this.im = r * Math.sin(phi);
  }
  return this;
}

Complex.prototype.setFromArrayPolar = function(array, offset=0) {
  var r = array[offset];
  var phi = array[offset+1] ?? 0;
  if (phi === Infinity || phi !== phi) {
    this.re = NaN;
    this.im = NaN;
  } else if (phi === -Infinity) {
    this.re = NaN;
    this.im = -NaN;
  } else if (phi === 0) {
    this.re = r;
    this.im = phi;
  } else {
    this.re = r * Math.cos(phi);
    this.im = r * Math.sin(phi);
  }
  return this;
}

// ──── GETTERS ────────────────────────────────────────────────────────────────
// get the euclidean real/imaginary components of a complex number
// static getReal(z)
// static getImage(z)
// static getComponent(z,i)
// static complexReal(z)
// static complexImag(z)
// getReal()
// getImage()
// getComponent(i)
// getComplexReal()
// getComplexImag()

Complex.getReal = function(z) {
  return z.re;
}

Complex.getImag = function(z) {
  return z.im;
}

Complex.getComponent = function(z, i) {
  if (i === 0 || i === -2) return z.re;
  else if (i === 1 || i === -1) return z.im;
  return z;
}

Complex.getComplexReal = function(z) {
  return new Complex(z.re, 0);
}

Complex.getComplexImag = function(z) {
  return new Complex(z.im, 0);
}

Complex.prototype.getReal = function() {
  return this.re;
}

Complex.prototype.getImag = function() {
  return this.im;
}

Complex.prototype.getComponent = function(i) {
  if (i === 0 || i === -2) return this.re;
  else if (i === 1 || i === -1) return this.im;
}

Complex.prototype.getComplexReal = function() {
  return new Complex(this.re, 0);
}

Complex.prototype.getComplexImag = function() {
  return new Complex(this.im, 0);
}

// get the polar magnitude/argument components of a complex number
// static getAbs(z)
// static getArg(z)
// static getComplexAbs(z)
// static getComplexArg(z)
// getAbs()
// getArg()
// getComplexAbs()
// getComplexArg()

Complex.getAbs = function(z) {
  return hypot(z.re, z.im);
}

Complex.getArg = function(z) {
  return Math.atan2(z.im, z.re);
}

Complex.getComplexAbs = function(z) {
  return new Complex(hypot(z.re, z.im), 0);
}

Complex.getComplexArg = function(z) {
  return new Complex(Math.atan2(z.im, z.re), 0);
}

Complex.prototype.getAbs = function() {
  return hypot(this.re, this.im);
}

Complex.prototype.getArg = function() {
  return Math.atan2(this.im, this.re);
}

Complex.prototype.getComplexAbs = function() {
  return new Complex(hypot(this.re, this.im), 0);
}

Complex.prototype.getComplexArg = function() {
  return new Complex(Math.atan2(this.im, this.re), 0);
}


// ──── PSEUDO PROPERTIES ──────────────────────────────────────────────────────
// get/set the components of a complex number via pseudo-properties
// get real()
// set real(x)
// get imag()
// set imag(y)
// get r()
// set r(r)
// get arg()
// set arg(phi)

Object.defineProperties(Complex.prototype, {
  real: {
    get: Complex.prototype.getReal,
    set: Complex.prototype.setReal,
    enumerable: false,
    configurable: true,
  },
  imag: {
    get: Complex.prototype.getImag,
    set: Complex.prototype.setImag,
    enumerable: false,
    configurable: true,
  },
  r: {
    get: Complex.prototype.getAbs,
    set: Complex.prototype.setAbs,
    enumerable: false,
    configurable: true,
  },
  theta: {
    get: Complex.prototype.getArg,
    set: Complex.prototype.setArg,
    enumerable: false,
    configurable: true,
  },
});

Complex.real = Complex.getReal;
Complex.imag = Complex.getImag;
Complex.abs = Complex.getAbs;
Complex.arg = Complex.getArg;

// ──── ADDITIVE OPERATIONS ────────────────────────────────────────────────────
// add complex numbers
// static add(z,w)
// add(z)
// addEq(z)
//
// out of 1_000_000 random inputs:
// off by 0ULP: 100%
// off by 1ULP+: 0%
//
// average error: +/- 0 ULP

Complex.add = function(z,w) {
  return new Complex(
    z.re + w.re,
    z.im + w.im,
  );
}

Complex.prototype.add = function(z) {
  return new Complex(
    this.re + z.re,
    this.im + z.im,
  );
}

Complex.prototype.addEq = function(z) {
  this.re += z.re;
  this.im += z.im;
  return this;
}

// add a scalar to a complex number
// static addScalar(z,n)
// addScalar(n)
// addScalarEq(n)
//
// out of 1_000_000 random inputs:
// off by 0ULP: 100%
// off by 1ULP+: 0%
//
// average error: +/- 0 ULP

Complex.addScalar = function(z,n) {
  return new Complex(
    z.re + n,
    z.im,
  );
}

Complex.prototype.addScalar = function(n) {
  return new Complex(
    this.re + n,
    this.im,
  );
}

Complex.prototype.addScalarEq = function(n) {
  this.re += n;
  return this;
}

// subtract complex numbers
// static sub(z,w)
// sub(z)
// subEq(z)
//
// out of 1_000_000 random inputs:
// off by 0ULP: 100%
// off by 1ULP+: 0%
//
// average error: +/- 0 ULP

Complex.sub = function(z,w) {
  return new Complex(
    z.re - w.re,
    z.im - w.im,
  );
}

Complex.prototype.sub = function(z) {
  return new Complex(
    this.re - z.re,
    this.im - z.im,
  );
}

Complex.prototype.subEq = function(z) {
  this.re -= z.re;
  this.im -= z.im;
  return this;
}

// subtract a scalar from a complex number
// static subScalar(z,n)
// subScalar(n)
// subScalarEq(n)
//
// out of 1_000_000 random inputs:
// off by 0ULP: 100%
// off by 1ULP+: 0%
//
// average error: +/- 0 ULP

Complex.subScalar = function(z,n) {
  return new Complex(
    z.re - n,
    z.im,
  );
}

Complex.prototype.subScalar = function(n) {
  return new Complex(
    this.re - n,
    this.im,
  );
}

Complex.prototype.subScalarEq = function(n) {
  this.re -= n;
  return this;
}

// subtract a complex number from a scalar
// static scalarSub(n,z)
// scalarSub(n)
// scalarSubEq(n)
//
// out of 1_000_000 random inputs:
// off by 0ULP: 100%
// off by 1ULP+: 0%
//
// average error: +/- 0 ULP

Complex.scalarSub = function(n,z) {
  return new Complex(
    n - z.re,
    -z.im,
  );
}

Complex.prototype.scalarSub = function(n) {
  return new Complex(
    n - this.re,
    -this.im,
  );
}

Complex.prototype.scalarSubEq = function(n) {
  this.re = n - this.re;
  this.im *= -1;
  return this;
}

// negate a complex number
// static neg(z)
// neg()
// negEq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 100%
// off by 1ULP+: 0%
//
// average error: +/- 0 ULP

Complex.neg = function(z) {
  return new Complex(
    -z.re,
    -z.im,
  );
}

Complex.prototype.neg = function() {
  return new Complex(
    -this.re,
    -this.im,
  );
}

Complex.prototype.negEq = function() {
  this.re *= -1;
  this.im *= -1;
  return this;
}

// take the conjugate of a complex number
// static conj(z)
// conj()
// conjEq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 100%
// off by 1ULP+: 0%
//
// average error: +/- 0 ULP

Complex.conj = function(z) {
  return new Complex(
    z.re,
    -z.im,
  );
}

Complex.prototype.conj = function() {
  return new Complex(
    this.re,
    -this.im,
  );
}

Complex.prototype.conjEq = function() {
  this.im *= -1;
  return this;
}

// ──── MULTIPLICATIVE OPERATIONS ──────────────────────────────────────────────
// multiply complex numbers
// static mul(z,w)
// mul(z)
// mulEq(z)
//
// out of 1_000_000 random inputs:
// off by 0ULP: 99.20215%
// off by 1ULP: 0.79755%
// off by 2ULP: 0.00015%
// off by 3ULP: 0%
// off by 4ULP: 0.00005%
// off by 5ULP+: 0%
//
// average error: +/- 0.0080005 ULP

Complex.mul = function(z,w) {
  var a = z.re;
  var b = z.im;
  var c = w.re;
  var d = w.im;
  var ac = a * c;
  var bd = b * d;
  var ad = a * d;
  var bc = b * c;
  var x = ac - bd;
  var y = ad + bc;

  if (x !== x || y !== y || x === 0 || y === 0) {
    if (x !== x) {
      x = Infinity * (a/b - d/c);
      if (c < 0) x *= -1;
      if (b < 0) x *= -1;
    }
    if (y !== y) {
      y = Infinity * (a/b + c/d);
      if (d < 0) y *= -1;
      if (b < 0) y *= -1;
    }
    if (x === 0) {
      x = 0 * (a/b - d/c);
      if (c < 0) x *= -1;
      if (b < 0) x *= -1;
    }
    if (y === 0) {
      y = 0 * (a/b + c/d);
      if (d < 0) y *= -1;
      if (b < 0) y *= -1;
    }
  } else if (Math.abs(x) < (Math.abs(ac) + Math.abs(bd)) ||
             Math.abs(y) < (Math.abs(ad) + Math.abs(bc))) {
    var ta = a * M_2_27P1;
    var tb = b * M_2_27P1;
    var tc = c * M_2_27P1;
    var td = d * M_2_27P1;
    var ahi = ta - (ta - a);
    var bhi = tb - (tb - b);
    var chi = tc - (tc - c);
    var dhi = td - (td - d);
    var alo = a - ahi;
    var blo = b - bhi;
    var clo = c - chi;
    var dlo = d - dhi;
    var ace = ahi*chi - ac + ahi*clo + alo*chi + alo*clo;
    var bde = bhi*dhi - bd + bhi*dlo + blo*dhi + blo*dlo;
    var ade = ahi*dhi - ad + ahi*dlo + alo*dhi + alo*dlo;
    var bce = bhi*chi - bc + bhi*clo + blo*chi + blo*clo;
    var xe = ace - bde;
    var ye = ade + bce;
    if (xe === xe) x += ace - bde;
    if (ye === ye) y += ade + bce;
  }

  return new Complex(x,y);
}

Complex.prototype.mul = function(z) {
  var a = this.re;
  var b = this.im;
  var c = z.re;
  var d = z.im;
  var ac = a * c;
  var bd = b * d;
  var ad = a * d;
  var bc = b * c;
  var x = ac - bd;
  var y = ad + bc;

  if (x !== x || y !== y || x === 0 || y === 0) {
    if (x !== x) {
      x = Infinity * (a/b - d/c);
      if (c < 0) x *= -1;
      if (b < 0) x *= -1;
    }
    if (y !== y) {
      y = Infinity * (a/b + c/d);
      if (d < 0) y *= -1;
      if (b < 0) y *= -1;
    }
    if (x === 0) {
      x = 0 * (a/b - d/c);
      if (c < 0) x *= -1;
      if (b < 0) x *= -1;
    }
    if (y === 0) {
      y = 0 * (a/b + c/d);
      if (d < 0) y *= -1;
      if (b < 0) y *= -1;
    }
  } else if (Math.abs(x) < (Math.abs(ac) + Math.abs(bd)) ||
             Math.abs(y) < (Math.abs(ad) + Math.abs(bc))) {
    var ta = a * M_2_27P1;
    var tb = b * M_2_27P1;
    var tc = c * M_2_27P1;
    var td = d * M_2_27P1;
    var ahi = ta - (ta - a);
    var bhi = tb - (tb - b);
    var chi = tc - (tc - c);
    var dhi = td - (td - d);
    var alo = a - ahi;
    var blo = b - bhi;
    var clo = c - chi;
    var dlo = d - dhi;
    var ace = ahi*chi - ac + ahi*clo + alo*chi + alo*clo;
    var bde = bhi*dhi - bd + bhi*dlo + blo*dhi + blo*dlo;
    var ade = ahi*dhi - ad + ahi*dlo + alo*dhi + alo*dlo;
    var bce = bhi*chi - bc + bhi*clo + blo*chi + blo*clo;
    var xe = ace - bde;
    var ye = ade + bce;
    if (xe === xe) x += ace - bde;
    if (ye === ye) y += ade + bce;
  }

  return new Complex(x,y);
}

Complex.prototype.mulEq = function(z) {
  var a = this.re;
  var b = this.im;
  var c = z.re;
  var d = z.im;
  var ac = a * c;
  var bd = b * d;
  var ad = a * d;
  var bc = b * c;
  var x = ac - bd;
  var y = ad + bc;

  if (x !== x || y !== y || x === 0 || y === 0) {
    if (x !== x) {
      x = Infinity * (a/b - d/c);
      if (c < 0) x *= -1;
      if (b < 0) x *= -1;
    }
    if (y !== y) {
      y = Infinity * (a/b + c/d);
      if (d < 0) y *= -1;
      if (b < 0) y *= -1;
    }
    if (x === 0) {
      x = 0 * (a/b - d/c);
      if (c < 0) x *= -1;
      if (b < 0) x *= -1;
    }
    if (y === 0) {
      y = 0 * (a/b + c/d);
      if (d < 0) y *= -1;
      if (b < 0) y *= -1;
    }
  } else if (Math.abs(x) < (Math.abs(ac) + Math.abs(bd)) ||
             Math.abs(y) < (Math.abs(ad) + Math.abs(bc))) {
    var ta = a * M_2_27P1;
    var tb = b * M_2_27P1;
    var tc = c * M_2_27P1;
    var td = d * M_2_27P1;
    var ahi = ta - (ta - a);
    var bhi = tb - (tb - b);
    var chi = tc - (tc - c);
    var dhi = td - (td - d);
    var alo = a - ahi;
    var blo = b - bhi;
    var clo = c - chi;
    var dlo = d - dhi;
    var ace = ahi*chi - ac + ahi*clo + alo*chi + alo*clo;
    var bde = bhi*dhi - bd + bhi*dlo + blo*dhi + blo*dlo;
    var ade = ahi*dhi - ad + ahi*dlo + alo*dhi + alo*dlo;
    var bce = bhi*chi - bc + bhi*clo + blo*chi + blo*clo;
    var xe = ace - bde;
    var ye = ade + bce;
    if (xe === xe) x += ace - bde;
    if (ye === ye) y += ade + bce;
  }

  this.re = x;
  this.im = y;
  return this;
}

// multiply a complex number by a scalar
// static mulScalar(z,n)
// mulScalar(n)
// mulScalarEq(n)
//
// out of 1_000_000 random inputs:
// off by 0ULP: 99.9876%
// off by 1ULP: 0.0124%
// off by 2ULP+: 0%
//
// average error: +/- 0.000124 ULP

Complex.mulScalar = function(z,n) {
  return new Complex(
    z.re * n,
    z.im * n,
  );
}

Complex.prototype.mulScalar = function(n) {
  return new Complex(
    this.re * n,
    this.im * n,
  );
}

Complex.prototype.mulScalarEq = function(n) {
  this.re *= n;
  this.im *= n;
  return this;
}

// divide complex numbers
// static div(z,w)
// div(z)
// divEq(z)

Complex.div = function(z,w) {
  var a = z.re;
  var b = z.im;
  var c = w.re;
  var d = w.im;
  var denom;
  var ratio;
  var x;
  var y;

  if (fabs(c) < fabs(d)) {
    ratio = c / d;
    denom = (c * ratio) + d;
    x = ((a * ratio) + b) / denom;
    y = ((b * ratio) - a) / denom;
  } else {
    ratio = d / c;
    denom = (d * ratio) + c;
    x = ((b * ratio) + a) / denom;
    y = (b - (a * ratio)) / denom;
  }

  if (x !== x) {
    x = Infinity * (a/b - d/c);
    if (c < 0 !== b < 0) x *= -1;
  }
  if (y !== y) {
    y = Infinity * (a/b + c/d);
    if (d < 0 !== b < 0) y *= -1;
  }
  if (x === 0) {
    x = 0 * (a/b - d/c);
    if (c < 0 !== b < 0) x *= -1;
  }
  if (y === 0) {
    y = 0 * (a/b + c/d);
    if (d < 0 !== b < 0) y *= -1;
  }

  return new Complex(x,y);
}

Complex.prototype.div = function(z) {
  var a = this.re;
  var b = this.im;
  var c = z.re;
  var d = z.im;
  var denom;
  var ratio;
  var x;
  var y;

  if (fabs(c) < fabs(d)) {
    ratio = c / d;
    denom = (c * ratio) + d;
    x = ((a * ratio) + b) / denom;
    y = ((b * ratio) - a) / denom;
  } else {
    ratio = d / c;
    denom = (d * ratio) + c;
    x = ((b * ratio) + a) / denom;
    y = (b - (a * ratio)) / denom;
  }

  if (x !== x) {
    x = Infinity * (a/b - d/c);
    if (c < 0 !== b < 0) x *= -1;
  }
  if (y !== y) {
    y = Infinity * (a/b + c/d);
    if (d < 0 !== b < 0) y *= -1;
  }
  if (x === 0) {
    x = 0 * (a/b - d/c);
    if (c < 0 !== b < 0) x *= -1;
  }
  if (y === 0) {
    y = 0 * (a/b + c/d);
    if (d < 0 !== b < 0) y *= -1;
  }

  return new Complex(x,y);
}

Complex.prototype.divEq = function(z) {
  var a = this.re;
  var b = this.im;
  var c = z.re;
  var d = z.im;
  var denom;
  var ratio;
  var x;
  var y;

  if (fabs(c) < fabs(d)) {
    ratio = c / d;
    denom = (c * ratio) + d;
    x = ((a * ratio) + b) / denom;
    y = ((b * ratio) - a) / denom;
  } else {
    ratio = d / c;
    denom = (d * ratio) + c;
    x = ((b * ratio) + a) / denom;
    y = (b - (a * ratio)) / denom;
  }

  if (x !== x) {
    x = Infinity * (a/b - d/c);
    if (c < 0 !== b < 0) x *= -1;
  }
  if (y !== y) {
    y = Infinity * (a/b + c/d);
    if (d < 0 !== b < 0) y *= -1;
  }
  if (x === 0) {
    x = 0 * (a/b - d/c);
    if (c < 0 !== b < 0) x *= -1;
  }
  if (y === 0) {
    y = 0 * (a/b + c/d);
    if (d < 0 !== b < 0) y *= -1;
  }

  this.re = x;
  this.im = y
  return this;
}

// divide a complex number by a scalar
// static divScalar(z,n)
// divScalar(n)
// divScalarEq(n)

Complex.divScalar = function(z,n) {
  return new Complex(
    z.re / n,
    z.im / n,
  );
}

Complex.prototype.divScalar = function(n) {
  return new Complex(
    this.re * n,
    this.im * n,
  );
}

Complex.prototype.divScalarEq = function(n) {
  this.re /= n;
  this.im /= n;
  return this;
}

// divide a complex number into a scalar
// static scalarDiv(z,n)
// scalarDiv(n)
// scalarDivEq(n)

Complex.scalarDiv = function(n,z) {
  var c = z.re;
  var d = z.im;
  var denom;
  var ratio;
  var x;
  var y;

  if (fabs(c) < fabs(d)) {
    ratio = c / d;
    denom = (c * ratio) + d;
    x = n*ratio / denom;
    y = -n / denom;
  } else {
    ratio = d / c;
    denom = (d * ratio) + c;
    x = n / denom;
    y = -n*ratio / denom;
  }

  if (x !== x && y !== y) {
    if (denom === 0) {
      x = copysign(Infinity, c) * n;
      y = copysign(Infinity, c) * 0;
    } else if (isinf(n) && isfinite(c) && isfinite(d)) {
      n = copysign(1, n);
      x = Infinity * (n * c + 0 * d);
      y = Infinity * (0 * c - n * d);
    } else if ((isinf(c) || isinf(d)) && isfinite(n)) {
      c = copysign(isinf(c) ? 1 : 0, c);
      d = copysign(isinf(d) ? 1 : 0, d);
      x = 0 * (n * c + 0 * d);
      y = 0 * (0 * c - n * d);
    }
  }

  return new Complex(x,y);
}

Complex.prototype.scalarDiv = function(n) {
  var c = this.re;
  var d = this.im;
  var denom;
  var ratio;
  var x;
  var y;

  if (fabs(c) < fabs(d)) {
    ratio = c / d;
    denom = (c * ratio) + d;
    x = n*ratio / denom;
    y = -n / denom;
  } else {
    ratio = d / c;
    denom = (d * ratio) + c;
    x = n / denom;
    y = -n*ratio / denom;
  }

  if (x !== x && y !== y) {
    if (denom === 0) {
      x = copysign(Infinity, c) * n;
      y = copysign(Infinity, c) * 0;
    } else if (isinf(n) && isfinite(c) && isfinite(d)) {
      n = copysign(1, n);
      x = Infinity * (n * c + 0 * d);
      y = Infinity * (0 * c - n * d);
    } else if ((isinf(c) || isinf(d)) && isfinite(n)) {
      c = copysign(isinf(c) ? 1 : 0, c);
      d = copysign(isinf(d) ? 1 : 0, d);
      x = 0 * (n * c + 0 * d);
      y = 0 * (0 * c - n * d);
    }
  }

  return new Complex(x,y);
}

Complex.prototype.scalarDivEq = function(n) {
  var c = this.re;
  var d = this.im;
  var denom;
  var ratio;
  var x;
  var y;

  if (fabs(c) < fabs(d)) {
    ratio = c / d;
    denom = (c * ratio) + d;
    x = n*ratio / denom;
    y = -n / denom;
  } else {
    ratio = d / c;
    denom = (d * ratio) + c;
    x = n / denom;
    y = -n*ratio / denom;
  }

  if (x !== x && y !== y) {
    if (denom === 0) {
      x = copysign(Infinity, c) * n;
      y = copysign(Infinity, c) * 0;
    } else if (isinf(n) && isfinite(c) && isfinite(d)) {
      n = copysign(1, n);
      x = Infinity * (n * c + 0 * d);
      y = Infinity * (0 * c - n * d);
    } else if ((isinf(c) || isinf(d)) && isfinite(n)) {
      c = copysign(isinf(c) ? 1 : 0, c);
      d = copysign(isinf(d) ? 1 : 0, d);
      x = 0 * (n * c + 0 * d);
      y = 0 * (0 * c - n * d);
    }
  }

  this.re = x;
  this.im = y;
  return this;
}

// take the reciprocal of a complex number
// static inv(z)
// inv()
// invEq()

Complex.inv = function(z) {
  var c = z.re;
  var d = z.im;
  var denom;
  var ratio;
  var x;
  var y;

  if (fabs(c) < fabs(d)) {
    ratio = c / d;
    denom = (c * ratio) + d;
    x = ratio / denom;
    y = -1 / denom;
  } else {
    ratio = d / c;
    denom = (d * ratio) + c;
    x = 1 / denom;
    y = -ratio / denom;
  }

  if (x !== x && y !== y) {
    if (denom === 0) {
      x = copysign(Infinity, c) * 1;
      y = copysign(Infinity, c) * 0;
    } else if (isinf(c) || isinf(d)) {
      c = copysign(isinf(c) ? 1 : 0, c);
      d = copysign(isinf(d) ? 1 : 0, d);
      x = 0 * (1 * c + 0 * d);
      y = 0 * (0 * c - 1 * d);
    }
  }

  return new Complex(x,y);
}

Complex.prototype.inv = function() {
  var c = this.re;
  var d = this.im;
  var denom;
  var ratio;
  var x;
  var y;

  if (fabs(c) < fabs(d)) {
    ratio = c / d;
    denom = (c * ratio) + d;
    x = ratio / denom;
    y = -1 / denom;
  } else {
    ratio = d / c;
    denom = (d * ratio) + c;
    x = 1 / denom;
    y = -ratio / denom;
  }

  if (x !== x && y !== y) {
    if (denom === 0) {
      x = copysign(Infinity, c) * 1;
      y = copysign(Infinity, c) * 0;
    } else if (isinf(c) || isinf(d)) {
      c = copysign(isinf(c) ? 1 : 0, c);
      d = copysign(isinf(d) ? 1 : 0, d);
      x = 0 * (1 * c + 0 * d);
      y = 0 * (0 * c - 1 * d);
    }
  }

  return new Complex(x,y);
}

Complex.prototype.invEq = function() {
  var c = this.re;
  var d = this.im;
  var denom;
  var ratio;
  var x;
  var y;

  if (fabs(c) < fabs(d)) {
    ratio = c / d;
    denom = (c * ratio) + d;
    x = ratio / denom;
    y = -1 / denom;
  } else {
    ratio = d / c;
    denom = (d * ratio) + c;
    x = 1 / denom;
    y = -ratio / denom;
  }

  if (x !== x && y !== y) {
    if (denom === 0) {
      x = copysign(Infinity, c) * 1;
      y = copysign(Infinity, c) * 0;
    } else if (isinf(c) || isinf(d)) {
      c = copysign(isinf(c) ? 1 : 0, c);
      d = copysign(isinf(d) ? 1 : 0, d);
      x = 0 * (1 * c + 0 * d);
      y = 0 * (0 * c - 1 * d);
    }
  }

  this.re = x;
  this.im = y;
  return this;
}

// take the directed sign of a complex number
// static sgn(z)
// sgn()
// sgnEq()

Complex.sgn = function(z) {
  var a = z.re;
  var b = z.im;
  var n = hypot(a, b);

  var x;
  var y;

  x = a / n;
  y = b / n;

  if (x !== x && y !== y) {
    if (n === 0 && (a === a || b === b)) {
      x = copysign(0, a);
      y = copysign(0, b);
    } else if ((isinf(a) || isinf(b)) && isfinite(n)) {
      a = copysign(isinf(a) ? 1 : 0, a);
      b = copysign(isinf(b) ? 1 : 0, b);
      x = Infinity * (a * n + b * 0);
      y = Infinity * (b * n - a * 0);
    } else if (isinf(n) && isfinite(a) && isfinite(b)) {
      n = copysign(isinf(n) ? 1 : 0, n);
      x = 0 * (a * n + b * 0);
      y = 0 * (b * n - a * 0);
    }
  }

  return new Complex(x,y);
}

Complex.prototype.sgn = function() {
  var a = this.re;
  var b = this.im;
  var n = hypot(a, b);

  var x;
  var y;

  x = a / n;
  y = b / n;

  if (x !== x && y !== y) {
    if (n === 0 && (a === a || b === b)) {
      x = copysign(0, a);
      y = copysign(0, b);
    } else if ((isinf(a) || isinf(b)) && isfinite(n)) {
      a = copysign(isinf(a) ? 1 : 0, a);
      b = copysign(isinf(b) ? 1 : 0, b);
      x = Infinity * (a * n + b * 0);
      y = Infinity * (b * n - a * 0);
    } else if (isinf(n) && isfinite(a) && isfinite(b)) {
      n = copysign(isinf(n) ? 1 : 0, n);
      x = 0 * (a * n + b * 0);
      y = 0 * (b * n - a * 0);
    }
  }

  return new Complex(x,y);
}

Complex.prototype.sgnEq = function() {
  var a = this.re;
  var b = this.im;
  var n = hypot(a, b);

  var x;
  var y;

  x = a / n;
  y = b / n;

  if (x !== x && y !== y) {
    if (n === 0 && (a === a || b === b)) {
      x = copysign(0, a);
      y = copysign(0, b);
    } else if ((isinf(a) || isinf(b)) && isfinite(n)) {
      a = copysign(isinf(a) ? 1 : 0, a);
      b = copysign(isinf(b) ? 1 : 0, b);
      x = Infinity * (a * n + b * 0);
      y = Infinity * (b * n - a * 0);
    } else if (isinf(n) && isfinite(a) && isfinite(b)) {
      n = copysign(isinf(n) ? 1 : 0, n);
      x = 0 * (a * n + b * 0);
      y = 0 * (b * n - a * 0);
    }
  }

  this.re = x;
  this.im = y;
  return this;
}

// take the absolute value of a complex number
// * defined above
//
// out of 1_000_000 random inputs:
// off by 0ULP: 99.9989%
// off by 1ULP: 0.0011%
// off by 2ULP+: 0%
//
// average error: +/- 0.000011 ULP

// take the argument of a complex number
// *defined above
//
// out of 1_000_000 random inputs:
// off by 0ULP: 99.5471%
// off by 1ULP: 0.4529%
// off by 2ULP+: 0%
//
// average error: +/- 0.004529 ULP

// project a complex number onto the riemann sphere
// static proj(z)
// prog()
// projEq()

Complex.proj = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  if (xre === Infinity || xim === Infinity) {
    return new Complex(Infinity, copysign(0, im));
  } else {
    return new Complex(re, im);
  }
}

Complex.prototype.proj = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  if (xre === Infinity || xim === Infinity) {
    return new Complex(Infinity, copysign(0, im));
  } else {
    return new Complex(re, im);
  }
}

Complex.prototype.projEq = function() {
  var xre = this.re;
  var xim = this.im;
  if (xre === Infinity || xim === Infinity) {
    this.re = Infinity;
    this.im = copysign(0, im);
  }
  return this;
}

// ──── EXPONENTIALS & LOGARITHMS ──────────────────────────────────────────────
// take e to the power of a complex number
// static exp(z)
// exp()
// expEq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 98.3926%
// off by 1ULP: 1.60035%
// off by 2ULP: 0.00705%
// off by 3ULP+: 0%
//
// average error: +/- 0.0161445 ULP

Complex.exp = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  if (xre !== xre) {
    // exp(NaN + 0i) = NaN + 0i
    // exp(NaN - 0i) = NaN - 0i
    // exp(NaN + non-zero i) = NaN + NaN i
    re = NaN;
    im = xim === 0 ? xim : NaN;
  } else if (xre === Infinity) {
    if (isFinite(xim)) {
      if (xim === 0) {
        // exp(+Infinity + 0i) = +Infinity + 0i
        // exp(+Infinity - 0i) = +Infinity - 0i
        re = Infinity;
        im = xim;
      } else {
        if (xim > M_MIN || xim < -M_MIN) {
          // re = +Infinity, |im| > M_MIN = infinities with signs of cos/sin
          re = Infinity/Math.cos(xim);
          im = Infinity/Math.sin(xim);
        } else {
          // same as above but |im| < M_MIN, meaning sin/cos is linear
          re = Infinity;
          im = Infinity/xim;
        }
      }
    } else {
      // re = +Infinity, im = NaN or +/-Infinity = Infinity + NaN*i
      re = Infinity;
      im = NaN;
    }
  } else if (xre === -Infinity) {
    if (isFinite(xim)) {
      if (xim === 0) {
        // exp(-Infinity + 0i) = +0 + 0i
        // exp(-Infinity - 0i) = +0 - 0i
        re = 0;
        im = xim;
      } else {
        if (xim > M_MIN || xim < -M_MIN) {
          // re = -Infinity, |im| > M_MIN = zeros with signs of cos/sin
          re = 0*Math.cos(xim);
          im = 0*Math.sin(xim);
        } else {
          // same as above but |im| < M_MIN, meaning sin/cos is linear
          re = 0;
          im = 0*xim;
        }
      }
    } else {
      // re = -Infinity, im = NaN or +/-Infinity = 0 + NaN*i
      re = 0;
      im = NaN;
    }
  } else {
    if (isFinite(xim)) {
      var sinix;
      var cosix;
      if (xim > M_MIN || xim < -M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        // sin/cos is linear
        sinix = xim;
        cosix = 1;
      }
      if (xre > 709) {
        // very small chance this prevents overflow in the output
        xre -= 709;
        sinix *= M_EXP709;
        cosix *= M_EXP709;
      }
      if (xre > 709) {
        // these could be infinity or zero, * M_MAX preserve the signs
        re = M_MAX * cosix;
        im = M_MAX * sinix;
      } else {
        // otherwise Math.exp(xre) wont overflow,
        var expre = Math.exp(xre);
        // although these may still overflow, but will have the correct signs
        // and wont produce nans in case of zeros since expre is finite
        re = expre * cosix;
        im = expre * sinix;
      }
    } else {
      re = NaN;
      im = NaN;
    }
  }
  return new Complex(re,im);
}

Complex.prototype.exp = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  if (xre !== xre) {
    // exp(NaN + 0i) = NaN + 0i
    // exp(NaN - 0i) = NaN - 0i
    // exp(NaN + non-zero i) = NaN + NaN i
    re = NaN;
    im = xim === 0 ? xim : NaN;
  } else if (xre === Infinity) {
    if (isFinite(xim)) {
      if (xim === 0) {
        // exp(+Infinity + 0i) = +Infinity + 0i
        // exp(+Infinity - 0i) = +Infinity - 0i
        re = Infinity;
        im = xim;
      } else {
        if (xim > M_MIN || xim < -M_MIN) {
          // re = +Infinity, |im| > M_MIN = infinities with signs of cos/sin
          re = Infinity/Math.cos(xim);
          im = Infinity/Math.sin(xim);
        } else {
          // same as above but |im| < M_MIN, meaning sin/cos is linear
          re = Infinity;
          im = Infinity/xim;
        }
      }
    } else {
      // re = +Infinity, im = NaN or +/-Infinity = Infinity + NaN*i
      re = Infinity;
      im = NaN;
    }
  } else if (xre === -Infinity) {
    if (isFinite(xim)) {
      if (xim === 0) {
        // exp(-Infinity + 0i) = +0 + 0i
        // exp(-Infinity - 0i) = +0 - 0i
        re = 0;
        im = xim;
      } else {
        if (xim > M_MIN || xim < -M_MIN) {
          // re = -Infinity, |im| > M_MIN = zeros with signs of cos/sin
          re = 0*Math.cos(xim);
          im = 0*Math.sin(xim);
        } else {
          // same as above but |im| < M_MIN, meaning sin/cos is linear
          re = 0;
          im = 0*xim;
        }
      }
    } else {
      // re = -Infinity, im = NaN or +/-Infinity = 0 + NaN*i
      re = 0;
      im = NaN;
    }
  } else {
    if (isFinite(xim)) {
      var sinix;
      var cosix;
      if (xim > M_MIN || xim < -M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        // sin/cos is linear
        sinix = xim;
        cosix = 1;
      }
      if (xre > 709) {
        // very small chance this prevents overflow in the output
        xre -= 709;
        sinix *= M_EXP709;
        cosix *= M_EXP709;
      }
      if (xre > 709) {
        // these could be infinity or zero, * M_MAX preserve the signs
        re = M_MAX * cosix;
        im = M_MAX * sinix;
      } else {
        // otherwise Math.exp(xre) wont overflow,
        var expre = Math.exp(xre);
        // although these may still overflow, but will have the correct signs
        // and wont produce nans in case of zeros since expre is finite
        re = expre * cosix;
        im = expre * sinix;
      }
    } else {
      re = NaN;
      im = NaN;
    }
  }
  return new Complex(re,im);
}

Complex.prototype.expEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  if (xre !== xre) {
    // exp(NaN + 0i) = NaN + 0i
    // exp(NaN - 0i) = NaN - 0i
    // exp(NaN + non-zero i) = NaN + NaN i
    this.re = NaN;
    this.im = xim === 0 ? xim : NaN;
  } else if (xre === Infinity) {
    if (isFinite(xim)) {
      if (xim === 0) {
        // exp(+Infinity + 0i) = +Infinity + 0i
        // exp(+Infinity - 0i) = +Infinity - 0i
        this.re = Infinity;
        this.im = xim;
      } else {
        if (xim > M_MIN || xim < -M_MIN) {
          // re = +Infinity, |im| > M_MIN = infinities with signs of cos/sin
          this.re = Infinity/Math.cos(xim);
          this.im = Infinity/Math.sin(xim);
        } else {
          // same as above but |im| < M_MIN, meaning sin/cos is linear
          this.re = Infinity;
          this.im = Infinity/xim;
        }
      }
    } else {
      // re = +Infinity, im = NaN or +/-Infinity = Infinity + NaN*i
      this.re = Infinity;
      this.im = NaN;
    }
  } else if (xre === -Infinity) {
    if (isFinite(xim)) {
      if (xim === 0) {
        // exp(-Infinity + 0i) = +0 + 0i
        // exp(-Infinity - 0i) = +0 - 0i
        this.re = 0;
        this.im = xim;
      } else {
        if (xim > M_MIN || xim < -M_MIN) {
          // re = -Infinity, |im| > M_MIN = zeros with signs of cos/sin
          this.re = 0*Math.cos(xim);
          this.im = 0*Math.sin(xim);
        } else {
          // same as above but |im| < M_MIN, meaning sin/cos is linear
          this.re = 0;
          this.im = 0*xim;
        }
      }
    } else {
      // re = -Infinity, im = NaN or +/-Infinity = 0 + NaN*i
      this.re = 0;
      this.im = NaN;
    }
  } else {
    if (isFinite(xim)) {
      var sinix;
      var cosix;
      if (xim > M_MIN || xim < -M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        // sin/cos is linear
        sinix = xim;
        cosix = 1;
      }
      if (xre > 709) {
        // very small chance this prevents overflow in the output
        xre -= 709;
        sinix *= M_EXP709;
        cosix *= M_EXP709;
      }
      if (xre > 709) {
        // these could be infinity or zero, * M_MAX preserve the signs
        this.re = M_MAX * cosix;
        this.im = M_MAX * sinix;
      } else {
        // otherwise Math.exp(xre) wont overflow,
        var expre = Math.exp(xre);
        // although these may still overflow, but will have the correct signs
        // and wont produce nans in case of zeros since expre is finite
        this.re = expre * cosix;
        this.im = expre * sinix;
      }
    } else {
      this.re = NaN;
      this.im = NaN;
    }
  }
  return this;
}

// take the natural logarithm of a complex number
// static log(z)
// log()
// logEq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 99.455%
// off by 1ULP: 0.5438%
// off by 2ULP: 0.00115%
// off by 3ULP: 0.00005%
// off by 4ULP+: 0%
//
// average error: +/- 0.0054625 ULP

Complex.log = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  if (xre !== xre || xim !== xim) {
    re = NaN;
    im = NaN;
  } else if (xre === Infinity) {
    // atan2(im,+Infinity)
    // atan2(-Infinity,+Infinity) = -pi/4
    // atan2(-0,+Infinity)        = -0
    // atan2(+0,+Infinity)        = +0
    // atan2(+Infinity,+Infinity) = +pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_PI_4;
    else if (xim === -Infinity) im = -M_PI_4;
    else if (1/xim > 0)         im = 0;
    else                        im = -0;
  } else if (xre === -Infinity) {
    // atan2(im,-Infinity)
    // atan2(-Infinity,-Infinity) = -3pi/4
    // atan2(-0,-Infinity)        = -pi
    // atan2(+0,-Infinity)        = +pi
    // atan2(+Infinity,-Infinity) = +3pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_3PI_4;
    else if (xim === -Infinity) im = -M_3PI_4;
    else if (1/xim > 0)         im = M_PI;
    else                        im = -M_PI;
  } else if (xim === Infinity) {
    // atan2(Infinity,re) = pi/2 for finite re
    re = Infinity;
    im = M_PI_2;
  } else if (xim === -Infinity) {
    // atan2(-Infinity,re) = -pi/2 for finite re
    re = Infinity;
    im = -M_PI_2;
  } else if (xim === 0) {
    if (xre === 0) {
      // atan2(im,re)
      // atan2(+0,+0) = +0
      // atan2(+0,-0) = +pi
      // atan2(-0,+0) = -0
      // atan2(-0,-0) = -pi
      re = -Infinity;
      if (1/xre > 0) im = 1/xim > 0 ? 0 : -0;
      else           im = 1/xim > 0 ? M_PI : -M_PI;
    } else {
      if (xre > 0) im = 1/xim > 0 ? 0 : -0;
      else         im = 1/xim > 0 ? M_PI : -M_PI;
      var absx = Math.abs(xre);
      if (absx > M_MAX_2) {
        // prevent overflow
        // log(x/2) + ln2
        re = Math.log(absx * 0.5) + M_LN2;
      } else if (absx < M_MIN) {
        // prevent underflow
        // log(x*2^53) - 53*ln2
        re = Math.log(absx * M_2_EPSILON) - M_53_LN2;
      } else {
        if (absx === 1) {
          // x = 1
          // log(1) = 0
          re = 0;
        } else if (absx > 0.5 && absx < 2) {
          // 1 < x < 2
          // log1p((x - 1)(x + 1))/2
          re = Math.log1p((absx - 1) * (absx + 1)) * 0.5
        } else {
          // general
          // log(x)
          re = Math.log(absx);
        }
      }
    }
  } else {
    var absx = xre > 0 ? xre : -xre;
    var absy = xim > 0 ? xim : -xim;
    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }
    im = Math.atan2(xim, xre);
    if (absx > M_MAX_2) {
      // prevent overflow
      // log(hypot(x/2, y/2)) + ln2
      absx *= 0.5;
      if (absy >= M_2_MIN) absy *= 0.5;
      else                 absy = 0; // y is negligible
      re = Math.log(hypot(absx, absy)) + M_LN2;
    } else if (absx < M_MIN && absy < M_MIN) {
      // prevent underflow
      // log(hypot(x*2^53, y*2^53)) - 53*ln2
      absx *= M_2_EPSILON;
      absy *= M_2_EPSILON;
      re = Math.log(hypot(absx, absy)) - M_53_LN2;
    } else {
      if (absx === 1) {
        // x = 1
        // log1p(y^2)/2
        re = Math.log1p(absy * absy) * 0.5;
      } else if (absx > 1 && absx < 2 && absy < 1) {
        // 1 < x < 2, y < 1
        // log1p((x - 1)(x + 1))/2
        var d2m1 = (absx - 1) * (absx + 1);
        // log1p((x - 1)(x + 1) + y^2)/2
        if (absy >= M_EPSILON) d2m1 += absy * absy;
        re = Math.log1p(d2m1) / 2;
      } else if (absx < 1
              && absx >= 0.5
              && absy < M_EPSILON_2) {
        // 0.5 < x < 1, y is negligible
        // log1p((x - 1)(x + 1))/2
        re = Math.log1p((absx - 1) * (absx + 1)) * 0.5;
      } else if (absx < 1
              && absx >= 0.5
              && absx * absx + absy * absy >= 0.5) {
        // 0.5 < x < 1, y is not negligible
        // log1p(x^2 + y^2 - 1)/2
        re = Math.log1p(x2y2m1(absx, absy)) * 0.5;
      } else {
        // general
        // log(hypot(x, y))
        re = Math.log(hypot(absx, absy));
      }
    }
  }
  return new Complex(re,im);
}

Complex.prototype.log = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  if (xre !== xre || xim !== xim) {
    re = NaN;
    im = NaN;
  } else if (xre === Infinity) {
    // atan2(im,+Infinity)
    // atan2(-Infinity,+Infinity) = -pi/4
    // atan2(-0,+Infinity)        = -0
    // atan2(+0,+Infinity)        = +0
    // atan2(+Infinity,+Infinity) = +pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_PI_4;
    else if (xim === -Infinity) im = -M_PI_4;
    else if (1/xim > 0)         im = 0;
    else                        im = -0;
  } else if (xre === -Infinity) {
    // atan2(im,-Infinity)
    // atan2(-Infinity,-Infinity) = -3pi/4
    // atan2(-0,-Infinity)        = -pi
    // atan2(+0,-Infinity)        = +pi
    // atan2(+Infinity,-Infinity) = +3pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_3PI_4;
    else if (xim === -Infinity) im = -M_3PI_4;
    else if (1/xim > 0)         im = M_PI;
    else                        im = -M_PI;
  } else if (xim === Infinity) {
    // atan2(Infinity,re) = pi/2 for finite re
    re = Infinity;
    im = M_PI_2;
  } else if (xim === -Infinity) {
    // atan2(-Infinity,re) = -pi/2 for finite re
    re = Infinity;
    im = -M_PI_2;
  } else if (xim === 0) {
    if (xre === 0) {
      // atan2(im,re)
      // atan2(+0,+0) = +0
      // atan2(+0,-0) = +pi
      // atan2(-0,+0) = -0
      // atan2(-0,-0) = -pi
      re = -Infinity;
      if (1/xre > 0) im = 1/xim > 0 ? 0 : -0;
      else           im = 1/xim > 0 ? M_PI : -M_PI;
    } else {
      if (xre > 0) im = 1/xim > 0 ? 0 : -0;
      else         im = 1/xim > 0 ? M_PI : -M_PI;
      var absx = Math.abs(xre);
      if (absx > M_MAX_2) {
        // prevent overflow
        // log(x/2) + ln2
        re = Math.log(absx * 0.5) + M_LN2;
      } else if (absx < M_MIN) {
        // prevent underflow
        // log(x*2^53) - 53*ln2
        re = Math.log(absx * M_2_EPSILON) - M_53_LN2;
      } else {
        if (absx === 1) {
          // x = 1
          // log(1) = 0
          re = 0;
        } else if (absx > 0.5 && absx < 2) {
          // 1 < x < 2
          // log1p((x - 1)(x + 1))/2
          re = Math.log1p((absx - 1) * (absx + 1)) * 0.5
        } else {
          // general
          // log(x)
          re = Math.log(absx);
        }
      }
    }
  } else {
    var absx = xre > 0 ? xre : -xre;
    var absy = xim > 0 ? xim : -xim;
    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }
    im = Math.atan2(xim, xre);
    if (absx > M_MAX_2) {
      // prevent overflow
      // log(hypot(x/2, y/2)) + ln2
      absx *= 0.5;
      if (absy >= M_2_MIN) absy *= 0.5;
      else                 absy = 0; // y is negligible
      re = Math.log(hypot(absx, absy)) + M_LN2;
    } else if (absx < M_MIN && absy < M_MIN) {
      // prevent underflow
      // log(hypot(x*2^53, y*2^53)) - 53*ln2
      absx *= M_2_EPSILON;
      absy *= M_2_EPSILON;
      re = Math.log(hypot(absx, absy)) - M_53_LN2;
    } else {
      if (absx === 1) {
        // x = 1
        // log1p(y^2)/2
        re = Math.log1p(absy * absy) * 0.5;
      } else if (absx > 1 && absx < 2 && absy < 1) {
        // 1 < x < 2, y < 1
        // log1p((x - 1)(x + 1))/2
        var d2m1 = (absx - 1) * (absx + 1);
        // log1p((x - 1)(x + 1) + y^2)/2
        if (absy >= M_EPSILON) d2m1 += absy * absy;
        re = Math.log1p(d2m1) / 2;
      } else if (absx < 1
              && absx >= 0.5
              && absy < M_EPSILON_2) {
        // 0.5 < x < 1, y is negligible
        // log1p((x - 1)(x + 1))/2
        re = Math.log1p((absx - 1) * (absx + 1)) * 0.5;
      } else if (absx < 1
              && absx >= 0.5
              && absx * absx + absy * absy >= 0.5) {
        // 0.5 < x < 1, y is not negligible
        // log1p(x^2 + y^2 - 1)/2
        re = Math.log1p(x2y2m1(absx, absy)) * 0.5;
      } else {
        // general
        // log(hypot(x, y))
        re = Math.log(hypot(absx, absy));
      }
    }
  }
  return new Complex(re,im);
}

Complex.prototype.logEq = function() {
  var xre = this.re;
  var xim = this.im;
  if (xre !== xre || xim !== xim) {
    this.re = NaN;
    this.im = NaN;
  } else if (xre === Infinity) {
    // atan2(im,+Infinity)
    // atan2(-Infinity,+Infinity) = -pi/4
    // atan2(-0,+Infinity)        = -0
    // atan2(+0,+Infinity)        = +0
    // atan2(+Infinity,+Infinity) = +pi/4
    this.re = Infinity;
    if (xim === Infinity)       this.im = M_PI_4;
    else if (xim === -Infinity) this.im = -M_PI_4;
    else if (1/xim > 0)         this.im = 0;
    else                        this.im = -0;
  } else if (xre === -Infinity) {
    // atan2(im,-Infinity)
    // atan2(-Infinity,-Infinity) = -3pi/4
    // atan2(-0,-Infinity)        = -pi
    // atan2(+0,-Infinity)        = +pi
    // atan2(+Infinity,-Infinity) = +3pi/4
    this.re = Infinity;
    if (xim === Infinity)       this.im = M_3PI_4;
    else if (xim === -Infinity) this.im = -M_3PI_4;
    else if (1/xim > 0)         this.im = M_PI;
    else                        this.im = -M_PI;
  } else if (xim === Infinity) {
    // atan2(Infinity,re) = pi/2 for finite re
    this.re = Infinity;
    this.im = M_PI_2;
  } else if (xim === -Infinity) {
    // atan2(-Infinity,re) = -pi/2 for finite re
    this.re = Infinity;
    this.im = -M_PI_2;
  } else if (xim === 0) {
    if (xre === 0) {
      // atan2(im,re)
      // atan2(+0,+0) = +0
      // atan2(+0,-0) = +pi
      // atan2(-0,+0) = -0
      // atan2(-0,-0) = -pi
      this.re = -Infinity;
      if (1/xre > 0) this.im = 1/xim > 0 ? 0 : -0;
      else           this.im = 1/xim > 0 ? M_PI : -M_PI;
    } else {
      if (xre > 0) this.im = 1/xim > 0 ? 0 : -0;
      else         this.im = 1/xim > 0 ? M_PI : -M_PI;
      var absx = Math.abs(xre);
      if (absx > M_MAX_2) {
        // prevent overflow
        // log(x/2) + ln2
        this.re = Math.log(absx * 0.5) + M_LN2;
      } else if (absx < M_MIN) {
        // prevent underflow
        // log(x*2^53) - 53*ln2
        this.re = Math.log(absx * M_2_EPSILON) - M_53_LN2;
      } else {
        if (absx === 1) {
          // x = 1
          // log(1) = 0
          this.re = 0;
        } else if (absx > 0.5 && absx < 2) {
          // 1 < x < 2
          // log1p((x - 1)(x + 1))/2
          this.re = Math.log1p((absx - 1) * (absx + 1)) * 0.5
        } else {
          // general
          // log(x)
          this.re = Math.log(absx);
        }
      }
    }
  } else {
    var absx = xre > 0 ? xre : -xre;
    var absy = xim > 0 ? xim : -xim;
    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }
    this.im = Math.atan2(xim, xre);
    if (absx > M_MAX_2) {
      // prevent overflow
      // log(hypot(x/2, y/2)) + ln2
      absx *= 0.5;
      if (absy >= M_2_MIN) absy *= 0.5;
      else                 absy = 0; // y is negligible
      this.re = Math.log(hypot(absx, absy)) + M_LN2;
    } else if (absx < M_MIN && absy < M_MIN) {
      // prevent underflow
      // log(hypot(x*2^53, y*2^53)) - 53*ln2
      absx *= M_2_EPSILON;
      absy *= M_2_EPSILON;
      this.re = Math.log(hypot(absx, absy)) - M_53_LN2;
    } else {
      if (absx === 1) {
        // x = 1
        // log1p(y^2)/2
        this.re = Math.log1p(absy * absy) * 0.5;
      } else if (absx > 1 && absx < 2 && absy < 1) {
        // 1 < x < 2, y < 1
        // log1p((x - 1)(x + 1))/2
        var d2m1 = (absx - 1) * (absx + 1);
        // log1p((x - 1)(x + 1) + y^2)/2
        if (absy >= M_EPSILON) d2m1 += absy * absy;
        this.re = Math.log1p(d2m1) / 2;
      } else if (absx < 1
              && absx >= 0.5
              && absy < M_EPSILON_2) {
        // 0.5 < x < 1, y is negligible
        // log1p((x - 1)(x + 1))/2
        this.re = Math.log1p((absx - 1) * (absx + 1)) * 0.5;
      } else if (absx < 1
              && absx >= 0.5
              && absx * absx + absy * absy >= 0.5) {
        // 0.5 < x < 1, y is not negligible
        // log1p(x^2 + y^2 - 1)/2
        this.re = Math.log1p(x2y2m1(absx, absy)) * 0.5;
      } else {
        // general
        // log(hypot(x, y))
        this.re = Math.log(hypot(absx, absy));
      }
    }
  }
  return this;
}

// take the base-10 logarithm of a complex number
// static log10(z)
// log10()
// log10Eq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 61.0621%
// off by 1ULP: 38.8279%
// off by 2ULP: 0.00995%
// off by 3ULP: 0.00005%
// off by 4ULP+: 0%
//
// average error: +/- 0.3884795 ULP

Complex.log10 = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  if (xre !== xre || xim !== xim) {
    re = NaN;
    im = NaN;
  } else if (xre === Infinity) {
    // atan2(im,+Infinity)
    // atan2(-Infinity,+Infinity) = -pi/4
    // atan2(-0,+Infinity)        = -0
    // atan2(+0,+Infinity)        = +0
    // atan2(+Infinity,+Infinity) = +pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_PI_4_LN10;
    else if (xim === -Infinity) im = -M_PI_4_LN10;
    else if (1/xim > 0)         im = 0;
    else                        im = -0;
  } else if (xre === -Infinity) {
    // atan2(im,-Infinity)
    // atan2(-Infinity,-Infinity) = -3pi/4
    // atan2(-0,-Infinity)        = -pi
    // atan2(+0,-Infinity)        = +pi
    // atan2(+Infinity,-Infinity) = +3pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_3PI4_LN10;
    else if (xim === -Infinity) im = -M_3PI4_LN10;
    else if (1/xim > 0)         im = M_PI_LN10;
    else                        im = -M_PI_LN10;
  } else if (xim === Infinity) {
    // atan2(Infinity,re) = pi/2 for finite re
    re = Infinity;
    im = M_PI_2_LN10;
  } else if (xim === -Infinity) {
    // atan2(-Infinity,re) = -pi/2 for finite re
    re = Infinity;
    im = -M_PI_2_LN10;
  } else if (xim === 0) {
    if (xre === 0) {
      // atan2(im,re)
      // atan2(+0,+0) = +0
      // atan2(+0,-0) = +pi
      // atan2(-0,+0) = -0
      // atan2(-0,-0) = -pi
      re = -Infinity;
      if (1/xre > 0) im = 1/xim > 0 ? 0 : -0;
      else           im = 1/xim > 0 ? M_PI_LN10 : -M_PI_LN10;
    } else {
      if (xre > 0) im = 1/xim > 0 ? 0 : -0;
      else         im = 1/xim > 0 ? M_PI_LN10 : -M_PI_LN10;
      var absx = Math.abs(xre);
      if (absx > M_MAX_2) {
        // prevent overflow
        // log10(x/2) + log10(2)
        re = Math.log10(absx * 0.5) + M_LOG2;
      } else if (absx < M_MIN) {
        // prevent underflow
        // log10(x*2^53) - 53*log10(2)
        re = Math.log10(absx * M_2_EPSILON) - M_53_LOG2;
      } else {
        if (absx === 1) {
          // x = 1
          // log10(1) = 0
          re = 0;
        } else if (absx > 0.5 && absx < 2) {
          // 1 < x < 2
          // log1p((x - 1)(x + 1))/2 / ln10
          re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN10;
        } else {
          // general
          // log10(x)
          re = Math.log10(absx);
        }
      }
    }
  } else {
    var absx = xre > 0 ? xre : -xre;
    var absy = xim > 0 ? xim : -xim;
    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }
    im = Math.atan2(xim, xre) * M_1_LN10;
    if (absx > M_MAX_2) {
      // prevent overflow
      // log10(hypot(x/2, y/2)) + log10(2)
      absx *= 0.5;
      if (absy >= M_2_MIN) absy *= 0.5;
      else                 absy = 0; // y is negligible
      re = Math.log10(hypot(absx, absy)) + M_LOG2;
    } else if (absx < M_MIN && absy < M_MIN) {
      // prevent underflow
      // log10(hypot(x*2^53, y*2^53)) - 53*log10(2)
      absx *= M_2_EPSILON;
      absy *= M_2_EPSILON;
      re = Math.log10(hypot(absx, absy)) - M_53_LOG2;
    } else {
      if (absx === 1) {
        // x = 1
        // log1p(y^2)/2 / ln10
        re = Math.log1p(absy * absy) * M_1_2_LN10;
      } else if (absx > 1 && absx < 2 && absy < 1) {
        // 1 < x < 2, y < 1
        // log1p((x - 1)(x + 1))/2 / ln10
        var d2m1 = (absx - 1) * (absx + 1);
        // log1p((x - 1)(x + 1) + y^2)/2
        if (absy >= M_EPSILON) d2m1 += absy * absy;
        re = Math.log1p(d2m1) * M_1_2_LN10;
      } else if (absx < 1
              && absx >= 0.5
              && absy < M_EPSILON_2) {
        // 0.5 < x < 1, y is negligible
        // log1p((x - 1)(x + 1))/2 / ln10
        re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN10;
      } else if (absx < 1
              && absx >= 0.5
              && absx * absx + absy * absy >= 0.5) {
        // 0.5 < x < 1, y is not negligible
        // log1p(x^2 + y^2 - 1)/2 / ln10
        re = Math.log1p(x2y2m1(absx, absy)) * M_1_2_LN10;
      } else {
        // general
        // log10(hypot(x, y))
        re = Math.log10(hypot(absx, absy));
      }
    }
  }
  return new Complex(re,im);
}

Complex.prototype.log10 = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  if (xre !== xre || xim !== xim) {
    re = NaN;
    im = NaN;
  } else if (xre === Infinity) {
    // atan2(im,+Infinity)
    // atan2(-Infinity,+Infinity) = -pi/4
    // atan2(-0,+Infinity)        = -0
    // atan2(+0,+Infinity)        = +0
    // atan2(+Infinity,+Infinity) = +pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_PI_4_LN10;
    else if (xim === -Infinity) im = -M_PI_4_LN10;
    else if (1/xim > 0)         im = 0;
    else                        im = -0;
  } else if (xre === -Infinity) {
    // atan2(im,-Infinity)
    // atan2(-Infinity,-Infinity) = -3pi/4
    // atan2(-0,-Infinity)        = -pi
    // atan2(+0,-Infinity)        = +pi
    // atan2(+Infinity,-Infinity) = +3pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_3PI4_LN10;
    else if (xim === -Infinity) im = -M_3PI4_LN10;
    else if (1/xim > 0)         im = M_PI_LN10;
    else                        im = -M_PI_LN10;
  } else if (xim === Infinity) {
    // atan2(Infinity,re) = pi/2 for finite re
    re = Infinity;
    im = M_PI_2_LN10;
  } else if (xim === -Infinity) {
    // atan2(-Infinity,re) = -pi/2 for finite re
    re = Infinity;
    im = -M_PI_2_LN10;
  } else if (xim === 0) {
    if (xre === 0) {
      // atan2(im,re)
      // atan2(+0,+0) = +0
      // atan2(+0,-0) = +pi
      // atan2(-0,+0) = -0
      // atan2(-0,-0) = -pi
      re = -Infinity;
      if (1/xre > 0) im = 1/xim > 0 ? 0 : -0;
      else           im = 1/xim > 0 ? M_PI_LN10 : -M_PI_LN10;
    } else {
      if (xre > 0) im = 1/xim > 0 ? 0 : -0;
      else         im = 1/xim > 0 ? M_PI_LN10 : -M_PI_LN10;
      var absx = Math.abs(xre);
      if (absx > M_MAX_2) {
        // prevent overflow
        // log10(x/2) + log10(2)
        re = Math.log10(absx * 0.5) + M_LOG2;
      } else if (absx < M_MIN) {
        // prevent underflow
        // log10(x*2^53) - 53*log10(2)
        re = Math.log10(absx * M_2_EPSILON) - M_53_LOG2;
      } else {
        if (absx === 1) {
          // x = 1
          // log10(1) = 0
          re = 0;
        } else if (absx > 0.5 && absx < 2) {
          // 1 < x < 2
          // log1p((x - 1)(x + 1))/2 / ln10
          re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN10;
        } else {
          // general
          // log10(x)
          re = Math.log10(absx);
        }
      }
    }
  } else {
    var absx = xre > 0 ? xre : -xre;
    var absy = xim > 0 ? xim : -xim;
    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }
    im = Math.atan2(xim, xre) * M_1_LN10;
    if (absx > M_MAX_2) {
      // prevent overflow
      // log10(hypot(x/2, y/2)) + log10(2)
      absx *= 0.5;
      if (absy >= M_2_MIN) absy *= 0.5;
      else                 absy = 0; // y is negligible
      re = Math.log10(hypot(absx, absy)) + M_LOG2;
    } else if (absx < M_MIN && absy < M_MIN) {
      // prevent underflow
      // log10(hypot(x*2^53, y*2^53)) - 53*log10(2)
      absx *= M_2_EPSILON;
      absy *= M_2_EPSILON;
      re = Math.log10(hypot(absx, absy)) - M_53_LOG2;
    } else {
      if (absx === 1) {
        // x = 1
        // log1p(y^2)/2 / ln10
        re = Math.log1p(absy * absy) * M_1_2_LN10;
      } else if (absx > 1 && absx < 2 && absy < 1) {
        // 1 < x < 2, y < 1
        // log1p((x - 1)(x + 1))/2 / ln10
        var d2m1 = (absx - 1) * (absx + 1);
        // log1p((x - 1)(x + 1) + y^2)/2
        if (absy >= M_EPSILON) d2m1 += absy * absy;
        re = Math.log1p(d2m1) * M_1_2_LN10;
      } else if (absx < 1
              && absx >= 0.5
              && absy < M_EPSILON_2) {
        // 0.5 < x < 1, y is negligible
        // log1p((x - 1)(x + 1))/2 / ln10
        re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN10;
      } else if (absx < 1
              && absx >= 0.5
              && absx * absx + absy * absy >= 0.5) {
        // 0.5 < x < 1, y is not negligible
        // log1p(x^2 + y^2 - 1)/2 / ln10
        re = Math.log1p(x2y2m1(absx, absy)) * M_1_2_LN10;
      } else {
        // general
        // log10(hypot(x, y))
        re = Math.log10(hypot(absx, absy));
      }
    }
  }
  return new Complex(re,im);
}

Complex.prototype.log10Eq = function() {
  var xre = this.re;
  var xim = this.im;
  if (xre !== xre || xim !== xim) {
    this.re = NaN;
    this.im = NaN;
  } else if (xre === Infinity) {
    // atan2(im,+Infinity)
    // atan2(-Infinity,+Infinity) = -pi/4
    // atan2(-0,+Infinity)        = -0
    // atan2(+0,+Infinity)        = +0
    // atan2(+Infinity,+Infinity) = +pi/4
    this.re = Infinity;
    if (xim === Infinity)       this.im = M_PI_4_LN10;
    else if (xim === -Infinity) this.im = -M_PI_4_LN10;
    else if (1/xim > 0)         this.im = 0;
    else                        this.im = -0;
  } else if (xre === -Infinity) {
    // atan2(im,-Infinity)
    // atan2(-Infinity,-Infinity) = -3pi/4
    // atan2(-0,-Infinity)        = -pi
    // atan2(+0,-Infinity)        = +pi
    // atan2(+Infinity,-Infinity) = +3pi/4
    this.re = Infinity;
    if (xim === Infinity)       this.im = M_3PI4_LN10;
    else if (xim === -Infinity) this.im = -M_3PI4_LN10;
    else if (1/xim > 0)         this.im = M_PI_LN10;
    else                        this.im = -M_PI_LN10;
  } else if (xim === Infinity) {
    // atan2(Infinity,re) = pi/2 for finite re
    this.re = Infinity;
    this.im = M_PI_2_LN10;
  } else if (xim === -Infinity) {
    // atan2(-Infinity,re) = -pi/2 for finite re
    this.re = Infinity;
    this.im = -M_PI_2_LN10;
  } else if (xim === 0) {
    if (xre === 0) {
      // atan2(im,re)
      // atan2(+0,+0) = +0
      // atan2(+0,-0) = +pi
      // atan2(-0,+0) = -0
      // atan2(-0,-0) = -pi
      this.re = -Infinity;
      if (1/xre > 0) this.im = 1/xim > 0 ? 0 : -0;
      else           this.im = 1/xim > 0 ? M_PI_LN10 : -M_PI_LN10;
    } else {
      if (xre > 0) this.im = 1/xim > 0 ? 0 : -0;
      else         this.im = 1/xim > 0 ? M_PI_LN10 : -M_PI_LN10;
      var absx = Math.abs(xre);
      if (absx > M_MAX_2) {
        // prevent overflow
        // log10(x/2) + log10(2)
        this.re = Math.log10(absx * 0.5) + M_LOG2;
      } else if (absx < M_MIN) {
        // prevent underflow
        // log10(x*2^53) - 53*log10(2)
        this.re = Math.log10(absx * M_2_EPSILON) - M_53_LOG2;
      } else {
        if (absx === 1) {
          // x = 1
          // log10(1) = 0
          this.re = 0;
        } else if (absx > 0.5 && absx < 2) {
          // 1 < x < 2
          // log1p((x - 1)(x + 1))/2 / ln10
          this.re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN10;
        } else {
          // general
          // log10(x)
          this.re = Math.log10(absx);
        }
      }
    }
  } else {
    var absx = xre > 0 ? xre : -xre;
    var absy = xim > 0 ? xim : -xim;
    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }
    this.im = Math.atan2(xim, xre) * M_1_LN10;
    if (absx > M_MAX_2) {
      // prevent overflow
      // log10(hypot(x/2, y/2)) + log10(2)
      absx *= 0.5;
      if (absy >= M_2_MIN) absy *= 0.5;
      else                 absy = 0; // y is negligible
      this.re = Math.log10(hypot(absx, absy)) + M_LOG2;
    } else if (absx < M_MIN && absy < M_MIN) {
      // prevent underflow
      // log10(hypot(x*2^53, y*2^53)) - 53*log10(2)
      absx *= M_2_EPSILON;
      absy *= M_2_EPSILON;
      this.re = Math.log10(hypot(absx, absy)) - M_53_LOG2;
    } else {
      if (absx === 1) {
        // x = 1
        // log1p(y^2)/2 / ln10
        this.re = Math.log1p(absy * absy) * M_1_2_LN10;
      } else if (absx > 1 && absx < 2 && absy < 1) {
        // 1 < x < 2, y < 1
        // log1p((x - 1)(x + 1))/2 / ln10
        var d2m1 = (absx - 1) * (absx + 1);
        // log1p((x - 1)(x + 1) + y^2)/2
        if (absy >= M_EPSILON) d2m1 += absy * absy;
        this.re = Math.log1p(d2m1) * M_1_2_LN10;
      } else if (absx < 1
              && absx >= 0.5
              && absy < M_EPSILON_2) {
        // 0.5 < x < 1, y is negligible
        // log1p((x - 1)(x + 1))/2 / ln10
        this.re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN10;
      } else if (absx < 1
              && absx >= 0.5
              && absx * absx + absy * absy >= 0.5) {
        // 0.5 < x < 1, y is not negligible
        // log1p(x^2 + y^2 - 1)/2 / ln10
        this.re = Math.log1p(x2y2m1(absx, absy)) * M_1_2_LN10;
      } else {
        // general
        // log10(hypot(x, y))
        this.re = Math.log10(hypot(absx, absy));
      }
    }
  }
  return this;
}

// take the base-2 logarithm of a complex number
// static log2(z)
// log2()
// log2Eq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 61.3282%
// off by 1ULP: 38.56645%
// off by 2ULP: 0.0052%
// off by 3ULP: 0.00015%
// off by 4ULP+: 0%
//
// average error: +/- 0.385773 ULP

Complex.log2 = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  if (xre !== xre || xim !== xim) {
    re = NaN;
    im = NaN;
  } else if (xre === Infinity) {
    // atan2(im,+Infinity)
    // atan2(-Infinity,+Infinity) = -pi/4
    // atan2(-0,+Infinity)        = -0
    // atan2(+0,+Infinity)        = +0
    // atan2(+Infinity,+Infinity) = +pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_PI_4_LN2;
    else if (xim === -Infinity) im = -M_PI_4_LN2;
    else if (1/xim > 0)         im = 0;
    else                        im = -0;
  } else if (xre === -Infinity) {
    // atan2(im,-Infinity)
    // atan2(-Infinity,-Infinity) = -3pi/4
    // atan2(-0,-Infinity)        = -pi
    // atan2(+0,-Infinity)        = +pi
    // atan2(+Infinity,-Infinity) = +3pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_3PI4_LN2;
    else if (xim === -Infinity) im = -M_3PI4_LN2;
    else if (1/xim > 0)         im = M_PI_LN2;
    else                        im = -M_PI_LN2;
  } else if (xim === Infinity) {
    // atan2(Infinity,re) = pi/2 for finite re
    re = Infinity;
    im = M_PI_2_LN2;
  } else if (xim === -Infinity) {
    // atan2(-Infinity,re) = -pi/2 for finite re
    re = Infinity;
    im = -M_PI_2_LN2;
  } else if (xim === 0) {
    if (xre === 0) {
      // atan2(im,re)
      // atan2(+0,+0) = +0
      // atan2(+0,-0) = +pi
      // atan2(-0,+0) = -0
      // atan2(-0,-0) = -pi
      re = -Infinity;
      if (1/xre > 0) im = 1/xim > 0 ? 0 : -0;
      else           im = 1/xim > 0 ? M_PI_LN2 : -M_PI_LN2;
    } else {
      if (xre > 0) im = 1/xim > 0 ? 0 : -0;
      else         im = 1/xim > 0 ? M_PI_LN2 : -M_PI_LN2;
      var absx = Math.abs(xre);
      if (absx > M_MAX_2) {
        // prevent overflow
        // log2(x/2) + log2(2) = log2(x/2) + 1
        re = Math.log2(absx * 0.5) + 1;
      } else if (absx < M_MIN) {
        // prevent underflow
        // log2(x*2^53) - 53*log2(2) = log2(x*2^53) - 53 
        re = Math.log2(absx * M_2_EPSILON) - 53;
      } else {
        if (absx === 1) {
          // x = 1
          // log2(1) = 0
          re = 0;
        } else if (absx > 0.5 && absx < 2) {
          // 1 < x < 2
          // log1p((x - 1)(x + 1))/2 / ln2
          re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN2;
        } else {
          // general
          // log2(x)
          re = Math.log2(absx);
        }
      }
    }
  } else {
    var absx = xre > 0 ? xre : -xre;
    var absy = xim > 0 ? xim : -xim;
    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }
    im = Math.atan2(xim, xre) * M_1_LN2;
    if (absx > M_MAX_2) {
      // prevent overflow
      // log2(hypot(x/2, y/2)) + log2(2) = log2(hypot(x/2, y/2)) + 1
      absx *= 0.5;
      if (absy >= M_2_MIN) absy *= 0.5;
      else                 absy = 0; // y is negligible
      re = Math.log2(hypot(absx, absy)) + 1;
    } else if (absx < M_MIN && absy < M_MIN) {
      // prevent underflow
      // log2(hypot(x*2^53, y*2^53)) - 53*log2(2)
      // = log2(hypot(x*2^53, y*2^53)) - 53
      absx *= M_2_EPSILON;
      absy *= M_2_EPSILON;
      re = Math.log2(hypot(absx, absy)) - 53;
    } else {
      if (absx === 1) {
        // x = 1
        // log1p(y^2)/2 / ln2
        re = Math.log1p(absy * absy) * M_1_2_LN2;
      } else if (absx > 1 && absx < 2 && absy < 1) {
        // 1 < x < 2, y < 1
        // log1p((x - 1)(x + 1))/2 / ln2
        var d2m1 = (absx - 1) * (absx + 1);
        // log1p((x - 1)(x + 1) + y^2)/2 / ln2
        if (absy >= M_EPSILON) d2m1 += absy * absy;
        re = Math.log1p(d2m1) * M_1_2_LN2;
      } else if (absx < 1
              && absx >= 0.5
              && absy < M_EPSILON_2) {
        // 0.5 < x < 1, y is negligible
        // log1p((x - 1)(x + 1))/2 / ln2
        re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN2;
      } else if (absx < 1
              && absx >= 0.5
              && absx * absx + absy * absy >= 0.5) {
        // 0.5 < x < 1, y is not negligible
        // log1p(x^2 + y^2 - 1)/2 / ln2
        re = Math.log1p(x2y2m1(absx, absy)) * M_1_2_LN2;
      } else {
        // general
        // log2(hypot(x, y))
        re = Math.log2(hypot(absx, absy));
      }
    }
  }
  return new Complex(re,im);
}

Complex.prototype.log2 = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  if (xre !== xre || xim !== xim) {
    re = NaN;
    im = NaN;
  } else if (xre === Infinity) {
    // atan2(im,+Infinity)
    // atan2(-Infinity,+Infinity) = -pi/4
    // atan2(-0,+Infinity)        = -0
    // atan2(+0,+Infinity)        = +0
    // atan2(+Infinity,+Infinity) = +pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_PI_4_LN2;
    else if (xim === -Infinity) im = -M_PI_4_LN2;
    else if (1/xim > 0)         im = 0;
    else                        im = -0;
  } else if (xre === -Infinity) {
    // atan2(im,-Infinity)
    // atan2(-Infinity,-Infinity) = -3pi/4
    // atan2(-0,-Infinity)        = -pi
    // atan2(+0,-Infinity)        = +pi
    // atan2(+Infinity,-Infinity) = +3pi/4
    re = Infinity;
    if (xim === Infinity)       im = M_3PI4_LN2;
    else if (xim === -Infinity) im = -M_3PI4_LN2;
    else if (1/xim > 0)         im = M_PI_LN2;
    else                        im = -M_PI_LN2;
  } else if (xim === Infinity) {
    // atan2(Infinity,re) = pi/2 for finite re
    re = Infinity;
    im = M_PI_2_LN2;
  } else if (xim === -Infinity) {
    // atan2(-Infinity,re) = -pi/2 for finite re
    re = Infinity;
    im = -M_PI_2_LN2;
  } else if (xim === 0) {
    if (xre === 0) {
      // atan2(im,re)
      // atan2(+0,+0) = +0
      // atan2(+0,-0) = +pi
      // atan2(-0,+0) = -0
      // atan2(-0,-0) = -pi
      re = -Infinity;
      if (1/xre > 0) im = 1/xim > 0 ? 0 : -0;
      else           im = 1/xim > 0 ? M_PI_LN2 : -M_PI_LN2;
    } else {
      if (xre > 0) im = 1/xim > 0 ? 0 : -0;
      else         im = 1/xim > 0 ? M_PI_LN2 : -M_PI_LN2;
      var absx = Math.abs(xre);
      if (absx > M_MAX_2) {
        // prevent overflow
        // log2(x/2) + log2(2) = log2(x/2) + 1
        re = Math.log2(absx * 0.5) + 1;
      } else if (absx < M_MIN) {
        // prevent underflow
        // log2(x*2^53) - 53*log2(2) = log2(x*2^53) - 53 
        re = Math.log2(absx * M_2_EPSILON) - 53;
      } else {
        if (absx === 1) {
          // x = 1
          // log2(1) = 0
          re = 0;
        } else if (absx > 0.5 && absx < 2) {
          // 1 < x < 2
          // log1p((x - 1)(x + 1))/2 / ln2
          re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN2;
        } else {
          // general
          // log2(x)
          re = Math.log2(absx);
        }
      }
    }
  } else {
    var absx = xre > 0 ? xre : -xre;
    var absy = xim > 0 ? xim : -xim;
    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }
    im = Math.atan2(xim, xre) * M_1_LN2;
    if (absx > M_MAX_2) {
      // prevent overflow
      // log2(hypot(x/2, y/2)) + log2(2) = log2(hypot(x/2, y/2)) + 1
      absx *= 0.5;
      if (absy >= M_2_MIN) absy *= 0.5;
      else                 absy = 0; // y is negligible
      re = Math.log2(hypot(absx, absy)) + 1;
    } else if (absx < M_MIN && absy < M_MIN) {
      // prevent underflow
      // log2(hypot(x*2^53, y*2^53)) - 53*log2(2)
      // = log2(hypot(x*2^53, y*2^53)) - 53
      absx *= M_2_EPSILON;
      absy *= M_2_EPSILON;
      re = Math.log2(hypot(absx, absy)) - 53;
    } else {
      if (absx === 1) {
        // x = 1
        // log1p(y^2)/2 / ln2
        re = Math.log1p(absy * absy) * M_1_2_LN2;
      } else if (absx > 1 && absx < 2 && absy < 1) {
        // 1 < x < 2, y < 1
        // log1p((x - 1)(x + 1))/2 / ln2
        var d2m1 = (absx - 1) * (absx + 1);
        // log1p((x - 1)(x + 1) + y^2)/2 / ln2
        if (absy >= M_EPSILON) d2m1 += absy * absy;
        re = Math.log1p(d2m1) * M_1_2_LN2;
      } else if (absx < 1
              && absx >= 0.5
              && absy < M_EPSILON_2) {
        // 0.5 < x < 1, y is negligible
        // log1p((x - 1)(x + 1))/2 / ln2
        re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN2;
      } else if (absx < 1
              && absx >= 0.5
              && absx * absx + absy * absy >= 0.5) {
        // 0.5 < x < 1, y is not negligible
        // log1p(x^2 + y^2 - 1)/2 / ln2
        re = Math.log1p(x2y2m1(absx, absy)) * M_1_2_LN2;
      } else {
        // general
        // log2(hypot(x, y))
        re = Math.log2(hypot(absx, absy));
      }
    }
  }
  return new Complex(re,im);
}

Complex.prototype.log2Eq = function() {
  var xre = this.re;
  var xim = this.im;
  if (xre !== xre || xim !== xim) {
    this.re = NaN;
    this.im = NaN;
  } else if (xre === Infinity) {
    // atan2(im,+Infinity)
    // atan2(-Infinity,+Infinity) = -pi/4
    // atan2(-0,+Infinity)        = -0
    // atan2(+0,+Infinity)        = +0
    // atan2(+Infinity,+Infinity) = +pi/4
    this.re = Infinity;
    if (xim === Infinity)       this.im = M_PI_4_LN2;
    else if (xim === -Infinity) this.im = -M_PI_4_LN2;
    else if (1/xim > 0)         this.im = 0;
    else                        this.im = -0;
  } else if (xre === -Infinity) {
    // atan2(im,-Infinity)
    // atan2(-Infinity,-Infinity) = -3pi/4
    // atan2(-0,-Infinity)        = -pi
    // atan2(+0,-Infinity)        = +pi
    // atan2(+Infinity,-Infinity) = +3pi/4
    this.re = Infinity;
    if (xim === Infinity)       this.im = M_3PI4_LN2;
    else if (xim === -Infinity) this.im = -M_3PI4_LN2;
    else if (1/xim > 0)         this.im = M_PI_LN2;
    else                        this.im = -M_PI_LN2;
  } else if (xim === Infinity) {
    // atan2(Infinity,re) = pi/2 for finite re
    this.re = Infinity;
    this.im = M_PI_2_LN2;
  } else if (xim === -Infinity) {
    // atan2(-Infinity,re) = -pi/2 for finite re
    this.re = Infinity;
    this.im = -M_PI_2_LN2;
  } else if (xim === 0) {
    if (xre === 0) {
      // atan2(im,re)
      // atan2(+0,+0) = +0
      // atan2(+0,-0) = +pi
      // atan2(-0,+0) = -0
      // atan2(-0,-0) = -pi
      this.re = -Infinity;
      if (1/xre > 0) this.im = 1/xim > 0 ? 0 : -0;
      else           this.im = 1/xim > 0 ? M_PI_LN2 : -M_PI_LN2;
    } else {
      if (xre > 0) this.im = 1/xim > 0 ? 0 : -0;
      else         this.im = 1/xim > 0 ? M_PI_LN2 : -M_PI_LN2;
      var absx = Math.abs(xre);
      if (absx > M_MAX_2) {
        // prevent overflow
        // log2(x/2) + log2(2) = log2(x/2) + 1
        this.re = Math.log2(absx * 0.5) + 1;
      } else if (absx < M_MIN) {
        // prevent underflow
        // log2(x*2^53) - 53*log2(2) = log2(x*2^53) - 53 
        this.re = Math.log2(absx * M_2_EPSILON) - 53;
      } else {
        if (absx === 1) {
          // x = 1
          // log2(1) = 0
          this.re = 0;
        } else if (absx > 0.5 && absx < 2) {
          // 1 < x < 2
          // log1p((x - 1)(x + 1))/2 / ln2
          this.re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN2;
        } else {
          // general
          // log2(x)
          this.re = Math.log2(absx);
        }
      }
    }
  } else {
    var absx = xre > 0 ? xre : -xre;
    var absy = xim > 0 ? xim : -xim;
    if (absx < absy) {
      var t = absx;
      absx = absy;
      absy = t;
    }
    this.im = Math.atan2(xim, xre) * M_1_LN2;
    if (absx > M_MAX_2) {
      // prevent overflow
      // log2(hypot(x/2, y/2)) + log2(2) = log2(hypot(x/2, y/2)) + 1
      absx *= 0.5;
      if (absy >= M_2_MIN) absy *= 0.5;
      else                 absy = 0; // y is negligible
      this.re = Math.log2(hypot(absx, absy)) + 1;
    } else if (absx < M_MIN && absy < M_MIN) {
      // prevent underflow
      // log2(hypot(x*2^53, y*2^53)) - 53*log2(2)
      // = log2(hypot(x*2^53, y*2^53)) - 53
      absx *= M_2_EPSILON;
      absy *= M_2_EPSILON;
      this.re = Math.log2(hypot(absx, absy)) - 53;
    } else {
      if (absx === 1) {
        // x = 1
        // log1p(y^2)/2 / ln2
        this.re = Math.log1p(absy * absy) * M_1_2_LN2;
      } else if (absx > 1 && absx < 2 && absy < 1) {
        // 1 < x < 2, y < 1
        // log1p((x - 1)(x + 1))/2 / ln2
        var d2m1 = (absx - 1) * (absx + 1);
        // log1p((x - 1)(x + 1) + y^2)/2 / ln2
        if (absy >= M_EPSILON) d2m1 += absy * absy;
        this.re = Math.log1p(d2m1) * M_1_2_LN2;
      } else if (absx < 1
              && absx >= 0.5
              && absy < M_EPSILON_2) {
        // 0.5 < x < 1, y is negligible
        // log1p((x - 1)(x + 1))/2 / ln2
        this.re = Math.log1p((absx - 1) * (absx + 1)) * M_1_2_LN2;
      } else if (absx < 1
              && absx >= 0.5
              && absx * absx + absy * absy >= 0.5) {
        // 0.5 < x < 1, y is not negligible
        // log1p(x^2 + y^2 - 1)/2 / ln2
        this.re = Math.log1p(x2y2m1(absx, absy)) * M_1_2_LN2;
      } else {
        // general
        // log2(hypot(x, y))
        this.re = Math.log2(hypot(absx, absy));
      }
    }
  }
  return this;
}

// take the complex power of a complex number
// static pow(x,c)
// pow(c)
// powEq(c)

Complex.pow = function(x,c) {
  return Complex.log(x).mul(c).exp();
}

Complex.prototype.pow = function(c) {
  return this.log().mul(c).exp();
}

Complex.prototype.powEq = function(c) {
  return this.logEq().mulEq(c).expEq();
}

// take the principal square root of a complex number
// static sqrt(z)
// sqrt()
// sqrtEq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 81.74115%
// off by 1ULP: 18.258%
// off by 2ULP: 0.00085%
// off by 3ULP+: 0%
//
// average error: +/- 0.182597 ULP

Complex.sqrt = function(x) {
  var xre = x.re;
  var xim = x.im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);
  var re;
  var im;

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls === FP_INFINITE) {
      re = Infinity;
      im = xim;
    } else if (rcls === FP_INFINITE) {
      if (xre < 0) {
        re = icls == FP_NAN ? NaN : 0;
        im = copysign(Infinity, xim);
      } else {
        re = xre;
        im = (icls === FP_NAN
                      ? NaN : copysign(0, xim));
      }
    } else {
      re = NaN;
      im = NaN;
    }
  } else {
    if (icls === FP_ZERO) {
      if (xre < 0) {
        re = 0;
        im = copysign(Math.sqrt(-xre), xim);
      } else {
        re = fabs(Math.sqrt(xre));
        im = copysign(0, xim);
      }
    } else if (rcls === FP_ZERO) {
      var r;
      if (fabs(xim) >= 2 * M_MIN)
        r = Math.sqrt(0.5 * fabs(xim));
      else
        r = 0.5 * Math.sqrt(2 * fabs(xim));

      re = r;
      im = copysign(r, xim);
    } else {
      var d, r, s;
      var scale = 0;
      var absx = xre > 0 ? xre : -xre;
      var absy = xim > 0 ? xim : -xim;

      if (absx > M_MAX / 4) {
        scale = 1;
        xre = scalbn(xre, -2 * scale);
        xim = scalbn(xim, -2 * scale);
      } else if (absy > M_MAX / 4) {
        scale = 1;
        if (absx >= 4 * M_MIN)
          xre = scalbn(xre, -2 * scale);
        else
          xre = 0;
        xim = scalbn(xim, -2 * scale);
      } else if (absx < 2 * M_MIN
              && absy < 2 * M_MIN) {
        scale = -((M_MANT_DIG + 1) / 2);
        xre = scalbn(xre, -2 * scale);
        xim = scalbn(xim, -2 * scale);
      }

      d = hypot(xre, xim);
      /* Use the identity   2  Re res  Im res = Im x
         to avoid cancellation error in  d +/- Re x.  */
      if (xre > 0) {
        r = Math.sqrt(0.5 * d + 0.5 * xre);
        if (scale == 1 && fabs(xim) < 1) {
          /* Avoid possible intermediate underflow.  */
          s = xim / r;
          r = scalbn(r, scale);
          scale = 0;
        } else
        s = 0.5 * (xim / r);
      } else {
        s = Math.sqrt(0.5 * d - 0.5 * xre);
        if (scale == 1 && fabs(xim) < 1) {
          /* Avoid possible intermediate underflow.  */
          r = abs(xim / s);
          s = scalbn(s, scale);
          scale = 0;
        } else
        r = abs(0.5 * (xim / s));
      }

      if (scale) {
        r = scalbn(r, scale);
        s = scalbn(s, scale);
      }

      re = r;
      im = copysign(s, xim);
    }
  }

  return new Complex(re,im);
}

Complex.prototype.sqrt = function() {
  var xre = this.re;
  var xim = this.im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);
  var re;
  var im;

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls === FP_INFINITE) {
      re = Infinity;
      im = xim;
    } else if (rcls === FP_INFINITE) {
      if (xre < 0) {
        re = icls == FP_NAN ? NaN : 0;
        im = copysign(Infinity, xim);
      } else {
        re = xre;
        im = (icls === FP_NAN
                      ? NaN : copysign(0, xim));
      }
    } else {
      re = NaN;
      im = NaN;
    }
  } else {
    if (icls === FP_ZERO) {
      if (xre < 0) {
        re = 0;
        im = copysign(Math.sqrt(-xre), xim);
      } else {
        re = fabs(Math.sqrt(xre));
        im = copysign(0, xim);
      }
    } else if (rcls === FP_ZERO) {
      var r;
      if (fabs(xim) >= 2 * M_MIN)
        r = Math.sqrt(0.5 * fabs(xim));
      else
        r = 0.5 * Math.sqrt(2 * fabs(xim));

      re = r;
      im = copysign(r, xim);
    } else {
      var d, r, s;
      var scale = 0;
      var absx = xre > 0 ? xre : -xre;
      var absy = xim > 0 ? xim : -xim;

      if (absx > M_MAX / 4) {
        scale = 1;
        xre = scalbn(xre, -2 * scale);
        xim = scalbn(xim, -2 * scale);
      } else if (absy > M_MAX / 4) {
        scale = 1;
        if (absx >= 4 * M_MIN)
          xre = scalbn(xre, -2 * scale);
        else
          xre = 0;
        xim = scalbn(xim, -2 * scale);
      } else if (absx < 2 * M_MIN
              && absy < 2 * M_MIN) {
        scale = -((M_MANT_DIG + 1) / 2);
        xre = scalbn(xre, -2 * scale);
        xim = scalbn(xim, -2 * scale);
      }

      d = hypot(xre, xim);
      /* Use the identity   2  Re res  Im res = Im x
         to avoid cancellation error in  d +/- Re x.  */
      if (xre > 0) {
        r = Math.sqrt(0.5 * d + 0.5 * xre);
        if (scale == 1 && fabs(xim) < 1) {
          /* Avoid possible intermediate underflow.  */
          s = xim / r;
          r = scalbn(r, scale);
          scale = 0;
        } else
        s = 0.5 * (xim / r);
      } else {
        s = Math.sqrt(0.5 * d - 0.5 * xre);
        if (scale == 1 && fabs(xim) < 1) {
          /* Avoid possible intermediate underflow.  */
          r = abs(xim / s);
          s = scalbn(s, scale);
          scale = 0;
        } else
        r = abs(0.5 * (xim / s));
      }

      if (scale) {
        r = scalbn(r, scale);
        s = scalbn(s, scale);
      }

      re = r;
      im = copysign(s, xim);
    }
  }

  return new Complex(re,im);
}

Complex.prototype.sqrtEq = function() {
  var xre = this.re;
  var xim = this.im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);
  var re;
  var im;

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls === FP_INFINITE) {
      re = Infinity;
      im = xim;
    } else if (rcls === FP_INFINITE) {
      if (xre < 0) {
        re = icls == FP_NAN ? NaN : 0;
        im = copysign(Infinity, xim);
      } else {
        re = xre;
        im = (icls === FP_NAN
                      ? NaN : copysign(0, xim));
      }
    } else {
      re = NaN;
      im = NaN;
    }
  } else {
    if (icls === FP_ZERO) {
      if (xre < 0) {
        re = 0;
        im = copysign(Math.sqrt(-xre), xim);
      } else {
        re = fabs(Math.sqrt(xre));
        im = copysign(0, xim);
      }
    } else if (rcls === FP_ZERO) {
      var r;
      if (fabs(xim) >= 2 * M_MIN)
        r = Math.sqrt(0.5 * fabs(xim));
      else
        r = 0.5 * Math.sqrt(2 * fabs(xim));

      re = r;
      im = copysign(r, xim);
    } else {
      var d, r, s;
      var scale = 0;
      var absx = xre > 0 ? xre : -xre;
      var absy = xim > 0 ? xim : -xim;

      if (absx > M_MAX / 4) {
        scale = 1;
        xre = scalbn(xre, -2 * scale);
        xim = scalbn(xim, -2 * scale);
      } else if (absy > M_MAX / 4) {
        scale = 1;
        if (absx >= 4 * M_MIN)
          xre = scalbn(xre, -2 * scale);
        else
          xre = 0;
        xim = scalbn(xim, -2 * scale);
      } else if (absx < 2 * M_MIN
              && absy < 2 * M_MIN) {
        scale = -((M_MANT_DIG + 1) / 2);
        xre = scalbn(xre, -2 * scale);
        xim = scalbn(xim, -2 * scale);
      }

      d = hypot(xre, xim);
      /* Use the identity   2  Re res  Im res = Im x
         to avoid cancellation error in  d +/- Re x.  */
      if (xre > 0) {
        r = Math.sqrt(0.5 * d + 0.5 * xre);
        if (scale == 1 && fabs(xim) < 1) {
          /* Avoid possible intermediate underflow.  */
          s = xim / r;
          r = scalbn(r, scale);
          scale = 0;
        } else
        s = 0.5 * (xim / r);
      } else {
        s = Math.sqrt(0.5 * d - 0.5 * xre);
        if (scale == 1 && fabs(xim) < 1) {
          /* Avoid possible intermediate underflow.  */
          r = abs(xim / s);
          s = scalbn(s, scale);
          scale = 0;
        } else
        r = abs(0.5 * (xim / s));
      }

      if (scale) {
        r = scalbn(r, scale);
        s = scalbn(s, scale);
      }

      re = r;
      im = copysign(s, xim);
    }
  }

  this.re = re;
  this.im = im;
  return this;
}

// ──── TRIGONOMETRY ───────────────────────────────────────────────────────────
// take the sine of a complex number
// static sin(z)
// sin()
// sinEq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 96.00545%
// off by 1ULP: 3.98515%
// off by 2ULP: 0.0094%
// off by 3ULP+: 0%
//
// average error: +/- 0.0400395 ULP

Complex.sin = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  var negate = signbit(xre);
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  xre = fabs(xre);

  if (icls >= FP_ZERO) {
    /* Imaginary part is finite.  */
    if (rcls >= FP_ZERO) {
      /* Real part is finite.  */
      var t = 709;
      var sinix, cosix;

      if (xre > M_MIN) {
        sinix = Math.sin(xre);
        cosix = Math.cos(xre);
      } else {
        sinix = xre;
        cosix = 1;
      }

      if (negate)
        sinix = -sinix;

      if (fabs(xim) > t) {
        var exp_t = Math.exp(t);
        var ix = fabs(xim);
        if (signbit(xim))
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
          re = M_MAX * sinix;
          im = M_MAX * cosix;
        } else {
          var exp_val = Math.exp(ix);
          re = exp_val * sinix;
          im = exp_val * cosix;
        }
      } else {
        re = Math.cosh(xim) * sinix;
        im = Math.sinh(xim) * cosix;
      }
    } else {
      if (icls == FP_ZERO) {
        /* Imaginary part is 0.0.  */
        re = xre - xre;
        im = xim;
      } else {
        re = NaN;
        im = NaN;
      }
    }
  } else if (icls == FP_INFINITE) {
    /* Imaginary part is infinite.  */
    if (rcls == FP_ZERO) {
      /* Real part is 0.0.  */
      re = copysign(0, negate ? -1 : 1);
      im = xim;
    } else if (rcls > FP_ZERO) {
      /* Real part is finite.  */
      var sinix, cosix;

      if (xre > M_MIN) {
        sinix = Math.sin(xre);
        cosix = Math.cos(xre);
      } else {
        sinix = xre;
        cosix = 1;
      }

      re = copysign(Infinity, sinix);
      im = copysign(Infinity, cosix);

      if (negate)
        re = -re;
      if (signbit(xim))
        im = -im;
    } else {
      re = xre - xre;
      im = Infinity;
    }
  } else {
    if (rcls == FP_ZERO)
      re = copysign(0, negate ? -1 : 1);
    else
      re = NaN;
    im = NaN;
  }

  return new Complex(re,im);
}

Complex.prototype.sin = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var negate = signbit(xre);
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  xre = fabs(xre);

  if (icls >= FP_ZERO) {
    /* Imaginary part is finite.  */
    if (rcls >= FP_ZERO) {
      /* Real part is finite.  */
      var t = 709;
      var sinix, cosix;

      if (xre > M_MIN) {
        sinix = Math.sin(xre);
        cosix = Math.cos(xre);
      } else {
        sinix = xre;
        cosix = 1;
      }

      if (negate)
        sinix = -sinix;

      if (fabs(xim) > t) {
        var exp_t = Math.exp(t);
        var ix = fabs(xim);
        if (signbit(xim))
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
          re = M_MAX * sinix;
          im = M_MAX * cosix;
        } else {
          var exp_val = Math.exp(ix);
          re = exp_val * sinix;
          im = exp_val * cosix;
        }
      } else {
        re = Math.cosh(xim) * sinix;
        im = Math.sinh(xim) * cosix;
      }
    } else {
      if (icls == FP_ZERO) {
        /* Imaginary part is 0.0.  */
        re = xre - xre;
        im = xim;
      } else {
        re = NaN;
        im = NaN;
      }
    }
  } else if (icls == FP_INFINITE) {
    /* Imaginary part is infinite.  */
    if (rcls == FP_ZERO) {
      /* Real part is 0.0.  */
      re = copysign(0, negate ? -1 : 1);
      im = xim;
    } else if (rcls > FP_ZERO) {
      /* Real part is finite.  */
      var sinix, cosix;

      if (xre > M_MIN) {
        sinix = Math.sin(xre);
        cosix = Math.cos(xre);
      } else {
        sinix = xre;
        cosix = 1;
      }

      re = copysign(Infinity, sinix);
      im = copysign(Infinity, cosix);

      if (negate)
        re = -re;
      if (signbit(xim))
        im = -im;
    } else {
      re = xre - xre;
      im = Infinity;
    }
  } else {
    if (rcls == FP_ZERO)
      re = copysign(0, negate ? -1 : 1);
    else
      re = NaN;
    im = NaN;
  }

  return new Complex(re,im);
}

Complex.prototype.sinEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var negate = signbit(xre);
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  xre = fabs(xre);

  if (icls >= FP_ZERO) {
    /* Imaginary part is finite.  */
    if (rcls >= FP_ZERO) {
      /* Real part is finite.  */
      var t = 709;
      var sinix, cosix;

      if (xre > M_MIN) {
        sinix = Math.sin(xre);
        cosix = Math.cos(xre);
      } else {
        sinix = xre;
        cosix = 1;
      }

      if (negate)
        sinix = -sinix;

      if (fabs(xim) > t) {
        var exp_t = Math.exp(t);
        var ix = fabs(xim);
        if (signbit(xim))
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
          re = M_MAX * sinix;
          im = M_MAX * cosix;
        } else {
          var exp_val = Math.exp(ix);
          re = exp_val * sinix;
          im = exp_val * cosix;
        }
      } else {
        re = Math.cosh(xim) * sinix;
        im = Math.sinh(xim) * cosix;
      }
    } else {
      if (icls == FP_ZERO) {
        /* Imaginary part is 0.0.  */
        re = xre - xre;
        im = xim;
      } else {
        re = NaN;
        im = NaN;
      }
    }
  } else if (icls == FP_INFINITE) {
    /* Imaginary part is infinite.  */
    if (rcls == FP_ZERO) {
      /* Real part is 0.0.  */
      re = copysign(0, negate ? -1 : 1);
      im = xim;
    } else if (rcls > FP_ZERO) {
      /* Real part is finite.  */
      var sinix, cosix;

      if (xre > M_MIN) {
        sinix = Math.sin(xre);
        cosix = Math.cos(xre);
      } else {
        sinix = xre;
        cosix = 1;
      }

      re = copysign(Infinity, sinix);
      im = copysign(Infinity, cosix);

      if (negate)
        re = -re;
      if (signbit(xim))
        im = -im;
    } else {
      re = xre - xre;
      im = Infinity;
    }
  } else {
    if (rcls == FP_ZERO)
      re = copysign(0, negate ? -1 : 1);
    else
      re = NaN;
    im = NaN;
  }

  this.re = re;
  this.im = im;
  return this;
}

// take the cosine of a complex number
// static cos(z)
// cos()
// cosEq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 96.0604%
// off by 1ULP: 3.9282%
// off by 2ULP: 0.0114%
// off by 3ULP+: 0%
//
// average error: +/- 0.03951 ULP

Complex.cos = function(x) {
  var y = new Complex(-x.im, x.re);
  return y.coshEq();
}

Complex.prototype.cos = function() {
  var y = new Complex(-this.im, this.re);
  return y.coshEq();
}

Complex.prototype.cosEq = function() {
  var re = this.re;
  this.re = -this.im;
  this.im = re;
  return this.coshEq();
}

// take the tangent of a complex number
// static tan(z)
// tan()
// tanEq()

Complex.tan = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;

  if (!isfinite(xre) || !isfinite(xim)) {
    if (isinf(xim)) {
      if (isfinite(xre) && fabs(xre) > 1) {
        var sinrx, cosrx;
        sinrx = Math.sin(xre);
        cosrx = Math.cos(xre);
        re = copysign(0, sinrx * cosrx);
      } else
        re = copysign(0, xre);
      im = copysign(1, xim);
    } else if (xre == 0) {
      re = xre;
      im = xim;
    } else {
      re = NaN;
      if (xim == 0)
        im = xim;
      else
        im = NaN;
    }
  } else {
    var sinrx, cosrx;
    var den;
    var t = ((M_MAX_EXP - 1) * M_LN2 / 2) | 0;

    /* tan(x+iy) = (sin(2x) + i*sinh(2y))/(cos(2x) + cosh(2y))
       = (sin(x)*cos(x) + i*sinh(y)*cosh(y)/(cos(x)^2 + sinh(y)^2). */

    if (fabs(xre) > M_MIN) {
      sinrx = Math.sin(xre);
      cosrx = Math.cos(xre);
    } else {
      sinrx = xre;
      cosrx = 1;
    }

    if (fabs(xim) > t) {
      /* Avoid intermediate overflow when the real part of the
         result may be subnormal.  Ignoring negligible terms, the
         imaginary part is +/- 1, the real part is
         sin(x)*cos(x)/sinh(y)^2 = 4*sin(x)*cos(x)/exp(2y).  */
      var exp_2t = Math.exp(2 * t);

      im = copysign(1, xim);
      re = 4 * sinrx * cosrx;
      xim = fabs(xim);
      xim -= t;
      re /= exp_2t;
      if (xim > t) {
        /* Underflow (original imaginary part of x has absolute
           value > 2t).  */
        re /= exp_2t;
      } else
        re /= Math.exp(2 * xim);
    } else {
      var sinhix, coshix;
      if (fabs(xim) > M_MIN) {
        sinhix = Math.sinh(xim);
        coshix = Math.cosh(xim);
      } else {
        sinhix = xim;
        coshix = 1;
      }

      if (fabs(sinhix) > fabs(cosrx) * M_EPSILON)
        den = cosrx * cosrx + sinhix * sinhix;
      else
        den = cosrx * cosrx;
      re = sinrx * cosrx / den;
      im = sinhix * coshix / den;
    }
  }

  return new Complex(re,im);
}

Complex.prototype.tan = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;

  if (!isfinite(xre) || !isfinite(xim)) {
    if (isinf(xim)) {
      if (isfinite(xre) && fabs(xre) > 1) {
        var sinrx, cosrx;
        sinrx = Math.sin(xre);
        cosrx = Math.cos(xre);
        re = copysign(0, sinrx * cosrx);
      } else
        re = copysign(0, xre);
      im = copysign(1, xim);
    } else if (xre == 0) {
      re = xre;
      im = xim;
    } else {
      re = NaN;
      if (xim == 0)
        im = xim;
      else
        im = NaN;
    }
  } else {
    var sinrx, cosrx;
    var den;
    var t = ((M_MAX_EXP - 1) * M_LN2 / 2) | 0;

    /* tan(x+iy) = (sin(2x) + i*sinh(2y))/(cos(2x) + cosh(2y))
       = (sin(x)*cos(x) + i*sinh(y)*cosh(y)/(cos(x)^2 + sinh(y)^2). */

    if (fabs(xre) > M_MIN) {
      sinrx = Math.sin(xre);
      cosrx = Math.cos(xre);
    } else {
      sinrx = xre;
      cosrx = 1;
    }

    if (fabs(xim) > t) {
      /* Avoid intermediate overflow when the real part of the
         result may be subnormal.  Ignoring negligible terms, the
         imaginary part is +/- 1, the real part is
         sin(x)*cos(x)/sinh(y)^2 = 4*sin(x)*cos(x)/exp(2y).  */
      var exp_2t = Math.exp(2 * t);

      im = copysign(1, xim);
      re = 4 * sinrx * cosrx;
      xim = fabs(xim);
      xim -= t;
      re /= exp_2t;
      if (xim > t) {
        /* Underflow (original imaginary part of x has absolute
           value > 2t).  */
        re /= exp_2t;
      } else
        re /= Math.exp(2 * xim);
    } else {
      var sinhix, coshix;
      if (fabs(xim) > M_MIN) {
        sinhix = Math.sinh(xim);
        coshix = Math.cosh(xim);
      } else {
        sinhix = xim;
        coshix = 1;
      }

      if (fabs(sinhix) > fabs(cosrx) * M_EPSILON)
        den = cosrx * cosrx + sinhix * sinhix;
      else
        den = cosrx * cosrx;
      re = sinrx * cosrx / den;
      im = sinhix * coshix / den;
    }
  }

  return new Complex(re,im);
}

Complex.prototype.tanEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;

  if (!isfinite(xre) || !isfinite(xim)) {
    if (isinf(xim)) {
      if (isfinite(xre) && fabs(xre) > 1) {
        var sinrx, cosrx;
        sinrx = Math.sin(xre);
        cosrx = Math.cos(xre);
        re = copysign(0, sinrx * cosrx);
      } else
        re = copysign(0, xre);
      im = copysign(1, xim);
    } else if (xre == 0) {
      re = xre;
      im = xim;
    } else {
      re = NaN;
      if (xim == 0)
        im = xim;
      else
        im = NaN;
    }
  } else {
    var sinrx, cosrx;
    var den;
    var t = ((M_MAX_EXP - 1) * M_LN2 / 2) | 0;

    /* tan(x+iy) = (sin(2x) + i*sinh(2y))/(cos(2x) + cosh(2y))
       = (sin(x)*cos(x) + i*sinh(y)*cosh(y)/(cos(x)^2 + sinh(y)^2). */

    if (fabs(xre) > M_MIN) {
      sinrx = Math.sin(xre);
      cosrx = Math.cos(xre);
    } else {
      sinrx = xre;
      cosrx = 1;
    }

    if (fabs(xim) > t) {
      /* Avoid intermediate overflow when the real part of the
         result may be subnormal.  Ignoring negligible terms, the
         imaginary part is +/- 1, the real part is
         sin(x)*cos(x)/sinh(y)^2 = 4*sin(x)*cos(x)/exp(2y).  */
      var exp_2t = Math.exp(2 * t);

      im = copysign(1, xim);
      re = 4 * sinrx * cosrx;
      xim = fabs(xim);
      xim -= t;
      re /= exp_2t;
      if (xim > t) {
        /* Underflow (original imaginary part of x has absolute
           value > 2t).  */
        re /= exp_2t;
      } else
        re /= Math.exp(2 * xim);
    } else {
      var sinhix, coshix;
      if (fabs(xim) > M_MIN) {
        sinhix = Math.sinh(xim);
        coshix = Math.cosh(xim);
      } else {
        sinhix = xim;
        coshix = 1;
      }

      if (fabs(sinhix) > fabs(cosrx) * M_EPSILON)
        den = cosrx * cosrx + sinhix * sinhix;
      else
        den = cosrx * cosrx;
      re = sinrx * cosrx / den;
      im = sinhix * coshix / den;
    }
  }

  this.re = re;
  this.im = im;
  return this;
}

// take the hyperbolic sine of a complex number
// static sinh(z)
// sinh()
// sinhEq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 96.00425%
// off by 1ULP: 3.9865%
// off by 2ULP: 0.00925%
// off by 3ULP+: 0%
//
// average error: +/- 0.04005 ULP

Complex.sinh = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  var negate = signbit(xre);
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  xre = fabs(xre);

  if (rcls >= FP_ZERO) {
    /* Real part is finite.  */
    if (icls >= FP_ZERO) {
      /* Imaginary part is finite.  */
      var t = 709;
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      if (negate)
        cosix = -cosix;

      if (fabs(xre) > t) {
        var exp_t = Math.exp(t);
        var rx = fabs(xre);
        if (signbit(xre))
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
          re = M_MAX * cosix;
          im = M_MAX * sinix;
        } else {
          var exp_val = Math.exp(rx);
          re = exp_val * cosix;
          im = exp_val * sinix;
        }
      } else {
        re = Math.sinh(xre) * cosix;
        im = Math.cosh(xre) * sinix;
      }
    } else {
      if (rcls == FP_ZERO) {
        /* Real part is 0.0.  */
        re = copysign(0, negate ? -1 : 1);
        im = xim - xim;
      } else {
        re = NaN;
        im = NaN;
      }
    }
  } else if (rcls == FP_INFINITE) {
    /* Real part is infinite.  */
    if (icls > FP_ZERO) {
      /* Imaginary part is finite.  */
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      re = copysign(Infinity, cosix);
      im = copysign(Infinity, sinix);

      if (negate)
        re = -re;
    } else if (icls == FP_ZERO) {
      /* Imaginary part is 0.0.  */
      re = negate ? -Infinity : Infinity;
      im = xim;
    } else {
      re = Infinity;
      im = xim - xim;
    }
  } else {
    re = NaN;
    im = xim == 0 ? xim : NaN;
  }

  return new Complex(re,im);
}

Complex.prototype.sinh = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var negate = signbit(xre);
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  xre = fabs(xre);

  if (rcls >= FP_ZERO) {
    /* Real part is finite.  */
    if (icls >= FP_ZERO) {
      /* Imaginary part is finite.  */
      var t = 709;
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      if (negate)
        cosix = -cosix;

      if (fabs(xre) > t) {
        var exp_t = Math.exp(t);
        var rx = fabs(xre);
        if (signbit(xre))
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
          re = M_MAX * cosix;
          im = M_MAX * sinix;
        } else {
          var exp_val = Math.exp(rx);
          re = exp_val * cosix;
          im = exp_val * sinix;
        }
      } else {
        re = Math.sinh(xre) * cosix;
        im = Math.cosh(xre) * sinix;
      }
    } else {
      if (rcls == FP_ZERO) {
        /* Real part is 0.0.  */
        re = copysign(0, negate ? -1 : 1);
        im = xim - xim;
      } else {
        re = NaN;
        im = NaN;
      }
    }
  } else if (rcls == FP_INFINITE) {
    /* Real part is infinite.  */
    if (icls > FP_ZERO) {
      /* Imaginary part is finite.  */
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      re = copysign(Infinity, cosix);
      im = copysign(Infinity, sinix);

      if (negate)
        re = -re;
    } else if (icls == FP_ZERO) {
      /* Imaginary part is 0.0.  */
      re = negate ? -Infinity : Infinity;
      im = xim;
    } else {
      re = Infinity;
      im = xim - xim;
    }
  } else {
    re = NaN;
    im = xim == 0 ? xim : NaN;
  }

  return new Complex(re,im);
}

Complex.prototype.sinhEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var negate = signbit(xre);
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  xre = fabs(xre);

  if (rcls >= FP_ZERO) {
    /* Real part is finite.  */
    if (icls >= FP_ZERO) {
      /* Imaginary part is finite.  */
      var t = 709;
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      if (negate)
        cosix = -cosix;

      if (fabs(xre) > t) {
        var exp_t = Math.exp(t);
        var rx = fabs(xre);
        if (signbit(xre))
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
          re = M_MAX * cosix;
          im = M_MAX * sinix;
        } else {
          var exp_val = Math.exp(rx);
          re = exp_val * cosix;
          im = exp_val * sinix;
        }
      } else {
        re = Math.sinh(xre) * cosix;
        im = Math.cosh(xre) * sinix;
      }
    } else {
      if (rcls == FP_ZERO) {
        /* Real part is 0.0.  */
        re = copysign(0, negate ? -1 : 1);
        im = xim - xim;
      } else {
        re = NaN;
        im = NaN;
      }
    }
  } else if (rcls == FP_INFINITE) {
    /* Real part is infinite.  */
    if (icls > FP_ZERO) {
      /* Imaginary part is finite.  */
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      re = copysign(Infinity, cosix);
      im = copysign(Infinity, sinix);

      if (negate)
        re = -re;
    } else if (icls == FP_ZERO) {
      /* Imaginary part is 0.0.  */
      re = negate ? -Infinity : Infinity;
      im = xim;
    } else {
      re = Infinity;
      im = xim - xim;
    }
  } else {
    re = NaN;
    im = xim == 0 ? xim : NaN;
  }

  this.re = re;
  this.im = im;
  return this;
}

// take the hyperbolic cosine of a complex number
// static cosh(z)
// cosh()
// coshEq()
//
// out of 1_000_000 random inputs:
// off by 0ULP: 96.0911%
// off by 1ULP: 3.8963%
// off by 2ULP: 0.0126%
// off by 3ULP+: 0%
//
// average error: +/- 0.039215 ULP

Complex.cosh = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls >= FP_ZERO) {
    /* Real part is finite.  */
    if (icls >= FP_ZERO) {
      /* Imaginary part is finite.  */
      var t = 709;
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      if (fabs(xre) > t) {
        var exp_t = Math.exp(t);
        var rx = fabs(xre);
        if (signbit(xre))
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
          re = M_MAX * cosix;
          im = M_MAX * sinix;
        } else {
          var exp_val = Math.exp(rx);
          re = exp_val * cosix;
          im = exp_val * sinix;
        }
      } else {
        re = Math.cosh(xre) * cosix;
        im = Math.sinh(xre) * sinix;
      }
    } else {
      im = xre == 0 ? 0 : NaN;
      re = xim - xim;
    }
  } else if (rcls == FP_INFINITE) {
    /* Real part is infinite.  */
    if (icls > FP_ZERO) {
      /* Imaginary part is finite.  */
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      re = copysign(Infinity, cosix);
      im = (copysign(Infinity, sinix)
                         * copysign(1, xre));
    } else if (icls == FP_ZERO) {
      /* Imaginary part is 0.0.  */
      re = Infinity;
      im = xim * copysign(1, xre);
    } else {
      re = Infinity;
      im = xim - xim;
    }
  } else {
    re = NaN;
    im = xim == 0 ? xim : NaN;
  }

  return new Complex(re,im);
}

Complex.prototype.cosh = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls >= FP_ZERO) {
    /* Real part is finite.  */
    if (icls >= FP_ZERO) {
      /* Imaginary part is finite.  */
      var t = 709;
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      if (fabs(xre) > t) {
        var exp_t = Math.exp(t);
        var rx = fabs(xre);
        if (signbit(xre))
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
          re = M_MAX * cosix;
          im = M_MAX * sinix;
        } else {
          var exp_val = Math.exp(rx);
          re = exp_val * cosix;
          im = exp_val * sinix;
        }
      } else {
        re = Math.cosh(xre) * cosix;
        im = Math.sinh(xre) * sinix;
      }
    } else {
      im = xre == 0 ? 0 : NaN;
      re = xim - xim;
    }
  } else if (rcls == FP_INFINITE) {
    /* Real part is infinite.  */
    if (icls > FP_ZERO) {
      /* Imaginary part is finite.  */
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      re = copysign(Infinity, cosix);
      im = (copysign(Infinity, sinix)
                         * copysign(1, xre));
    } else if (icls == FP_ZERO) {
      /* Imaginary part is 0.0.  */
      re = Infinity;
      im = xim * copysign(1, xre);
    } else {
      re = Infinity;
      im = xim - xim;
    }
  } else {
    re = NaN;
    im = xim == 0 ? xim : NaN;
  }

  return new Complex(re,im);
}

Complex.prototype.coshEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls >= FP_ZERO) {
    /* Real part is finite.  */
    if (icls >= FP_ZERO) {
      /* Imaginary part is finite.  */
      var t = 709;
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      if (fabs(xre) > t) {
        var exp_t = Math.exp(t);
        var rx = fabs(xre);
        if (signbit(xre))
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
          re = M_MAX * cosix;
          im = M_MAX * sinix;
        } else {
          var exp_val = Math.exp(rx);
          re = exp_val * cosix;
          im = exp_val * sinix;
        }
      } else {
        re = Math.cosh(xre) * cosix;
        im = Math.sinh(xre) * sinix;
      }
    } else {
      im = xre == 0 ? 0 : NaN;
      re = xim - xim;
    }
  } else if (rcls == FP_INFINITE) {
    /* Real part is infinite.  */
    if (icls > FP_ZERO) {
      /* Imaginary part is finite.  */
      var sinix, cosix;

      if (fabs(xim) > M_MIN) {
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
      } else {
        sinix = xim;
        cosix = 1;
      }

      re = copysign(Infinity, cosix);
      im = (copysign(Infinity, sinix)
                         * copysign(1, xre));
    } else if (icls == FP_ZERO) {
      /* Imaginary part is 0.0.  */
      re = Infinity;
      im = xim * copysign(1, xre);
    } else {
      re = Infinity;
      im = xim - xim;
    }
  } else {
    re = NaN;
    im = xim == 0 ? xim : NaN;
  }

  this.re = re;
  this.im = im;
  return this;
}

// take the hyperbolic tangent of a complex number
// static tanh(z)
// tanh()
// atanhEq()

Complex.tanh = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;

  if (!isfinite(xre) || !isfinite(xim)) {
    if (isinf(xre)) {
      re = copysign(1, xre);
      if (isfinite(xim) && fabs(xim) > 1) {
        var sinix, cosix;
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
        im = copysign(0, sinix * cosix);
      } else
        im = copysign(0, xim);
    } else if (xim == 0) {
      re = xre;
      im = xim;
    } else {
      if (xre == 0)
        re = xre;
      else
        re = NaN;
      im = NaN;
    }
  } else {
    var sinix, cosix;
    var den;
    var t = ((M_MAX_EXP - 1) * M_LN2 / 2) | 0;

    /* tanh(x+iy) = (sinh(2x) + i*sin(2y))/(cosh(2x) + cos(2y))
       = (sinh(x)*cosh(x) + i*sin(y)*cos(y))/(sinh(x)^2 + cos(y)^2).  */

    if (fabs(xim) > M_MIN) {
      sinix = Math.sin(xim);
      cosix = Math.cos(xim);
    } else {
      sinix = xim;
      cosix = 1;
    }

    if (fabs(xre) > t) {
      /* Avoid intermediate overflow when the imaginary part of
         the result may be subnormal.  Ignoring negligible terms,
         the real part is +/- 1, the imaginary part is
         sin(y)*cos(y)/sinh(x)^2 = 4*sin(y)*cos(y)/exp(2x).  */
      var exp_2t = Math.exp(2 * t);

      re = copysign(1, xre);
      im = 4 * sinix * cosix;
      xre = fabs(xre);
      xre -= t;
      im /= exp_2t;
      if (xre > t) {
        /* Underflow (original real part of x has absolute value
           > 2t).  */
        im /= exp_2t;
      } else
        im /= Math.exp(2 * xre);
    } else {
      var sinhrx, coshrx;
      if (fabs(xre) > M_MIN) {
        sinhrx = Math.sinh(xre);
        coshrx = Math.cosh(xre);
      } else {
        sinhrx = xre;
        coshrx = 1;
      }

      if (fabs(sinhrx) > fabs(cosix) * M_EPSILON)
        den = sinhrx * sinhrx + cosix * cosix;
      else
        den = cosix * cosix;
      re = sinhrx * coshrx / den;
      im = sinix * cosix / den;
    }
  }

  return new Complex(re,im);
}

Complex.prototype.tanh = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;

  if (!isfinite(xre) || !isfinite(xim)) {
    if (isinf(xre)) {
      re = copysign(1, xre);
      if (isfinite(xim) && fabs(xim) > 1) {
        var sinix, cosix;
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
        im = copysign(0, sinix * cosix);
      } else
        im = copysign(0, xim);
    } else if (xim == 0) {
      re = xre;
      im = xim;
    } else {
      if (xre == 0)
        re = xre;
      else
        re = NaN;
      im = NaN;
    }
  } else {
    var sinix, cosix;
    var den;
    var t = ((M_MAX_EXP - 1) * M_LN2 / 2) | 0;

    /* tanh(x+iy) = (sinh(2x) + i*sin(2y))/(cosh(2x) + cos(2y))
       = (sinh(x)*cosh(x) + i*sin(y)*cos(y))/(sinh(x)^2 + cos(y)^2).  */

    if (fabs(xim) > M_MIN) {
      sinix = Math.sin(xim);
      cosix = Math.cos(xim);
    } else {
      sinix = xim;
      cosix = 1;
    }

    if (fabs(xre) > t) {
      /* Avoid intermediate overflow when the imaginary part of
         the result may be subnormal.  Ignoring negligible terms,
         the real part is +/- 1, the imaginary part is
         sin(y)*cos(y)/sinh(x)^2 = 4*sin(y)*cos(y)/exp(2x).  */
      var exp_2t = Math.exp(2 * t);

      re = copysign(1, xre);
      im = 4 * sinix * cosix;
      xre = fabs(xre);
      xre -= t;
      im /= exp_2t;
      if (xre > t) {
        /* Underflow (original real part of x has absolute value
           > 2t).  */
        im /= exp_2t;
      } else
        im /= Math.exp(2 * xre);
    } else {
      var sinhrx, coshrx;
      if (fabs(xre) > M_MIN) {
        sinhrx = Math.sinh(xre);
        coshrx = Math.cosh(xre);
      } else {
        sinhrx = xre;
        coshrx = 1;
      }

      if (fabs(sinhrx) > fabs(cosix) * M_EPSILON)
        den = sinhrx * sinhrx + cosix * cosix;
      else
        den = cosix * cosix;
      re = sinhrx * coshrx / den;
      im = sinix * cosix / den;
    }
  }

  return new Complex(re,im);
}

Complex.prototype.tanhEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;

  if (!isfinite(xre) || !isfinite(xim)) {
    if (isinf(xre)) {
      re = copysign(1, xre);
      if (isfinite(xim) && fabs(xim) > 1) {
        var sinix, cosix;
        sinix = Math.sin(xim);
        cosix = Math.cos(xim);
        im = copysign(0, sinix * cosix);
      } else
        im = copysign(0, xim);
    } else if (xim == 0) {
      re = xre;
      im = xim;
    } else {
      if (xre == 0)
        re = xre;
      else
        re = NaN;
      im = NaN;
    }
  } else {
    var sinix, cosix;
    var den;
    var t = ((M_MAX_EXP - 1) * M_LN2 / 2) | 0;

    /* tanh(x+iy) = (sinh(2x) + i*sin(2y))/(cosh(2x) + cos(2y))
       = (sinh(x)*cosh(x) + i*sin(y)*cos(y))/(sinh(x)^2 + cos(y)^2).  */

    if (fabs(xim) > M_MIN) {
      sinix = Math.sin(xim);
      cosix = Math.cos(xim);
    } else {
      sinix = xim;
      cosix = 1;
    }

    if (fabs(xre) > t) {
      /* Avoid intermediate overflow when the imaginary part of
         the result may be subnormal.  Ignoring negligible terms,
         the real part is +/- 1, the imaginary part is
         sin(y)*cos(y)/sinh(x)^2 = 4*sin(y)*cos(y)/exp(2x).  */
      var exp_2t = Math.exp(2 * t);

      re = copysign(1, xre);
      im = 4 * sinix * cosix;
      xre = fabs(xre);
      xre -= t;
      im /= exp_2t;
      if (xre > t) {
        /* Underflow (original real part of x has absolute value
           > 2t).  */
        im /= exp_2t;
      } else
        im /= Math.exp(2 * xre);
    } else {
      var sinhrx, coshrx;
      if (fabs(xre) > M_MIN) {
        sinhrx = Math.sinh(xre);
        coshrx = Math.cosh(xre);
      } else {
        sinhrx = xre;
        coshrx = 1;
      }

      if (fabs(sinhrx) > fabs(cosix) * M_EPSILON)
        den = sinhrx * sinhrx + cosix * cosix;
      else
        den = cosix * cosix;
      re = sinhrx * coshrx / den;
      im = sinix * cosix / den;
    }
  }

  this.re = re;
  this.im = im;
  return this;
}

// ──── INVERSE TRIGONOMETRY ───────────────────────────────────────────────────
// sin(z) = -i*sinh(i*z)
// cos(z) =    cosh(i*z)
// tan(z) = -i*tanh(i*z)
// asin(z) = -i*asinh(i*z)
// atan(z) = -i*atanh(i*z)

// take the inverse sine of a complex number

Complex.asin = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;

  if (xre !== xre || xim !== xim) {
    if (xre == 0) {
      re = xre;
      im = xim;
    } else if (isinf(xre) || isinf(xim)) {
      re = NaN;
      im = copysign(Infinity, xim);
    } else {
      re = NaN;
      im = NaN;
    }
  } else {
    var y = new Complex(-xim,xre);

    y.asinhEq();

    re = y.im;
    im = -y.re;
  }

  return new Complex(re,im);
}

Complex.prototype.asin = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;

  if (xre !== xre || xim !== xim) {
    if (xre == 0) {
      re = xre;
      im = xim;
    } else if (isinf(xre) || isinf(xim)) {
      re = NaN;
      im = copysign(Infinity, xim);
    } else {
      re = NaN;
      im = NaN;
    }
  } else {
    var y = new Complex(-xim,xre);

    y.asinhEq();

    re = y.im;
    im = -y.re;
  }

  return new Complex(re,im);
}

Complex.prototype.asinEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;

  if (xre !== xre || xim !== xim) {
    if (xre == 0) {
      re = xre;
      im = xim;
    } else if (isinf(xre) || isinf(xim)) {
      re = NaN;
      im = copysign(Infinity, xim);
    } else {
      re = NaN;
      im = NaN;
    }
  } else {
    var y = new Complex(-xim,xre);

    y.asinhEq();

    re = y.im;
    im = -y.re;
  }

  this.re = re;
  this.im = im;
  return this;
}

// take the inverse cosine of a complex number

Complex.acos = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE
    || (rcls == FP_ZERO && icls == FP_ZERO)) {
    var y = x.asin();

    re = M_PI_2 - y.re;
    if (re == 0)
      re = 0;
    im = -y.im;
  } else {
    var y = kernel_casinh(new Complex(-xim,xre), 1);

    re = y.im;
    im = y.re;
  }

  return new Complex(re,im);
}

Complex.prototype.acos = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE
    || (rcls == FP_ZERO && icls == FP_ZERO)) {
    var y = this.asin();

    re = M_PI_2 - y.re;
    if (re == 0)
      re = 0;
    im = -y.im;
  } else {
    var y = kernel_casinh(new Complex(-xim,xre), 1);

    re = y.im;
    im = y.re;
  }

  return new Complex(re,im);
}

Complex.prototype.acosEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE
    || (rcls == FP_ZERO && icls == FP_ZERO)) {
    var y = this.asin();

    re = M_PI_2 - y.re;
    if (re == 0)
      re = 0;
    im = -y.im;
  } else {
    var y = kernel_casinh(new Complex(-xim,xre), 1);

    re = y.im;
    im = y.re;
  }

  this.re = re;
  this.im = im;
  return this;
}

// take the inverse tangent of a complex number

Complex.atan = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (rcls == FP_INFINITE) {
        re = copysign(M_PI_2, xre);
        im = copysign(0, xim);
      } else if (icls == FP_INFINITE) {
        if (rcls >= FP_ZERO)
          re = copysign(M_PI_2, xre);
        else
          re = NaN;
        im = copysign(0, xim);
      } else if (icls == FP_ZERO || icls == FP_INFINITE) {
        re = NaN;
        im = copysign(0, xim);
      } else {
        re = NaN;
        im = NaN;
      }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = xre;
    im = xim;
  } else {
    if (fabs(xre) >= M_16_EPSILON
      || fabs(xim) >= M_16_EPSILON) {
      re = copysign(M_PI_2, xre);
      if (fabs(xre) <= 1)
        im = 1 / xim;
      else if (fabs(xim) <= 1)
        im = xim / xre / xre;
      else {
          var h = hypot(xre / 2, xim / 2);
          im = xim / h / h / 4;
        }
    } else {
      var den, absx, absy;

      absx = fabs(xre);
      absy = fabs(xim);
      if (absx < absy) {
          var t = absx;
          absx = absy;
          absy = t;
        }

      if (absy < M_EPSILON_2) {
          den = (1 - absx) * (1 + absx);
          if (den == 0)
            den = 0;
        } else if (absx >= 1)
        den = (1 - absx) * (1 + absx) - absy * absy;
      else if (absx >= 0.75 || absy >= 0.5)
        den = -x2y2m1(absx, absy);
      else
        den = (1 - absx) * (1 + absx) - absy * absy;

      re = 0.5 * Math.atan2(2 * xre, den);

      if (fabs(xim) == 1
          && fabs(xre) < M_EPSILON_SQ)
        im = (copysign(0.5, xim)
                      * (M_LN2 - Math.log(fabs(xre))));
      else {
        var r2 = 0, num, f;

        if (fabs(xre) >= M_EPSILON_SQ)
          r2 = xre * xre;

        num = xim + 1;
        num = r2 + num * num;

        den = xim - 1;
        den = r2 + den * den;

        f = num / den;
        if (f < 0.5)
          im = 0.25 * Math.log(f);
        else {
          num = 4 * xim;
          im = 0.25 * Math.log1p(num / den);
        }
      }
    }
  }

  return new Complex(re,im);
}

Complex.prototype.atan = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (rcls == FP_INFINITE) {
        re = copysign(M_PI_2, xre);
        im = copysign(0, xim);
      } else if (icls == FP_INFINITE) {
        if (rcls >= FP_ZERO)
          re = copysign(M_PI_2, xre);
        else
          re = NaN;
        im = copysign(0, xim);
      } else if (icls == FP_ZERO || icls == FP_INFINITE) {
        re = NaN;
        im = copysign(0, xim);
      } else {
        re = NaN;
        im = NaN;
      }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = xre;
    im = xim;
  } else {
    if (fabs(xre) >= M_16_EPSILON
      || fabs(xim) >= M_16_EPSILON) {
      re = copysign(M_PI_2, xre);
      if (fabs(xre) <= 1)
        im = 1 / xim;
      else if (fabs(xim) <= 1)
        im = xim / xre / xre;
      else {
          var h = hypot(xre / 2, xim / 2);
          im = xim / h / h / 4;
        }
    } else {
      var den, absx, absy;

      absx = fabs(xre);
      absy = fabs(xim);
      if (absx < absy) {
          var t = absx;
          absx = absy;
          absy = t;
        }

      if (absy < M_EPSILON_2) {
          den = (1 - absx) * (1 + absx);
          if (den == 0)
            den = 0;
        } else if (absx >= 1)
        den = (1 - absx) * (1 + absx) - absy * absy;
      else if (absx >= 0.75 || absy >= 0.5)
        den = -x2y2m1(absx, absy);
      else
        den = (1 - absx) * (1 + absx) - absy * absy;

      re = 0.5 * Math.atan2(2 * xre, den);

      if (fabs(xim) == 1
          && fabs(xre) < M_EPSILON_SQ)
        im = (copysign(0.5, xim)
                      * (M_LN2 - Math.log(fabs(xre))));
      else {
        var r2 = 0, num, f;

        if (fabs(xre) >= M_EPSILON_SQ)
          r2 = xre * xre;

        num = xim + 1;
        num = r2 + num * num;

        den = xim - 1;
        den = r2 + den * den;

        f = num / den;
        if (f < 0.5)
          im = 0.25 * Math.log(f);
        else {
          num = 4 * xim;
          im = 0.25 * Math.log1p(num / den);
        }
      }
    }
  }

  return new Complex(re,im);
}

Complex.prototype.atanEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (rcls == FP_INFINITE) {
        re = copysign(M_PI_2, xre);
        im = copysign(0, xim);
      } else if (icls == FP_INFINITE) {
        if (rcls >= FP_ZERO)
          re = copysign(M_PI_2, xre);
        else
          re = NaN;
        im = copysign(0, xim);
      } else if (icls == FP_ZERO || icls == FP_INFINITE) {
        re = NaN;
        im = copysign(0, xim);
      } else {
        re = NaN;
        im = NaN;
      }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = xre;
    im = xim;
  } else {
    if (fabs(xre) >= M_16_EPSILON
      || fabs(xim) >= M_16_EPSILON) {
      re = copysign(M_PI_2, xre);
      if (fabs(xre) <= 1)
        im = 1 / xim;
      else if (fabs(xim) <= 1)
        im = xim / xre / xre;
      else {
          var h = hypot(xre / 2, xim / 2);
          im = xim / h / h / 4;
        }
    } else {
      var den, absx, absy;

      absx = fabs(xre);
      absy = fabs(xim);
      if (absx < absy) {
          var t = absx;
          absx = absy;
          absy = t;
        }

      if (absy < M_EPSILON_2) {
          den = (1 - absx) * (1 + absx);
          if (den == 0)
            den = 0;
        } else if (absx >= 1)
        den = (1 - absx) * (1 + absx) - absy * absy;
      else if (absx >= 0.75 || absy >= 0.5)
        den = -x2y2m1(absx, absy);
      else
        den = (1 - absx) * (1 + absx) - absy * absy;

      re = 0.5 * Math.atan2(2 * xre, den);

      if (fabs(xim) == 1
          && fabs(xre) < M_EPSILON_SQ)
        im = (copysign(0.5, xim)
                      * (M_LN2 - Math.log(fabs(xre))));
      else {
        var r2 = 0, num, f;

        if (fabs(xre) >= M_EPSILON_SQ)
          r2 = xre * xre;

        num = xim + 1;
        num = r2 + num * num;

        den = xim - 1;
        den = r2 + den * den;

        f = num / den;
        if (f < 0.5)
          im = 0.25 * Math.log(f);
        else {
          num = 4 * xim;
          im = 0.25 * Math.log1p(num / den);
        }
      }
    }
  }

  this.re = re;
  this.im = im;
  return this;
}

// take the inverse hyperbolic sine of a complex number

Complex.asinh = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      re = copysign(Infinity, xre);

      if (rcls == FP_NAN)
        im = NaN;
      else
        im = copysign((rcls >= FP_ZERO
                                    ? M_PI_2 : M_PI_4),
                                    im);
    } else if (rcls <= FP_INFINITE) {
      re = xre;
      if ((rcls == FP_INFINITE && icls >= FP_ZERO)
          || (rcls == FP_NAN && icls == FP_ZERO))
        im = copysign(0, xim);
      else
        im = NaN;
    } else {
      re = NaN;
      im = NaN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = xre;
    im = xim;
  } else {
    var y = kernel_casinh(x, 0);

    re = y.re;
    im = y.im;
  }

  return new Complex(re,im);
}

Complex.prototype.asinh = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      re = copysign(Infinity, xre);

      if (rcls == FP_NAN)
        im = NaN;
      else
        im = copysign((rcls >= FP_ZERO
                                    ? M_PI_2 : M_PI_4),
                                    im);
    } else if (rcls <= FP_INFINITE) {
      re = xre;
      if ((rcls == FP_INFINITE && icls >= FP_ZERO)
          || (rcls == FP_NAN && icls == FP_ZERO))
        im = copysign(0, xim);
      else
        im = NaN;
    } else {
      re = NaN;
      im = NaN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = xre;
    im = xim;
  } else {
    var y = kernel_casinh(this, 0);

    re = y.re;
    im = y.im;
  }

  return new Complex(re,im);
}

Complex.prototype.asinhEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      re = copysign(Infinity, xre);

      if (rcls == FP_NAN)
        im = NaN;
      else
        im = copysign((rcls >= FP_ZERO
                                    ? M_PI_2 : M_PI_4),
                                    im);
    } else if (rcls <= FP_INFINITE) {
      re = xre;
      if ((rcls == FP_INFINITE && icls >= FP_ZERO)
          || (rcls == FP_NAN && icls == FP_ZERO))
        im = copysign(0, xim);
      else
        im = NaN;
    } else {
      re = NaN;
      im = NaN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = xre;
    im = xim;
  } else {
    var y = kernel_casinh(this, 0);

    re = y.re;
    im = y.im;
  }

  this.re = re;
  this.im = im;
  return this;
}

// take the inverse hyperbolic cosine of a complex number

Complex.acosh = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      re = Infinity;

      if (rcls == FP_NAN)
        im = NaN;
      else
        im = copysign((rcls == FP_INFINITE
                                    ? (xre < 0
                                       ? M_PI - M_PI_4
                                       : M_PI_4)
                                    : M_PI_2), xim);
    } else if (rcls == FP_INFINITE) {
      re = Infinity;

      if (icls >= FP_ZERO)
        im = copysign(signbit(xre)
                                   ? M_PI : 0, xim);
      else
        im = NaN;
    } else {
      re = NaN;
      if (rcls == FP_ZERO)
        im = M_PI_2;
      else
        im = NaN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = 0;
    im = copysign(M_PI_2, xim);
  } else {
    var y = kernel_casinh(new Complex(-xim,xre), 1);

    if (signbit(xim)) {
      re = y.re;
      im = -y.im;
    } else {
      re = -y.re;
      im = y.im;
    }
  }

  return new Complex(re,im);
}

Complex.prototype.acosh = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      re = Infinity;

      if (rcls == FP_NAN)
        im = NaN;
      else
        im = copysign((rcls == FP_INFINITE
                                    ? (xre < 0
                                       ? M_PI - M_PI_4
                                       : M_PI_4)
                                    : M_PI_2), xim);
    } else if (rcls == FP_INFINITE) {
      re = Infinity;

      if (icls >= FP_ZERO)
        im = copysign(signbit(xre)
                                   ? M_PI : 0, xim);
      else
        im = NaN;
    } else {
      re = NaN;
      if (rcls == FP_ZERO)
        im = M_PI_2;
      else
        im = NaN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = 0;
    im = copysign(M_PI_2, xim);
  } else {
    var y = kernel_casinh(new Complex(-xim,xre), 1);

    if (signbit(xim)) {
      re = y.re;
      im = -y.im;
    } else {
      re = -y.re;
      im = y.im;
    }
  }

  return new Complex(re,im);
}

Complex.prototype.acoshEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      re = Infinity;

      if (rcls == FP_NAN)
        im = NaN;
      else
        im = copysign((rcls == FP_INFINITE
                                    ? (xre < 0
                                       ? M_PI - M_PI_4
                                       : M_PI_4)
                                    : M_PI_2), xim);
    } else if (rcls == FP_INFINITE) {
      re = Infinity;

      if (icls >= FP_ZERO)
        im = copysign(signbit(xre)
                                   ? M_PI : 0, xim);
      else
        im = NaN;
    } else {
      re = NaN;
      if (rcls == FP_ZERO)
        im = M_PI_2;
      else
        im = NaN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = 0;
    im = copysign(M_PI_2, xim);
  } else {
    var y = kernel_casinh(new Complex(-xim,xre), 1);

    if (signbit(xim)) {
      re = y.re;
      im = -y.im;
    } else {
      re = -y.re;
      im = y.im;
    }
  }

  this.re = re;
  this.im = im;
  return this;
}

// take the inverse hyperbolic tangent of a complex number

Complex.atanh = function(x) {
  var xre = x.re;
  var xim = x.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      re = copysign(0, xre);
      im = copysign(M_PI_2, xim);
    } else if (rcls == FP_INFINITE || rcls == FP_ZERO) {
      re = copysign(0, xre);
      if (icls >= FP_ZERO)
        im = copysign(M_PI_2, xim);
      else
        im = NaN;
    } else {
      re = NaN;
      im = NaN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = xre;
    im = xim;
  } else {
    if (fabs(xre) >= M_16_EPSILON
        || fabs(xim) >= M_16_EPSILON) {
      im = copysign(M_PI_2, xim);
      if (fabs(xim) <= 1)
        re = 1 / xre;
      else if (fabs(xre) <= 1)
        re = xre / xim / xim;
      else {
        var h = hypot(xre / 2, xim / 2);
        re = xre / h / h / 4;
      }
    } else {
      if (fabs(xre) == 1
          && fabs(xim) < M_EPSILON_SQ)
        re = (copysign(0.5, xre)
                      * (M_LN2 - Math.log(fabs(xim))));
      else {
        var i2 = 0;
        if (fabs(xim) >= M_EPSILON_SQ)
          i2 = xim * xim;

        var num = 1 + xre;
        num = i2 + num * num;

        var den = 1 - xre;
        den = i2 + den * den;

        var f = num / den;
        if (f < 0.5)
          re = 0.25 * Math.log(f);
        else {
          num = 4 * xre;
          re = 0.25 * Math.log1p(num / den);
        }
      }

      var absx, absy, den;

      absx = fabs(xre);
      absy = fabs(xim);
      if (absx < absy) {
        var t = absx;
        absx = absy;
        absy = t;
      }

      if (absy < M_EPSILON_2) {
        den = (1 - absx) * (1 + absx);
        if (den == 0)
          den = 0;
      } else if (absx >= 1)
      den = (1 - absx) * (1 + absx) - absy * absy;
      else if (absx >= 0.75 || absy >= 0.5)
        den = -x2y2m1(absx, absy);
      else
        den = (1 - absx) * (1 + absx) - absy * absy;

      im = 0.5 * Math.atan2(2 * xim, den);
    }
  }

  return new Complex(re,im);
}

Complex.prototype.atanh = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      re = copysign(0, xre);
      im = copysign(M_PI_2, xim);
    } else if (rcls == FP_INFINITE || rcls == FP_ZERO) {
      re = copysign(0, xre);
      if (icls >= FP_ZERO)
        im = copysign(M_PI_2, xim);
      else
        im = NaN;
    } else {
      re = NaN;
      im = NaN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = xre;
    im = xim;
  } else {
    if (fabs(xre) >= M_16_EPSILON
        || fabs(xim) >= M_16_EPSILON) {
      im = copysign(M_PI_2, xim);
      if (fabs(xim) <= 1)
        re = 1 / xre;
      else if (fabs(xre) <= 1)
        re = xre / xim / xim;
      else {
        var h = hypot(xre / 2, xim / 2);
        re = xre / h / h / 4;
      }
    } else {
      if (fabs(xre) == 1
          && fabs(xim) < M_EPSILON_SQ)
        re = (copysign(0.5, xre)
                      * (M_LN2 - Math.log(fabs(xim))));
      else {
        var i2 = 0;
        if (fabs(xim) >= M_EPSILON_SQ)
          i2 = xim * xim;

        var num = 1 + xre;
        num = i2 + num * num;

        var den = 1 - xre;
        den = i2 + den * den;

        var f = num / den;
        if (f < 0.5)
          re = 0.25 * Math.log(f);
        else {
          num = 4 * xre;
          re = 0.25 * Math.log1p(num / den);
        }
      }

      var absx, absy, den;

      absx = fabs(xre);
      absy = fabs(xim);
      if (absx < absy) {
        var t = absx;
        absx = absy;
        absy = t;
      }

      if (absy < M_EPSILON_2) {
        den = (1 - absx) * (1 + absx);
        if (den == 0)
          den = 0;
      } else if (absx >= 1)
      den = (1 - absx) * (1 + absx) - absy * absy;
      else if (absx >= 0.75 || absy >= 0.5)
        den = -x2y2m1(absx, absy);
      else
        den = (1 - absx) * (1 + absx) - absy * absy;

      im = 0.5 * Math.atan2(2 * xim, den);
    }
  }

  return new Complex(re,im);
}

Complex.prototype.atanhEq = function() {
  var xre = this.re;
  var xim = this.im;
  var re;
  var im;
  var rcls = fpclassify(xre);
  var icls = fpclassify(xim);

  if (rcls <= FP_INFINITE || icls <= FP_INFINITE) {
    if (icls == FP_INFINITE) {
      re = copysign(0, xre);
      im = copysign(M_PI_2, xim);
    } else if (rcls == FP_INFINITE || rcls == FP_ZERO) {
      re = copysign(0, xre);
      if (icls >= FP_ZERO)
        im = copysign(M_PI_2, xim);
      else
        im = NaN;
    } else {
      re = NaN;
      im = NaN;
    }
  } else if (rcls == FP_ZERO && icls == FP_ZERO) {
    re = xre;
    im = xim;
  } else {
    if (fabs(xre) >= M_16_EPSILON
        || fabs(xim) >= M_16_EPSILON) {
      im = copysign(M_PI_2, xim);
      if (fabs(xim) <= 1)
        re = 1 / xre;
      else if (fabs(xre) <= 1)
        re = xre / xim / xim;
      else {
        var h = hypot(xre / 2, xim / 2);
        re = xre / h / h / 4;
      }
    } else {
      if (fabs(xre) == 1
          && fabs(xim) < M_EPSILON_SQ)
        re = (copysign(0.5, xre)
                      * (M_LN2 - Math.log(fabs(xim))));
      else {
        var i2 = 0;
        if (fabs(xim) >= M_EPSILON_SQ)
          i2 = xim * xim;

        var num = 1 + xre;
        num = i2 + num * num;

        var den = 1 - xre;
        den = i2 + den * den;

        var f = num / den;
        if (f < 0.5)
          re = 0.25 * Math.log(f);
        else {
          num = 4 * xre;
          re = 0.25 * Math.log1p(num / den);
        }
      }

      var absx, absy, den;

      absx = fabs(xre);
      absy = fabs(xim);
      if (absx < absy) {
        var t = absx;
        absx = absy;
        absy = t;
      }

      if (absy < M_EPSILON_2) {
        den = (1 - absx) * (1 + absx);
        if (den == 0)
          den = 0;
      } else if (absx >= 1)
      den = (1 - absx) * (1 + absx) - absy * absy;
      else if (absx >= 0.75 || absy >= 0.5)
        den = -x2y2m1(absx, absy);
      else
        den = (1 - absx) * (1 + absx) - absy * absy;

      im = 0.5 * Math.atan2(2 * xim, den);
    }
  }

  this.re = re;
  this.im = im;
  return this;
}


// ──── TRIGONOMETRIC RECIPROCALS ──────────────────────────────────────────────
// take the secant of a complex number

Complex.sec = function(x) {
  return x.cos().inv();
}

Complex.prototype.sec = function() {
  return this.cos().inv();
}

Complex.prototype.secEq = function() {
  return this.cosEq().invEq();
}

// take the cosecant of a complex number

Complex.csc = function(x) {
  return x.sin().inv();
}

Complex.prototype.csc = function() {
  return this.sin().inv();
}

Complex.prototype.cscEq = function() {
  return this.sinEq().invEq();
}

// take the cotangent of a complex number

Complex.cot = function(x) {
  return x.tan().inv();
}

Complex.prototype.cot = function() {
  return this.tan().inv();
}

Complex.prototype.cotEq = function() {
  return this.tanEq().invEq();
}

// take the hyperbolic secant of a complex number

Complex.sech = function(x) {
  return x.cosh().inv();
}

Complex.prototype.sech = function(x) {
  return this.cosh().inv();
}

Complex.prototype.sechEq = function(x) {
  return this.coshEq().invEq();
}

// take the hyperbolic cosecant of a complex number

Complex.csch = function(x) {
  return x.sinh().inv();
}

Complex.prototype.csch = function(x) {
  return this.sinh().inv();
}

Complex.prototype.cschEq = function(x) {
  return this.sinhEq().invEq();
}

// take the hyperbolic cotangent of a complex number

Complex.coth = function(x) {
  return x.tanh().inv();
}

Complex.prototype.coth = function(x) {
  return this.tanh().inv();
}

Complex.prototype.cothEq = function(x) {
  return this.tanhEq().invEq();
}

// take the inverse secant of a complex number

Complex.asec = function(x) {
  return x.inv().acos();
}

Complex.prototype.asec = function() {
  return this.inv().acos();
}

Complex.prototype.asecEq = function() {
  return this.invEq().acosEq();
}

// take the inverse cosecant of a complex number

Complex.acsc = function(x) {
  return x.inv().asin();
}

Complex.prototype.acsc = function() {
  return this.inv().asin();
}

Complex.prototype.acscEq = function() {
  return this.invEq().asinEq();
}

// take the inverse cotangent of a complex number

Complex.acot = function(x) {
  return x.atan().scalarSub(M_PI_2);
}

Complex.prototype.acot = function() {
  return this.atan().scalarSub(M_PI_2);
}

Complex.prototype.acotEq = function() {
  return this.atanEq().scalarSubEq(M_PI_2);
}

// take the inverse hyperbolic secant of a complex number

Complex.asech = function(x) {
  return x.inv().acosh();
}

Complex.prototype.asech = function() {
  return this.inv().acosh();
}

Complex.prototype.asechEq = function() {
  return this.invEq().acoshEq();
}

// take the inverse hyperbolic cosecant of a complex number

Complex.acsch = function(x) {
  return x.inv().asinh();
}

Complex.prototype.acsch = function() {
  return this.inv().asinh();
}

Complex.prototype.acschEq = function() {
  return this.invEq().asinhEq();
}

// take the inverse hyperbolic cotangent of a complex number

Complex.acoth = function(x) {
  return x.inv().atanh();
}

Complex.prototype.acoth = function() {
  return this.inv().atanh();
}

Complex.prototype.acothEq = function() {
  return this.invEq().atanhEq();
}

// ──── COMPONENT OPERATIONS ───────────────────────────────────────────────────
/* Compute complex component floor */
Complex.floor = function(x) {
  return new Complex(
    Math.floor(x.re),
    Math.floor(x.im),
  );
}

Complex.prototype.floor = function() {
  return new Complex(
    Math.floor(this.re),
    Math.floor(this.im),
  );
}

Complex.prototype.floorEq = function() {
  this.re = Math.floor(this.re);
  this.im = Math.floor(this.im);
  return this;
}

/* Compute complex component ceiling */
Complex.ceil = function(x) {
  return new Complex(
    Math.ceil(x.re),
    Math.ceil(x.im),
  );
}

Complex.prototype.ceil = function() {
  return new Complex(
    Math.ceil(this.re),
    Math.ceil(this.im),
  );
}

Complex.prototype.ceilEq = function() {
  this.re = Math.ceil(this.re);
  this.im = Math.ceil(this.im);
  return this;
}

/* Compute complex component truncation */
Complex.trunc = function(x) {
  return new Complex(
    Math.trunc(x.re),
    Math.trunc(x.im),
  );
}

Complex.prototype.trunc = function() {
  return new Complex(
    Math.trunc(this.re),
    Math.trunc(this.im),
  );
}

Complex.prototype.truncEq = function() {
  this.re = Math.trunc(this.re);
  this.im = Math.trunc(this.im);
  return this;
}

/* Compute complex component round */
Complex.round = function(x) {
  return new Complex(
    Math.round(x.re),
    Math.round(x.im),
  );
}

Complex.prototype.round = function() {
  return new Complex(
    Math.round(this.re),
    Math.round(this.im),
  );
}

Complex.prototype.roundEq = function() {
  this.re = Math.round(this.re);
  this.im = Math.round(this.im);
  return this;
}

/* Compute complex component fround */
Complex.fround = function(x) {
  return new Complex(
    Math.fround(x.re),
    Math.fround(x.im),
  );
}

Complex.prototype.fround = function() {
  return new Complex(
    Math.fround(this.re),
    Math.fround(this.im),
  );
}

Complex.prototype.froundEq = function() {
  this.re = Math.fround(this.re);
  this.im = Math.fround(this.im);
  return this;
}

/* Compute complex mod */
// mod(z,w) = z - w*floor(z/w)
Complex.mod = function(z, w) {
  return z.sub(z.div(w).floor().mul(w));
}

Complex.prototype.mod = function(z) {
  return this.sub(this.div(z).floor().mul(z));
}

Complex.prototype.modEq = function(z) {
  return this.subEq(this.div(z).floor().mul(z));
}

/* Compute complex component random */
Complex.random = function() {
  return new Complex(Math.random(), Math.random());
}

Complex.setRandom = function(z) {
  z.re = Math.random();
  z.im = Math.random();
  return z;
}

Complex.prototype.setRandom = function() {
  this.re = Math.random();
  this.im = Math.random();
  return this;
}

/* Compute complex sum */
Complex.sum = function(...z) {
  return z.reduce(Complex.add, new Complex(0,0));
}

Complex.prototype.sum = function(...z) {
  return z.reduce(Complex.add, this.clone());
}

Complex.prototype.sumEq = function(...z) {
  return this.addEq(z.reduce(Complex.add, new Complex(0,0)));
}

/* Compute complex product */
Complex.prod = function(...z) {
  return z.reduce(Complex.mul, new Complex(1,0));
}

Complex.prototype.prod = function(...z) {
  return z.reduce(Complex.mul, this.clone());
}

Complex.prototype.prodEq = function(...z) {
  return this.mulEq(z.reduce(Complex.mul, new Complex(1,0)));
}

/* Compute complex component minimum */
Complex.min = function(...z) {
  var x = Infinity;
  var y = Infinity;
  var l = z.length;
  var i = 0;
  for (;i < l;) {
    if (z[i].re < x) x = z[i].re;
    if (z[i].im < y) y = z[i].im;
    i++;
  }
  return new Complex(x, y);
}

Complex.prototype.min = function(...z) {
  var x = this.re;
  var y = this.im;
  var l = z.length;
  var i = 0;
  for (;i < l;) {
    if (z[i].re < x) x = z[i].re;
    if (z[i].im < y) y = z[i].im;
    i++;
  }
  return new Complex(x, y);
}

Complex.prototype.minEq = function(...z) {
  var x = this.re;
  var y = this.im;
  var l = z.length;
  var i = 0;
  for (;i < l;) {
    if (z[i].re < x) x = z[i].re;
    if (z[i].im < y) y = z[i].im;
    i++;
  }
  this.re = x;
  this.im = y;
  return this;
}

/* Compute complex component maximum */
Complex.max = function(...z) {
  var x = -Infinity;
  var y = -Infinity;
  var l = z.length;
  var i = 0;
  for (;i < l;) {
    if (z[i].re > x) x = z[i].re;
    if (z[i].im > y) y = z[i].im;
    i++;
  }
  return new Complex(x, y);
}

Complex.prototype.max = function(...z) {
  var x = this.re;
  var y = this.im;
  var l = z.length;
  var i = 0;
  for (;i < l;) {
    if (z[i].re > x) x = z[i].re;
    if (z[i].im > y) y = z[i].im;
    i++;
  }
  return new Complex(x, y);
}

Complex.prototype.maxEq = function(...z) {
  var x = this.re;
  var y = this.im;
  var l = z.length;
  var i = 0;
  for (;i < l;) {
    if (z[i].re > x) x = z[i].re;
    if (z[i].im > y) y = z[i].im;
    i++;
  }
  this.re = x;
  this.im = y;
  return this;
}

// ──── BOOLEAN TESTS ──────────────────────────────────────────────────────────
/* Complex equality test */
Complex.equal = function(z, w, t=M_EPSILON) {
  if (z === w) return true;
  var dre = z.re - w.re;
  var dim = z.im - w.im;
  return dre*dre + dim*dim < t*t;
}

Complex.prototype.equal = function(z, t=M_EPSILON) {
  if (this === z) return true;
  var dre = this.re - z.re;
  var dim = this.im - z.im;
  return dre*dre + dim*dim < t*t;
}

/* Complex NaN test */
Complex.isNaN = function(z) {
  var rcls = fpclassify(z.re);
  var icls = fpclassify(z.im);
  return rcls === FP_NAN || icls === FP_NAN;
}

/* Complex finite test */
Complex.isFinite = function(z) {
  var rcls = fpclassify(z.re);
  var icls = fpclassify(z.im);
  return rcls > FP_INFINITE && icls > FP_INFINITE;
}

/* Complex zero test */
Complex.isZero = function(z) {
  var rcls = fpclassify(z.re);
  var icls = fpclassify(z.im);
  return rcls === FP_ZERO && icls === FP_ZERO;
}

/* Complex real test */
Complex.isReal = function(z) {
  var icls = fpclassify(z.im);
  return icls === FP_ZERO;
}

// ──── SPECIAL FUNCTIONS ──────────────────────────────────────────────────────

var G = 7;
var G0 =  0.99999999999980993;
var G1 =  676.5203681218851;
var G2 = -1259.1392167224028;
var G3 =  771.32342877765313;
var G4 = -176.61502916214059;
var G5 =  12.507343278686905;
var G6 = -0.13857109526572012;
var G7 =  9.9843695780195716e-6;
var G8 =  1.5056327351493116e-7;

function gammaLanczos(x, y) {
  var y2 = y * y;
  var w = x - 0.5;

  // t = z + g - 0.5
  var t = G + w;

  // log(t)
  var ltx = Math.log(hypot(t, y));
  var lty = Math.atan2(y, t);

  // t^(z - 0.5)
  var twp = w * lty + y * ltx;
  var twr = Math.exp(w * ltx - y * lty);
  var twx = twr * Math.cos(twp);
  var twy = twr * Math.sin(twp);

  // e^{-t}
  var et = Math.exp(-t);
  var etx = et * Math.cos(y);
  var ety = -et * Math.sin(y);

  // Sum: C[0] + Σ_{k=1}^{8} C[k] / (z + k - 1)
  var sx = G0;
  var sy = 0;
  var d = 1 / (x * x + y2); sx += G1 * x++ * d; sy -= G1 * y * d;
  d = 1 / (x * x + y2); sx += G2 * x++ * d; sy -= G2 * y * d;
  d = 1 / (x * x + y2); sx += G3 * x++ * d; sy -= G3 * y * d;
  d = 1 / (x * x + y2); sx += G4 * x++ * d; sy -= G4 * y * d;
  d = 1 / (x * x + y2); sx += G5 * x++ * d; sy -= G5 * y * d;
  d = 1 / (x * x + y2); sx += G6 * x++ * d; sy -= G6 * y * d;
  d = 1 / (x * x + y2); sx += G7 * x++ * d; sy -= G7 * y * d;
  d = 1 / (x * x + y2); sx += G8 * x++ * d; sy -= G8 * y * d;

  var Bx = twx * etx - twy * ety;
  var By = twx * ety + twy * etx;

  return new Complex(
    M_SQRT2PI * (Bx * sx - By * sy),
    M_SQRT2PI * (Bx * sy + By * sx),
  ) 
}

// take the gamma function of a complex number
// static gamma
// gamma
// gammaEq

Complex.gamma = function(z) {
  var x = z.re;
  var y = z.im;
  var x0 = Math.round(x);

  if (x <= 0 && fabs(y) < 1e-15 && fabs(x - x0) < 1e-15) {
      return new Complex(Infinity, 0);
  }

  if (x < 0.5) {
    // Γ(z) = π / ( sin(πz) * Γ(1-z) )
    var g = gammaLanczos(1 - x, -y);
    var gx = g.re;
    var gy = g.im;

    // sin(πz) = sin(πx)cosh(πy) + i cos(πx)sinh(πy)
    var pix = M_PI * x;
    var piy = M_PI * y;
    var spx = Math.sin(pix) * Math.cosh(piy);
    var spy = Math.cos(pix) * Math.sinh(piy);

    // π / (sin(πz) * Gamma(1-z))
    var dx = spx * gx - spy * gy;
    var dy = spx * gy + spy * gx;
    var d = dx * dx + dy * dy;
    var gx = M_PI * dx / d;
    var gy = -M_PI * dy / d;

    if (gx === Infinity || gx === -Infinity ||
      gy === Infinity || gy === -Infinity) {
      if (fabs(y) < 1e-12 && x <= 0 && fabs(x - x0) < 1e-12) {
        g.re = Infinity;
        g.im = 0;
        return g;
      }
    }

    g.re = gx;
    g.im = gy;
    return g;
  }

  return gammaLanczos(x, y);
}

Complex.prototype.gamma = function() {
  return Complex.gamma(this);
}

Complex.prototype.gammaEq = function() {
  var g = Complex.gamma(this);
  this.re = g.re;
  this.im = g.im;
  return this;
}

// take the factorial of a complex number
// static gamma
// gamma
// gammaEq

Complex.fact = function(z) {
  return Complex.gamma(z.addScalar(1));
}

Complex.prototype.fact = function() {
  return this.addScalar(1).gamma();
}

Complex.prototype.factEq = function() {
  return this.addScalarEq(1).gammaEq();
}

// take the beta function of a complex numbers
// static beta
// beta
// betaEq

var cadd = Complex.add;
var cmul = Complex.mul;
var cdiv = Complex.div;
var cgamma = Complex.gamma;

Complex.beta = function(z, w) {
  return cdiv(cmul(cgamma(z),cgamma(w)),cgamma(cadd(z,w)));
}

Complex.prototype.beta = function(z) {
  return cdiv(cmul(cgamma(this),cgamma(z)),cgamma(cadd(this,z)));
}

Complex.prototype.betaEq = function(z) {
  var b = cdiv(cmul(cgamma(this),cgamma(z)),cgamma(cadd(this,z)));
  this.re = b.re;
  this.im = b.im;
  return this;
}

// take the binomial coefficient of a complex numbers
// static binom
// binom
// binomEq

var cadd = Complex.add;
var caddn = Complex.addScalar;
var csub = Complex.sub;
var cmul = Complex.mul;
var cdiv = Complex.div;
var cgamma = Complex.gamma;

Complex.binom = function(z, w) {
  return cdiv(cgamma(caddn(z,1)),cmul(cgamma(caddn(w,1)),cgamma(caddn(csub(z,w),1))));
}

Complex.prototype.binom = function(z) {
  return cdiv(cgamma(caddn(this,1)),cmul(cgamma(caddn(z,1)),cgamma(caddn(csub(this,z),1))));
}

Complex.prototype.binomEq = function(z) {
  var b = cdiv(cgamma(caddn(this,1)),cmul(cgamma(caddn(z,1)),cgamma(caddn(csub(this,z),1))));
  this.re = b.re;
  this.im = b.im;
  return this;
}

var E0 = 1.119;
var E1 = 1 / 12;
var E2 = 7 / 480;
var E3 = 5 / 896;
var E4 = 787 / 276480;
function rerf(x) {
  var k = 1 - Math.exp(-x * x);
  return E0 * Math.sqrt(k) * (1 - k * (E1 + k * (E2 + k * (E3 + k * E4))));
}

function cerf_large(x, y, out) {
  var x2  = x * x;
  var y2  = y * y;
  var d   = 1 / (x2 + y2);
  var kx  = y * d;
  var ky  = x * d;
  var sx  = 1 + 0.5 * (kx * kx - ky * ky);
  var sy  = kx * ky;
  var cr  = (kx * sx - ky * sy) * M_1_SQRTPI;
  var ci  = (kx * sy + ky * sx) * M_1_SQRTPI;
  var er  = Math.exp(y2 - x2);
  var ep  = -2 * x * y;
  var ex  = er*Math.cos(ep);
  var ey  = er*Math.sin(ep);
  out.re  = 1 - ex * ci + ey * cr;
  out.im  = ex * cr - ey * ci;
}

function cerf_small(x, y, out) {
  x = Math.max(x, 1e-15);
  var x2     = x * x;
  var xy     = 2 * x * y;
  var K      = Math.exp(-x2) * M_1_PI;
  var q      = 4 * x2;
  var a      = Math.cos(xy);
  var b      = Math.sin(xy);
  var offset = y > x ? (y + x) * (y - x) : 0;
  var scale  = Math.exp(-offset);
  var sx = 0;
  var sy = 0;
  for (var n = 1; n <= 16; n++) {
    var kk  = n*n*0.25 + x2 + offset;
    var kz  = n * y;
    var e1  = Math.exp( kz - kk);
    var e2 = Math.exp(-kz - kk);
    var aa  = x*(e1 + e2);
    var bb = 0.5*n*(e1 - e2);
    var iq  = 1 / (n*n + q);
    sx += (2*x*Math.exp(-kk) - a*aa + b*bb) * iq;
    sy += (b*aa + a*bb) * iq;
  }
  var eoff = Math.exp(offset)*M_2_PI;
  var k2x  = K / (2 * x);
  out.re = scale*(rerf(x) + k2x*(1 - a)) + sx*eoff;
  out.im = scale*(k2x * b)               + sy*eoff;
}

Complex.erf = function(z) {
  out = new Complex(0, 0);
  var x  = z.re;
  var y  = z.im;
  var ax = fabs(x);
  var ay = fabs(y);
  if (ay > 5.5) cerf_large(ax, ay, out);
  else          cerf_small(ax, ay, out);
  if (y < 0) out.im *= -1;
  if (x < 0) out.re *= -1;
  return out;
}

Complex.prototype.erf = function() {
  out = this.clone();
  var x  = this.re;
  var y  = this.im;
  var ax = fabs(x);
  var ay = fabs(y);
  if (ay > 5.5) cerf_large(ax, ay, out);
  else          cerf_small(ax, ay, out);
  if (y < 0) out.im *= -1;
  if (x < 0) out.re *= -1;
  return out;
}

Complex.prototype.erfEq = function() {
  var x  = this.re;
  var y  = this.im;
  var ax = fabs(x);
  var ay = fabs(y);
  if (ay > 5.5) cerf_large(ax, ay, this);
  else          cerf_small(ax, ay, this);
  if (y < 0) this.im *= -1;
  if (x < 0) this.re *= -1;
  return this;
}

function cw(x, y, out) {
  out = out ?? new Complex(0, 0);

  var r   = Math.sqrt(x*x + y*y);
  var l1x = Math.log(r < 1e-15 ? 1e-15 : r);
  var l1y = Math.atan2(y, x);

  if (y < -8) { out.re = x; out.im = y; return out; }

  var zx = Math.abs(x), zy = Math.abs(y);
  var zr = Math.sqrt(zx*zx + zy*zy);

  var ex, ey;
  if (zr > 3) {
    var lr  = Math.sqrt(l1x*l1x + l1y*l1y);
    var llx = Math.log(lr < 1e-15 ? 1e-15 : lr);
    var lly = Math.atan2(l1y, l1x);
    ex = l1x - llx;
    ey = l1y - lly;
  } else {
    var sx = Math.E * zx + 1;
    var sy = Math.E * zy;
    var sr = Math.sqrt(Math.sqrt(sx*sx + sy*sy));
    var st = Math.atan2(sy, sx) * 0.5;
    ex = sr * Math.cos(st) - 1;
    ey = sr * Math.sin(st);
  }


  if (r < M_1_E) {
    for (var i = 0; i < 8; i++) {
      var er   = Math.exp(ex);
      var erx  = er * Math.cos(ey);
      var ery  = er * Math.sin(ey);
      var wewx = ex*erx - ey*ery;
      var wewy = ex*ery + ey*erx;
      var fx   = wewx - x;
      var fy   = wewy - y;
      var d1x  = (ex+1)*erx - ey*ery;
      var d1y  = (ex+1)*ery + ey*erx;
      var dn   = 1 / (d1x*d1x + d1y*d1y);
      ex -= (fx*d1x + fy*d1y) * dn;
      ey -= (fy*d1x - fx*d1y) * dn;
    }
  } else {
    for (var i = 0; i < 10; i++) {
      var er  = Math.sqrt(ex*ex + ey*ey);
      var lx  = Math.log(er < 1e-15 ? 1e-15 : er);
      var ly  = Math.atan2(ey, ex);
      var ax  = ex + lx - l1x;
      var ay  = ey + ly - l1y;
      var bx  = ex*ax - ey*ay;
      var by  = ex*ay + ey*ax;
      var d1x = 1 + ex, d1y = ey;
      var dn  = 1 / (d1x*d1x + d1y*d1y);
      ex -= (bx*d1x + by*d1y) * dn;
      ey -= (by*d1x - bx*d1y) * dn;
    }
  }

  out.re = ex;
  out.im = ey;
  return out;
}

Complex.lambertw = function(z) {
  return cw(z.re, z.im);
}

Complex.prototype.lambertw = function() {
  return cw(this.re, this.im);
}

Complex.prototype.lambertwEq = function() {
  return cw(this.re, this.im, this);
}

// ──── EXPORTS ────────────────────────────────────────────────────────────────

if (typeof module !== 'undefined') {
  var util = require('util');
  Complex.prototype[util.inspect.custom] = function(depth, options) {
    return options.stylize(this.toString(),'number');
  }
  module.exports = Complex;
}
else globalThis.Complex = Complex;

})();// ──── END IIFE ──────────────────────────────────────────────────────────
