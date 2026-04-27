# complex.js
complex number library

usage (node):

```javascript
Complex = require('./complex.js');
```

usage (web):

```html
<script src="complex.js"></script>
```

The base of library is a direct port of glibc complex number function implementations. These functions include:

- creal, cimag, cabs, carg, conj, cproj
- cexp, clog, clog10, cpow, csqrt
- csin, ccos, ctan, csinh, ccosh, ctanh
- casin, cacos, catan, casinh, cacosh, catanh

Each of which should behave identical to the original functions to the extent possible in JavaScript. Some things, such as floating-point exception flags, have been left out in favor of practicality. The corresponding functions have been left in the main file but are defined as no-op. This is done in attempt to leave the functions as un-altered as possible.

Each of the constants and helper functions have also been named as they appear in the original code. This introduces some redundacy and additional overhead, though provides a good working model to experiment with in a JavaScript environment.

Much additional behavior and functions have been inluded on top of the ones mentioned above. These are not glibc ports and many of which are implemented via naive methods and haven't been optimized much.

The data structure within the main file has complex numbers all stored in a single Float64Array buffer. Complex numbers are then represented by their index where they are stored in this buffer. All functions expect complex number arguments to be provided in this form. This README and related documentation are a work in progress and are intended on being updated further in the future. The following is an AI generated write-up that can be referenced for more details.

# complex_glibc.js
 
A complex number library for JavaScript, ported from glibc's math implementation with additional functions. Uses a typed array buffer for allocation-free arithmetic.
 
## Design
 
All complex numbers are stored as pairs of `Float64` values in a shared `Float64Array` buffer. Functions return an integer index into this buffer rather than allocating objects. This avoids GC pressure in tight loops such as domain coloring renderers or iterative algorithms.
 
```js
var z = Complex.create(1, 2);   // returns buffer index
var w = Complex.exp(z);         // returns buffer index
console.log(Complex.toString(w));
```
 
Call `Complex.resetBuffer()` between frames or computation passes to reclaim space.
 
## Installation
 
```html
<script src="complex_glibc.js"></script>
<!-- Complex is available as a global -->
```
 
```js
// Node
const Complex = require('./complex_glibc.js');
```
 
## API
 
### Buffer
 
| Function | Description |
|---|---|
| `Complex.resetBuffer()` | Reset buffer index to 0 |
| `Complex.BUFFER` | The underlying `Float64Array` |
| `Complex.BUFFER_SIZE` | Current buffer size (settable) |
| `Complex.INDEX` | Current buffer index (settable) |
 
### Create / Inspect
 
| Function | Description |
|---|---|
| `Complex.create(re, im)` | Create complex number |
| `Complex.real(z)` | Real part |
| `Complex.imag(z)` | Imaginary part |
| `Complex.abs(z)` | Absolute value (modulus) |
| `Complex.arg(z)` | Argument (angle) |
| `Complex.toString(z)` | String representation |
 
### Arithmetic
 
| Function | Description |
|---|---|
| `Complex.add(z, w)` | z + w |
| `Complex.sub(z, w)` | z - w |
| `Complex.mul(z, w)` | z * w |
| `Complex.div(z, w)` | z / w |
| `Complex.addScalar(z, n)` | z + n |
| `Complex.subScalar(z, n)` | z - n |
| `Complex.mulScalar(z, n)` | z * n |
| `Complex.divScalar(z, n)` | z / n |
| `Complex.neg(z)` | -z |
| `Complex.inv(z)` | 1/z |
| `Complex.conj(z)` | Conjugate |
| `Complex.sgn(z)` | Sign (unit vector) |
| `Complex.square(z)` | z² |
| `Complex.cube(z)` | z³ |
 
### Exponential / Logarithm
 
| Function | Description |
|---|---|
| `Complex.exp(z)` | eᶻ |
| `Complex.log(z)` | ln(z) |
| `Complex.log10(z)` | log₁₀(z) |
| `Complex.log2(z)` | log₂(z) |
| `Complex.pow(z, w)` | zʷ |
| `Complex.sqrt(z)` | √z |
 
### Trigonometric
 
| Function | Description |
|---|---|
| `Complex.sin(z)` | sin(z) |
| `Complex.cos(z)` | cos(z) |
| `Complex.tan(z)` | tan(z) |
| `Complex.sec(z)` | sec(z) |
| `Complex.csc(z)` | csc(z) |
| `Complex.cot(z)` | cot(z) |
| `Complex.asin(z)` | arcsin(z) |
| `Complex.acos(z)` | arccos(z) |
| `Complex.atan(z)` | arctan(z) |
| `Complex.asec(z)` | arcsec(z) |
| `Complex.acsc(z)` | arccsc(z) |
| `Complex.acot(z)` | arccot(z) |
 
### Hyperbolic
 
| Function | Description |
|---|---|
| `Complex.sinh(z)` | sinh(z) |
| `Complex.cosh(z)` | cosh(z) |
| `Complex.tanh(z)` | tanh(z) |
| `Complex.sech(z)` | sech(z) |
| `Complex.csch(z)` | csch(z) |
| `Complex.coth(z)` | coth(z) |
| `Complex.asinh(z)` | arcsinh(z) |
| `Complex.acosh(z)` | arccosh(z) |
| `Complex.atanh(z)` | arctanh(z) |
| `Complex.asech(z)` | arcsech(z) |
| `Complex.acsch(z)` | arccsch(z) |
| `Complex.acoth(z)` | arccoth(z) |
 
### Special Functions
 
| Function | Description |
|---|---|
| `Complex.gamma(z)` | Γ(z) via Lanczos approximation |
| `Complex.fact(z)` | z! = Γ(z+1) |
| `Complex.beta(z, w)` | B(z,w) = Γ(z)Γ(w)/Γ(z+w) |
| `Complex.binom(z, w)` | Binomial coefficient C(z,w) |
| `Complex.sinc(z)` | Normalized sinc |
| `Complex.cis(x)` | cos(x) + i·sin(x) |
| `Complex.mod(z, w)` | Complex modulo |
 
### Utility
 
| Function | Description |
|---|---|
| `Complex.derivative(f, h?)` | Numerical derivative, default h=1e-7 |
| `Complex.equal(z, w, t?)` | Equality within tolerance |
| `Complex.isNaN(z)` | NaN test |
| `Complex.isFinite(z)` | Finite test |
| `Complex.isZero(z)` | Zero test |
| `Complex.isReal(z)` | Real test (imaginary part is zero) |
| `Complex.random()` | Random complex in [0,1)² |
| `Complex.min(...z)` | Component-wise minimum |
| `Complex.max(...z)` | Component-wise maximum |
| `Complex.floor(z)` | Component-wise floor |
| `Complex.ceil(z)` | Component-wise ceil |
| `Complex.round(z)` | Component-wise round |
| `Complex.trunc(z)` | Component-wise truncate |
 
## Example — Domain Coloring
 
```js
function draw(f) {
  Complex.resetBuffer();
  var imageData = ctx.createImageData(canvas.width, canvas.height);
  var data = imageData.data;
  var w = canvas.width, h = canvas.height, s = Math.min(w, h);
  var i = 0;
  for (var cy = 0; cy < h; cy++) {
    for (var cx = 0; cx < w; cx++) {
      var z = Complex.create((cx - w/2) / s, (h/2 - cy) / s);
      var out = f(z);
      var hue = (Complex.arg(out) / (2*Math.PI) + 1) % 1;
      var lig = 1 - 1/(1 + Math.log1p(Complex.abs(out)));
      // HSL to RGB ...
      data[i++] = r; data[i++] = g; data[i++] = b; data[i++] = 255;
    }
  }
  ctx.putImageData(imageData, 0, 0);
}
 
draw(z => Complex.gamma(z));
```
 
## Notes
 
- Core trig, hyperbolic, exponential, and logarithm functions are ported directly from glibc and handle all IEEE 754 edge cases including signed zeros, infinities, and NaNs.
- Special functions (gamma, beta, etc.) use standard approximations and have not been exhaustively tested over all edge cases.
- The buffer is a flat `Float64Array`. If you run out of space, increase `Complex.BUFFER_SIZE` or call `resetBuffer()` more frequently.
