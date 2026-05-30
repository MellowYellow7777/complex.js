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

## Properties

### Static

todo

### Instance

data properties:
```javascript
re
im
```
pseudo-properties:
```javascript
get real()
set real(x)
get imag()
set imag(y)
get r()
set r(r)
get theta()
set theta(phi)
```

## Included Functions

### Constructors

create a new complex number given euclidean components:
```javascript
constructor(x, y)
static create(x=0, y=0)
static fromScalar(n=0)
atatic fromArray(array, offset=0)
```
create a new complex number given polar components:
```javascript
static polar(r=0,phi=0)
static cis(phi=0)
static fromArrayPolar(array,offset=0)
```
### Copy & Clone

create a clone of a complex number with the same components:
```javascript
static clone(z)
clone()
```
copy the components of a complex number to another:
```javascript
static clone(z)     // new Complex = z
static copy(z, w)   // z = w
static copyTo(z, w) // w = z
clone()             // new Complex = z
copy(z)             // this = z
copyTo(z)           // z = this
```
### Type Conversion

write a complex numbers components into a new or existing array:
```javascript
static toArray(z, array=[], offset=0)
static toArrayPolar(z,array=[], offset=0)
toArray(array=[], offset=0)
toArrayPolar(array=[], offset=0)
```
get a string representation of a complex number
```javascript
toString()
```
### Setters

set the euclidean real/imaginary components of a complex number:
```javascript
static setReal(z, x)
static setImag(z, y)
static set(z, x, y)
static setScalar(z, n)
static setComponent(z, i, n)
static setFromArray(z, array, offset=0)
setReal(x)
setImag(y)
set(x, y)
setScalar(n)
setComponent(i, n)
setFromArray(array, offset=0)
```
set the polar magnitude/argument components of a complex number:
```javascript
static setAbs(z, r)
static setArg(z, phi)
static setPolar(z, x, y)
static setFromArrayPolar(z, array, offset=0)
setAbs(r)
setArg(phi)
setPolar(r, phi)
setFromArrayPolar(array, offset=0)
```
### Setters

get the euclidean real/imaginary components of a complex number
```javascript
static getReal(z)
static getImage(z)
static getComponent(z, i)
static complexReal(z)
static complexImag(z)
getReal()
getImage()
getComponent(i)
getComplexReal()
getComplexImag()
```
get the polar magnitude/argument components of a complex number:
```javascript
static getAbs(z)
static getArg(z)
static getComplexAbs(z)
static getComplexArg(z)
getAbs()
getArg()
getComplexAbs()
getComplexArg()
```
### Operations

#### Additive Operations

add complex numbers:
```javascript
static add(z, w)
add(z)
addEq(z)
```
add a scalar to a complex number:
```javascript
static addScalar(z,n)
addScalar(n)
addScalarEq(n)
```

subtract complex numbers:
```javascript
static sub(z, w)
sub(z)
subEq(z)
```
subtract a scalar from a complex number:
```javascript
static subScalar(z, n)
subScalar(n)
subScalarEq(n)
```
subtract a complex number from a scalar:
```javascript
static scalarSub(n, z)
scalarSub(n)
scalarSubEq(n)
```
negate a complex number:
```javascript
static neg(z)
neg()
negEq()
```
take the conjugate of a complex number:
```javascript
static conj(z)
conj()
conjEq()
```
#### Multiplicative Operations

multiply complex numbers:
```javascript
static mul(z, w)
mul(z)
mulEq(z)
```
multiply a complex number by a scalar:
```javascript
static mulScalar(z, n)
mulScalar(n)
mulScalarEq(n)
```
divide complex numbers:
```javascript
static div(z, w)
div(z)
divEq(z)
```
divide a complex number by a scalar:
```javascript
static divScalar(z, n)
divScalar(n)
divScalarEq(n)
```
divide a scalar by a complex number:
```javascript
static scalarDiv(z, n)
scalarDiv(n)
scalarDivEq(n)
```
take the reciprocal of a complex number:
```javascript
static inv(z)
inv()
invEq()
```
take the directed sign of a complex number:
```javascript
static sgn(z)
sgn()
sgnEq()
```
take the projection a complex number on the riemann sphere:
```javascript
static proj(z)
prog()
projEq()
```
#### Exponential and Logarithmic Operations

take e to the power of a complex number:
```javascript
static exp(z)
exp()
expEq()
```
take the a logarithm of a complex number:
```javascript
static log(z)
static log10(z)
static log2(z)
log()
log10()
log2()
logEq()
log10Eq()
log2Eq()
```
take the complex power of a complex number:
```javascript
static pow(z, w)
pow(z)
powEq(z)
```
take the principal square root of a complex number:
```javascript
static sqrt(z)
sqrt()
sqrtEq()
```
#### Trigonometric Operations

take the sine, cosine, and tangent of a complex number:
```javascript
static sin(z)
static cos(z)
static tan(z)
sin()
cos()
tan()
sinEq()
cosEq()
tanEq()
```
take the hyperbolic sine, cosine, and tangent of a complex number:
```javascript
static sinh(z)
static cosh(z)
static tanh(z)
sinh()
cosh()
tanh()
sinhEq()
coshEq()
tanhEq()
```
take the inverse sine, cosine, and tangent of a complex number:
```javascript
static asin(z)
static acos(z)
static atan(z)
asin()
acos()
atan()
asinEq()
acosEq()
atanEq()
```
take the inverse hyperbolic sine, cosine, and tangent of a complex number:
```javascript
static asinh(z)
static acosh(z)
static atanh(z)
asinh()
acosh()
atanh()
asinhEq()
acoshEq()
atanhEq()
```
#### Trigonometric Reciprocals

take the secant, cosecant, and cotangent of a complex number:
```javascript
static sec(z)
static csc(z)
static cot(z)
sec()
csc()
cot()
secEq()
cscEq()
cotEq()
```
take the hyperbolic secant, cosecant, and cotangent of a complex number:
```javascript
static sech(z)
static csch(z)
static coth(z)
sech()
csch()
coth()
sechEq()
cschEq()
cothEq()
```
take the inverse secant, cosecant, and cotangent of a complex number:
```javascript
static asec(z)
static acsc(z)
static acot(z)
asec()
acsc()
acot()
asecEq()
acscEq()
acotEq()
```
take the inverse hyperbolic secant, cosecant, and cotangent of a complex number:
```javascript
static asech(z)
static acsch(z)
static acoth(z)
asech()
acsch()
acoth()
asechEq()
acschEq()
acothEq()
```
#### Component-Wise Operations

rounding functions:
```javascript
static floor(z)
static ceil(z)
static trunc(z)
static round(z)
static fround(z)
floor()
ceil()
trunc()
round()
fround()
floorEq()
ceilEq()
truncEq()
roundEq()
froundEq()
```
complex mod:
```javascript
static mod(z, w)
mod()
modEq()
```
random:
```javascript
static random()
static setRandom(z)
setRandom()
```
variadic operations:
```javascript
static sum(...z)
static prod(...z)
static min(...z)
static max(...z)
sum(...z)
prod(...z)
min(...z)
max(...z)
sumEq(...z)
prodEq(...z)
minEq(...z)
maxEq(...z)
```
#### Boolean Tests

check if two complex numbers are equal:
```javascript
static equal(z, w, t=EPSILON)
equal(z, t=EPSILON)
```
check the classification of a complex number:
```javascript
static isNaN(z)
static isFinite(z)
static isZero(z)
static isReal(z)
```
#### Special Functions

```javascript
static erf(z)
static gamma(z)
static fact(z)
static beta(z, w)
static binom(z, w)
static lambertw(z)
erf()
gamma()
fact()
beta(z)
binom(z)
lambertw()
erfEq()
gammaEq()
factEq()
betaEq(z)
binomEq(z)
lambertwEq()
```
