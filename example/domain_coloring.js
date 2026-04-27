canvas = document.createElement('canvas');
ctx = canvas.getContext('2d');
document.body.appendChild(canvas);
document.body.style.margin = '0';

function resize() {
  canvas.width = window.innerWidth;
  canvas.height = window.innerHeight;
}

resize();

window.onresize = resize;

function hue2rgb(p, q, t) {
  if (t < 0) t += 1;
  if (t > 1) t -= 1;
  if (t < 1/6) return p + (q-p)*6*t;
  if (t < 1/2) return q;
  if (t < 2/3) return p + (q-p)*(2/3-t)*6;
  return p;
}

function draw(scale = 1e-1) {
var buffer = Complex.BUFFER;
var w = canvas.width;
var h = canvas.height;
var s = Math.min(w,h);
var imageData = ctx.createImageData(w,h);
var data = imageData.data;
var i = 0;
for (var cy = 0; cy < h; cy++) {
for (var cx = 0; cx < w; cx++) {
  var x = cx - w/2;
  var y = h/2 - cy;
  x /= s * scale;
  y /= s * scale;
  Complex.resetBuffer();
  var z = Complex.create(x,y);
  var out = f(z);
  var arg = Complex.arg(out);
  var abs = Complex.abs(out);
  var hue = (arg / (2*Math.PI) + 1) % 1;
  var lig = 1 - 1/(1 + Math.log1p(abs));
  var q = lig < 0.5 ? lig * 2 : 2 - lig * 2;
  var p = 2 * lig - q;
  var r = hue2rgb(p, q, hue + 1/3);
  var g = hue2rgb(p, q, hue);
  var b = hue2rgb(p, q, hue - 1/3);
  data[i++] = r * 255 | 0;
  data[i++] = g * 255 | 0;
  data[i++] = b * 255 | 0;
  data[i++] = 255;
}}
ctx.putImageData(imageData,0,0);
}

f = Complex.acsch;

draw();
