Open up dev tools, and you can use the library to defined functions to draw to the screen via domain coloring. You can set a custom function by assigning to the 'f' variable. For example:

```javascript
f = z => Complex.addScalar(Complex.square(z),1);
```

to set the function to z^2 + 1. To draw, just call 'draw':

```javascript
draw();
```

and you will see the domain color plot of your function. You can also provide a scale argument to zoom in or out. This is a rudementary example but a good way to test out the library.
