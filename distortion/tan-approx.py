#!/usr/bin/env python3

import numpy as np
from scipy.optimize import least_squares

def approx(x,c):
  c1, c2 = c
  x2 = x * x
  return x * (c1 + c2 * x2) / (c1 + x2)

def residual(c):
  xs = np.linspace(0, 0.2501*np.pi, 10000)

  f = np.tan(xs)
  y = approx(xs, c)

  # avoid division issues at x=0
  rel_err = (y - f) / (f + 1e-30)

  return rel_err

initial_guess=(-2,0.1)
res = least_squares(residual, initial_guess)

c = res.x
print ("# c1=%.17g" % c[0])
print ("# c2=%.17g" % c[1])

hz = 20
while hz < 0.499*44100:
  v = hz / 44100 * np.pi
  if v < np.pi / 4:
    a = approx (v, c)
  else:
    a = 1 / approx ((np.pi / 2 - v), c)
  print (hz, "%.17g %.17g %.17g" % (v, np.tan (v), a))
  hz *= 1.001

