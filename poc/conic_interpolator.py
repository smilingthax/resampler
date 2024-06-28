#!/usr/bin/env python3

# https://zipcpu.com/dsp/2018/03/30/quadratic.html
# -> a simple quadratic interpolation (support 3) does not have a continuous impulse response,
#    but that does not mean that such a function could not exist...
# Other than a cubic (e.g. catmull-rom) spline interpolator with support 4,
# the improved quadratic interpolator (continuous, but not smooth) has a even larger support of size 5!

import sympy as sy

sy.init_printing(order='grlex')

t = sy.symbols('t', real = True)

(p_2, p_1, p0, p1, p2) = sy.symbols('p_-2 p_-1 p0 p1 p2', real = True)

p = (-28 * p0 + 16 * (p1 + p_1) - 2 * (p2 + p_2)) * t**2 / 16 + (10 * (p1 - p_1) - (p2 - p_2)) * t / 16 + p0
sy.pprint(p)

sy.pprint(p.factor([p_2, p_1, p0, p1, p2]))

pl = sy.poly(p, t)
sy.pprint(sy.horner(pl))

# -- impulse response

h = sy.Piecewise(
  (0, sy.Abs(t) >= 5/2),
  (p.subs({p_2: 0, p_1: 0, p0: 0, p1: 0, p2: 1, t: t + 2}), t < -3/2),
  (p.subs({p_2: 0, p_1: 0, p0: 0, p1: 1, p2: 0, t: t + 1}), t < -1/2),
  (p.subs({p_2: 0, p_1: 0, p0: 1, p1: 0, p2: 0}), t < 1/2),
  (p.subs({p_2: 0, p_1: 1, p0: 0, p1: 0, p2: 0, t: t - 1}), t < 3/2),
  (p.subs({p_2: 1, p_1: 0, p0: 0, p1: 0, p2: 0, t: t - 2}), t < 5/2)
)
sy.pprint(h)

sy.plot(h, (t, -3, 3))

