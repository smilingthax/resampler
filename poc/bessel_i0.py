#!/usr/bin/env python3

import numpy as np

# The special case of n=0 gives I_0(z) as the series
#  I_0(z) = sum_(k=0)^infty ((1/4 z^2)^k) / ((k!)^2) = sum_(k=0)^infty ( (z/2)^k / k! )^2 = ((z/2)^0 / 0!)^2 + ((z/2)^1 / 1!)^2 + ((z/2)^2 / 2!)^2 + ...

# (sox uses a similar implementation)
def _i0(z):
  z2 = z / 2
  term = 1
  k = 1
  sum = 1   # pre-sum first term (which needs 0! == 1)
  last_sum = 0
  while sum != last_sum:
    last_sum = sum
    term *= z2 / k
    sum += term * term
    k += 1
  return sum

i0 = np.vectorize(_i0, [float])

def kaiser(len, beta):
  return i0(beta * np.sqrt(1 - np.linspace(-1, 1, len)**2.0)) / i0(beta)

print(np.kaiser(11, 10.0))
print(kaiser(11, 10.0))

