#include "kaiser.h"
#include <math.h>

// The special case of n=0 gives I_0(z) as the series
// I_0(z) = sum_(k=0)^infty ((1/4 z^2)^k) / ((k!)^2) = sum_(k=0)^infty ( (z/2)^k / k! )^2 = ((z/2)^0 / 0!)^2 + ((z/2)^1 / 1!)^2 + ((z/2)^2 / 2!)^2 + ...
static float i0f(float x) {
  const float x2 = x / 2.0f;
  // pre-sum first term (which needs 0! == 1)
  float term = 1.0f,
        sum = 1.0f,
        last_sum = 0.0f;
  for (int k = 1; sum != last_sum; k++) {
    last_sum = sum;
    term *= x2 / k;
    sum += term * term;
  }
  return sum;
}

// pos = [-1, 1]
float kaiserf(float pos, float beta)
{
  return i0f(beta * sqrtf(1 - pos * pos)) / i0f(beta);
}

