#include "window.h"
#include <math.h>

#if __cpp_lib_math_constants >= 201907L
#include <numbers>
static constexpr float pi = std::numbers::pi<float>;
#else
#include <math.h>
static constexpr float pi = (float)M_PI;
#endif

// TODO? c++17 std::cyl_bessel_if(0.0f, x) ?

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

static float window_blackman(float pos, void *user)
{
  const float alpha = static_cast<window_u *>(user)->blackman.alpha;
  return 0.5f * (1.0f - alpha) + 0.5f * cosf(pi * pos) + 0.5f * alpha * cosf(2 * pi * pos);
}

window_func_t window_init_blackman(window_u &win, float alpha)
{
  win.blackman.alpha = alpha;
  return window_blackman;
}

static float window_kaiser(float pos, void *user)
{
  window_u *win = static_cast<window_u *>(user);
  return i0f(win->kaiser.beta * sqrtf(1 - pos * pos)) / win->kaiser._denom;
}

window_func_t window_init_kaiser(window_u &win, float beta)
{
  win.kaiser.beta = beta;
  win.kaiser._denom = i0f(beta);
  return window_kaiser;
}

