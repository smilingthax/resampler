#include "sincwin.h"
#include <stdlib.h>
#include <math.h>

struct _sincwin_t {
  uint32_t halflen;
  float freq;
  float (*window_fn)(float pos);
};

sincwin_t *sincwin_create(uint32_t halflen, uint32_t num_phases, float freq, float (*window_fn)(float pos))
{
  // assert(halflen > 0); // num_phases not used ...
  // assert(freq > 0.0f && freq <= 1.0f);

  sincwin_t *ret = (sincwin_t *)malloc(sizeof(sincwin_t));
  if (!ret) {
    return NULL;
  }

  ret->halflen = halflen;
  ret->freq = freq;
  ret->window_fn = window_fn;

  return ret;
}

void sincwin_destroy(sincwin_t *sw)
{
  free(sw);
}

void sincwin_calc_fir(sincwin_t *sw, float *dst, float phase)
{
  // assert(sw);
  // assert(phase >= 0.0f && phase < 1.0f);
  const uint32_t halflen = sw->halflen;

  const float freq = sw->freq;
  float (*window_fn)(float pos) = sw->window_fn;
  const float atten = freq; // attenuation to prevent clipping (-> all energy above freq could also end up in output)

  for (int32_t i = 0, len = 2 * halflen; i < len; i++) {
    const float pos = (i - ((int32_t)halflen - 1) - phase) / halflen;
    const float x = pos * halflen * freq * M_PI;
    if (pos == 0.0f) {
      dst[i] = atten * 1.0f;
    } else {
      dst[i] = atten * sinf(x) / x * window_fn(pos);
    }
  }
#if 1 // fir output is shifted left by 1, for table based method, and require outermost coefficients to be 0
  // we can't rely on window_fn(1) to be zero, but calc_fir needs/uses non-symmetric half-open interval [-1,1) ...
  // -> assume window to be infinitesimally smaller: [-1+eps,1-eps] -> (-1,1)
  if (phase == 0.0f) {     // FIXME/PROBLEM: w/ num_phases, interpolation to zero for phase nearly 0... (window does help)
    dst[2 * halflen - 1] = 0.0f; // used as last coeff in phase=0
  }
#endif
}

