#include "sincwin.h"
#include <stdlib.h>
#include <math.h>

struct _sincwin_t {
  uint32_t halflen, num_phases;
  float table[];
};

sincwin_t *sincwin_create(uint32_t halflen, uint32_t num_phases, float freq, float (*window_fn)(float pos))
{
  const uint32_t hlen = halflen * num_phases;
  // assert(hlen > 0);
  // assert(freq > 0.0f && freq <= 1.0f);

  sincwin_t *ret = malloc(sizeof(sincwin_t) + sizeof(ret->table[0]) * (hlen + 1));
  if (!ret) {
    return NULL;
  }

  ret->halflen = halflen;
  ret->num_phases = num_phases;

  const float atten = freq; // attenuation to prevent clipping (-> all energy above freq could also end up in output)

  float *t = ret->table;
  *t++ = atten * 1.0f * window_fn(0.0f); // (window_fn(0) should also be 1.0f ...)
  for (uint32_t i = 1; i < hlen; i++) {
    const float pos = (float)i / hlen;
    const float x = pos * halflen * freq * M_PI;
    *t++ = atten * sinf(x) / x * window_fn(pos);
  }
  // we can't rely on window_fn(1) to be zero, but calc_fir needs/uses non-symmetric half-open interval [-1,1) ...
  // -> assume window to be infinitesimally smaller: [-1+eps,1-eps] -> (-1,1)
  *t++ = 0.0f; // used as last coeff in phase=0, and for interpolation

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
  const uint32_t num_phases = sw->num_phases;

  phase *= num_phases;
  const uint32_t p = (uint32_t)phase;
  phase -= p;

  const float *s = sw->table + p;
  float *d = dst + halflen - 1;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = s[0] * (1.0f - phase) + s[1] * phase;
    s += num_phases;
    --d;
  }

  s = sw->table + (num_phases - 1 - p);
  d = dst + halflen;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = s[1] * (1.0f - phase) + s[0] * phase;
    s += num_phases;
    ++d;
  }
}

