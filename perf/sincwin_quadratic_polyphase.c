#include "sincwin.h"
#include <stdlib.h>
#include <math.h>

//#define INTERLEAVED_TABLE

#ifdef INTERLEAVED_TABLE
struct qcoff_t {
  float a, b, c;
};
#endif

static void calc_qcoeff(float *ret_a, float *ret_b, float x_2, float x_1, float x0, float x1, float x2)
{
  *ret_a = (-14.0f * x0 + 8.0f * (x1 + x_1) - x2 - x_2) / 8.0f;
  *ret_b = (10.0f * (x1 - x_1) - x2 + x_2) / 16.0f;
//  *ret_c = x0;
}

struct _sincwin_t {
  uint32_t halflen, num_phases;
#ifdef INTERLEAVED_TABLE
  struct qcoff_t table[];
#else
  float table[]; // can't have three flexible array members...
#endif
};

// FIXME...: num_phases >= 4
sincwin_t *sincwin_create(uint32_t halflen, uint32_t num_phases, float freq, float (*window_fn)(float pos, void *user), void *user)
{
  const uint32_t hlen = halflen * num_phases;
  // assert(hlen > 0);
  // assert(freq > 0.0f && freq <= 1.0f);

#ifdef INTERLEAVED_TABLE
  sincwin_t *ret = malloc(sizeof(sincwin_t) + sizeof(ret->table[0]) * (hlen + 1));
#else
  sincwin_t *ret = malloc(sizeof(sincwin_t) + 3 * sizeof(ret->table[0]) * (hlen + 1));
#endif
  if (!ret) {
    return NULL;
  }

  ret->halflen = halflen;
  ret->num_phases = num_phases;

  const float atten = freq; // attenuation to prevent clipping (-> all energy above freq could also end up in output)

  // Trick: c is always the direct value, compute these first
#ifdef INTERLEAVED_TABLE
  struct qcoff_t *t = ret->table;

  // calculate a
  for (uint32_t phase = num_phases - 1; phase > 0; phase--) {
    for (uint32_t j = 0; j < halflen; j++) {
      const float pos = (j * num_phases + phase) / (float)hlen;
      const float x = pos * halflen * freq * M_PI;
      t->c = atten * sinf(x) / x * window_fn(pos, user);
      ++t;
    }
  }
  // phase 0 is special
  t->c = atten * 1.0f * window_fn(0.0f, user); // (window_fn(0) should also be 1.0f ...)
  ++t;
  for (uint32_t j = 1; j < halflen; j++) {
    const float pos = (j * num_phases + 0) / (float)hlen;
    const float x = pos * halflen * freq * M_PI;
    t->c = atten * sinf(x) / x * window_fn(pos, user);
    ++t;
  }
  // we can't rely on window_fn(1) to be zero, but calc_fir needs/uses non-symmetric half-open interval [-1,1) ...
  // -> assume window to be infinitesimally smaller: [-1+eps,1-eps] -> (-1,1)
  t->c = 0.0f; // used as last coeff in phase=0
  //++t;

  // calculate b, c
  t = ret->table;
  const int32_t wrap = -(num_phases - 1) * halflen - 1;

  // phase (num_phase - 1)
  for (uint32_t j = 1; j < halflen; j++) {  // (actually: j = 0; j < halflen - 1)
    calc_qcoeff(&t->a, &t->b, t[2 * halflen].c, t[halflen].c, t[0].c, t[-wrap].c, t[-wrap - halflen].c);
    ++t;
  }
  calc_qcoeff(&t->a, &t->b, t[2 * halflen].c, t[halflen].c, t[0].c, t[-wrap].c, 0.0f); // [or: 0.0f, 0.0f); - last coeff was explicitly set to 0.0f ]
  ++t;

  // phase (num_phase - 2)
  for (uint32_t j = 0; j < halflen; j++) {
    calc_qcoeff(&t->a, &t->b, t[2 * halflen].c, t[halflen].c, t[0].c, t[-(int32_t)halflen].c, t[-wrap - halflen].c);
    ++t;
  }

  for (uint32_t phase = num_phases - 3; phase > 1; phase--) {
    for (uint32_t j = 0; j < halflen; j++) {
      calc_qcoeff(&t->a, &t->b, t[2 * halflen].c, t[halflen].c, t[0].c, t[-(int32_t)halflen].c, t[-2 * (int32_t)halflen].c);
      ++t;
    }
  }

  // phase 1
  calc_qcoeff(&t->a, &t->b, t[0].c, t[halflen].c, t[0].c, t[-(int32_t)halflen].c, t[-2 * (int32_t)halflen].c); // symmetric around table[0] !
  ++t;
  for (uint32_t j = 1; j < halflen; j++) {
    calc_qcoeff(&t->a, &t->b, t[wrap + (int32_t)halflen].c, t[halflen].c, t[0].c, t[-(int32_t)halflen].c, t[-2 * (int32_t)halflen].c);
    ++t;
  }

  // phase 0
  calc_qcoeff(&t->a, &t->b, t[-2 * (int32_t)halflen].c, t[-(int32_t)halflen].c, t[0].c, t[-(int32_t)halflen].c, t[-2 * (int32_t)halflen].c);
  ++t;
  for (uint32_t j = 1; j < halflen; j++) {
    calc_qcoeff(&t->a, &t->b, t[wrap + (int32_t)halflen].c, t[wrap].c, t[0].c, t[-(int32_t)halflen].c, t[-2 * (int32_t)halflen].c);
    ++t;
  }
  calc_qcoeff(&t->a, &t->b, t[wrap + (int32_t)halflen].c, t[wrap].c, t[0].c, 0.0f, 0.0f);
  //++t;

#else
  float *as = ret->table,
        *bs = ret->table + (hlen + 1),
        *cs = ret->table + 2 * (hlen + 1);

  // calculate a
  for (uint32_t phase = num_phases - 1; phase > 0; phase--) {
    for (uint32_t j = 0; j < halflen; j++) {
      const float pos = (j * num_phases + phase) / (float)hlen;
      const float x = pos * halflen * freq * M_PI;
      *cs++ = atten * sinf(x) / x * window_fn(pos, user);
    }
  }
  // phase 0 is special
  *cs++ = atten * 1.0f * window_fn(0.0f, user); // (window_fn(0) should also be 1.0f ...)
  for (uint32_t j = 1; j < halflen; j++) {
    const float pos = (j * num_phases + 0) / (float)hlen;
    const float x = pos * halflen * freq * M_PI;
    *cs++ = atten * sinf(x) / x * window_fn(pos, user);
  }
  // we can't rely on window_fn(1) to be zero, but calc_fir needs/uses non-symmetric half-open interval [-1,1) ...
  // -> assume window to be infinitesimally smaller: [-1+eps,1-eps] -> (-1,1)
  *cs++ = 0.0f; // used as last coeff in phase=0

  // calculate b, c
  cs = ret->table + 2 * (hlen + 1);
  const int32_t wrap = -(num_phases - 1) * halflen - 1;

  // phase (num_phase - 1)
  for (uint32_t j = 1; j < halflen; j++) {  // (actually: j = 0; j < halflen - 1)
    calc_qcoeff(as, bs, cs[2 * halflen], cs[halflen], cs[0], cs[-wrap], cs[-wrap - halflen]);
    ++as; ++bs; ++cs;
  }
  calc_qcoeff(as, bs, cs[2 * halflen], cs[halflen], cs[0], cs[-wrap], 0.0f); // [or: 0.0f, 0.0f); - last coeff was explicitly set to 0.0f ]
  ++as; ++bs; ++cs;

  // phase (num_phase - 2)
  for (uint32_t j = 0; j < halflen; j++) {
    calc_qcoeff(as, bs, cs[2 * halflen], cs[halflen], cs[0], cs[-(int32_t)halflen], cs[-wrap - halflen]);
    ++as; ++bs; ++cs;
  }

  for (uint32_t phase = num_phases - 3; phase > 1; phase--) {
    for (uint32_t j = 0; j < halflen; j++) {
      calc_qcoeff(as, bs, cs[2 * halflen], cs[halflen], cs[0], cs[-(int32_t)halflen], cs[-2 * (int32_t)halflen]);
      ++as; ++bs; ++cs;
    }
  }

  // phase 1
  calc_qcoeff(as, bs, cs[0], cs[halflen], cs[0], cs[-(int32_t)halflen], cs[-2 * (int32_t)halflen]); // symmetric around table[0] !
  ++as; ++bs; ++cs;
  for (uint32_t j = 1; j < halflen; j++) {
    calc_qcoeff(as, bs, cs[wrap + (int32_t)halflen], cs[halflen], cs[0], cs[-(int32_t)halflen], cs[-2 * (int32_t)halflen]);
    ++as; ++bs; ++cs;
  }

  // phase 0
  calc_qcoeff(as, bs, cs[-2 * (int32_t)halflen], cs[-(int32_t)halflen], cs[0], cs[-(int32_t)halflen], cs[-2 * (int32_t)halflen]);
  ++as; ++bs; ++cs;
  for (uint32_t j = 1; j < halflen; j++) {
    calc_qcoeff(as, bs, cs[wrap + (int32_t)halflen], cs[wrap], cs[0], cs[-(int32_t)halflen], cs[-2 * (int32_t)halflen]);
    ++as; ++bs; ++cs;
  }
  calc_qcoeff(as, bs, cs[wrap + (int32_t)halflen], cs[wrap], cs[0], 0.0f, 0.0f);
  // ++as; ++bs; ++cs;
#endif

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

  phase *= sw->num_phases;
  const uint32_t p = (uint32_t)phase;
  phase -= p;

#ifdef INTERLEAVED_TABLE
  const struct qcoff_t *s = sw->table + (sw->num_phases - 1 - p) * halflen;
  float *d = dst + halflen - 1;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = phase * (phase * s->a + s->b) + s->c;
    ++s;
    --d;
  }

  s = sw->table + ((p == 0) ? (sw->num_phases - 1) * halflen + 1 : (p - 1) * halflen);
  d = dst + halflen;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = phase * (phase * s->a - s->b) + s->c;
    ++s;
    ++d;
  }

#else
  const uint32_t hlen1 = halflen * sw->num_phases + 1;
  const float *as = sw->table + (sw->num_phases - 1 - p) * halflen,
              *bs = as + hlen1,
              *cs = as + 2 * hlen1;

  float *d = dst + halflen - 1;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = phase * (phase * *as + *bs) + *cs;
    ++as;
    ++bs;
    ++cs;
    --d;
  }

  as = sw->table + ((p == 0) ? (sw->num_phases - 1) * halflen + 1 : (p - 1) * halflen);
  bs = as + hlen1;
  cs = as + 2 * hlen1;
  d = dst + halflen;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = phase * (phase * *as - *bs) + *cs;
    ++as;
    ++bs;
    ++cs;
    ++d;
  }
#endif
}

