#include "sincwin.h"
#include <stdlib.h>
#include <math.h>

#define INTERLEAVED_TABLE

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

// FIXME...: hlen >= 2
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
  t->c = atten * 1.0f * window_fn(0.0f, user); // (window_fn(0) should also be 1.0f ...)
  ++t;
  for (uint32_t i = 1; i < hlen; i++) {
    const float pos = (float)i / hlen;
    const float x = pos * halflen * freq * M_PI;
    t->c = atten * sinf(x) / x * window_fn(pos, user);
    ++t;
  }
  // we can't rely on window_fn(1) to be zero, but calc_fir needs/uses non-symmetric half-open interval [-1,1) ...
  // -> assume window to be infinitesimally smaller: [-1+eps,1-eps] -> (-1,1)
  t->c = 0.0f; // used as last coeff in phase=0
  // ++t;

  t = ret->table;
  calc_qcoeff(&t->a, &t->b, t[2].c, t[1].c, t[0].c, t[1].c, t[2].c); // symmetric around table[0] !
  ++t;
  calc_qcoeff(&t->a, &t->b, t[0].c, t[-1].c, t[0].c, t[1].c, t[2].c);  // table[0] = t[-1]
  ++t;
  for (uint32_t i = 3; i < hlen; i++) { // (actually: i = 2; i < hlen - 1)
    calc_qcoeff(&t->a, &t->b, t[-2].c, t[-1].c, t[0].c, t[1].c, t[2].c);
    ++t;
  }
  calc_qcoeff(&t->a, &t->b, t[-2].c, t[-1].c, t[0].c, t[1].c, 0.0f);
  ++t;
  calc_qcoeff(&t->a, &t->b, t[-2].c, t[-1].c, t[0].c, 0.0f, 0.0f);
  //++t;

#else
  float *as = ret->table,
        *bs = ret->table + (hlen + 1),
        *cs = ret->table + 2 * (hlen + 1);

  *cs++ = atten * 1.0f * window_fn(0.0f, user); // (window_fn(0) should also be 1.0f ...)
  for (uint32_t i = 1; i < hlen; i++) {
    const float pos = (float)i / hlen;
    const float x = pos * halflen * freq * M_PI;
    *cs++ = atten * sinf(x) / x * window_fn(pos, user);
  }
  // we can't rely on window_fn(1) to be zero, but calc_fir needs/uses non-symmetric half-open interval [-1,1) ...
  // -> assume window to be infinitesimally smaller: [-1+eps,1-eps] -> (-1,1)
  *cs++ = 0.0f; // used as last coeff in phase=0

  cs = ret->table + 2 * (hlen + 1);
  calc_qcoeff(as, bs, cs[2], cs[1], cs[0], cs[1], cs[2]); // symmetric around table[0] !
  ++as; ++bs; ++cs;
  calc_qcoeff(as, bs, cs[0], cs[-1], cs[0], cs[1], cs[2]);  // table[0].cs = cs[-1]
  ++as; ++bs; ++cs;
  for (uint32_t i = 3; i < hlen; i++) { // (actually: i = 2; i < hlen - 1)
    calc_qcoeff(as, bs, cs[-2], cs[-1], cs[0], cs[1], cs[2]);
    ++as; ++bs; ++cs;
  }
  calc_qcoeff(as, bs, cs[-2], cs[-1], cs[0], cs[1], 0.0f);
  ++as; ++bs; ++cs;
  calc_qcoeff(as, bs, cs[-2], cs[-1], cs[0], 0.0f, 0.0f);
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
  const uint32_t num_phases = sw->num_phases;

  phase *= num_phases;
  const uint32_t p = (uint32_t)phase;
  phase -= p;

#ifdef INTERLEAVED_TABLE
  const struct qcoff_t *s = sw->table + p;
  float *d = dst + halflen - 1;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = phase * (phase * s->a + s->b) + s->c;
    s += num_phases;
    --d;
  }

  s = sw->table + (num_phases - p);
  d = dst + halflen;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = phase * (phase * s->a - s->b) + s->c;
    s += num_phases;
    ++d;
  }

#else
  const uint32_t hlen1 = halflen * num_phases + 1;
  const float *as = sw->table + p,
              *bs = sw->table + hlen1 + p,
              *cs = sw->table + 2 * hlen1 + p;

  float *d = dst + halflen - 1;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = phase * (phase * *as + *bs) + *cs;
    as += num_phases;
    bs += num_phases;
    cs += num_phases;
    --d;
  }

  const uint32_t q = num_phases - p;
  as = sw->table + q,
  bs = sw->table + hlen1 + q,
  cs = sw->table + 2 * hlen1 + q;
  d = dst + halflen;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = phase * (phase * *as - *bs) + *cs;
    as += num_phases;
    bs += num_phases;
    cs += num_phases;
    ++d;
  }
#endif
}

