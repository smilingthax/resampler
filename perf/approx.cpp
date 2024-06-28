#include "sincwin.h"
#include <stdio.h>

extern "C" {
#define sincwin_t sincwin_ref_t
#define _sincwin_t _sincwin_ref_t
#define sincwin_create sincwin_ref_create
#define sincwin_destroy sincwin_ref_destroy
#define sincwin_calc_fir sincwin_ref_calc_fir
typedef struct _sincwin_t sincwin_t;
#include "sincwin_ref.c"
#undef sincwin_calc_fir
#undef sincwin_destroy
#undef sincwin_create
#undef _sincwin_t
#undef sincwin_t
}

// real window is necessary to get good approximation values...
static float window(float pos)
{
  return 0.42f + 0.5f * cosf(M_PI * pos) + 0.08f * cosf(2 * M_PI * pos); // blackman
}

static float sumdiff(const float *a, const float *b, size_t len)
{
  float ret = 0.0f;

  for (size_t i = 0; i < len; i++) {
    ret += fabsf(a[i] - b[i]);
  }

  return ret;
}

int main(int argc, char **argv)
{
#define HLEN 30
//#define PHASES 120  // esp. linear
#define PHASES 40     // quadratic stores 3 coeffs per phase! (linear computes it on-the-fly)
//#define FREQ 1.0f
#define FREQ 0.95f

  sincwin_ref_t *sr = sincwin_ref_create(HLEN, PHASES, FREQ, window);
  if (!sr) {
    fprintf(stderr, "sincwin_ref_create failed\n");
    return 1;
  }

  float fir_ref[2*HLEN];

  sincwin_t *sw = sincwin_create(HLEN, PHASES, FREQ, window);
  if (!sw) {
    fprintf(stderr, "sincwin_create failed\n");
    return 1;
  }

  float fir[2*HLEN];

  const float phase_incr = 0.001f;
  float max_sd = 0.0f, sum_sd = 0.0f;
  for (float phase = 0.0f; phase <= 1.0f - phase_incr + 1e-6; phase += phase_incr) {
    sincwin_calc_fir(sw, fir, phase);
    sincwin_ref_calc_fir(sr, fir_ref, phase);
    const float sd = sumdiff(fir, fir_ref, 2*HLEN);
//    printf("%f: %g\n", phase, sd);
    if (max_sd < sd) {
      max_sd = sd;
    }
    sum_sd += sd;
  }
  const float max_rd = max_sd * phase_incr;
  printf("max diff: %g = %g%% (%g dB)\n", max_sd, 100.0f * max_rd, 20.0f * log10f(max_rd));
  printf("avg diff: %g (%g dB)\n", sum_sd * phase_incr, 20.0f * log10f(sum_sd * phase_incr));

  sincwin_destroy(sw);
  sincwin_ref_destroy(sr);

  return 0;
}

