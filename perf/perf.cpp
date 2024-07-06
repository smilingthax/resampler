#include <chrono>
#include <stdio.h>
#include "sincwin.h"

static float window(float pos)
{
  return 1.0f;   // FIXME
}

int main()
{
#define HLEN 32
//#define PHASES 256
#define PHASES 40
//#define FREQ 1.0f
#define FREQ 0.95f

  sincwin_t *sw = sincwin_create(HLEN, PHASES, FREQ, window);
  if (!sw) {
    fprintf(stderr, "sincwin_create failed\n");
    return 1;
  }

  float fir[2*HLEN];

  const auto t0 = std::chrono::high_resolution_clock::now();

  for (size_t i=0; i < 10000000; i++) {
    sincwin_calc_fir(sw, fir, 0.1f);
    asm volatile("" : : "g" (fir) );
  }

  const auto t1 = std::chrono::high_resolution_clock::now();
  const auto dt1 = t1 - t0;

  printf("%.3f ms\n", dt1 / std::chrono::microseconds(1) / 1000.f);

  sincwin_destroy(sw);

  return 0;
}
