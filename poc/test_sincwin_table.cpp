#include <math.h>

#include <stdint.h>
#include <memory>

//#define DEBUG_IOTA

#include <math.h>
constexpr float pi [[maybe_unused]] = M_PI;

static float window(float pos)
{
  return 1.0f;
}

void calc_fir_ref(float *dst, float phase)
{
const uint32_t num_phases = 10;
const uint32_t halflen = 3;
const float freq = 1.0f;

  for (int32_t i = 0, len = 2 * halflen; i < len; i++) {
#ifndef DEBUG_IOTA
    const float pos = (i - ((int32_t)halflen - 1) - phase) * num_phases / (2 * halflen);
    const float x = pos * halflen * freq * pi;
    if (pos == 0.0f) {
      dst[i] = 1.0f;
    } else {
      dst[i] = sinf(x) / x * window(pos);
    }
#else
    dst[i] = fabsf(i - ((int32_t)halflen - 1) - phase) * num_phases;
#endif
  }
}

// ---

class SincWinTable {
public:
  void gen_sinc_table(uint32_t halflen, uint32_t num_phases, float freq);

  // dst must have (2 * halflen) space
  // phase = [0, 1)
  void calc_fir(float *dst, float phase) const;

private:
  std::unique_ptr<float[]> table;
  uint32_t halflen, num_phases;
};

// ---

void SincWinTable::gen_sinc_table(uint32_t halflen, uint32_t num_phases, float freq)
{
  const uint32_t hlen = halflen * num_phases;
  // assert(hlen > 0);
  table.reset(new float[hlen + 1]); // uninitialized ("for overwrite")
#ifndef DEBUG_IOTA
  table[0] = 1.0f * window(0.0f); // (window(0) should also be 1.0f ...)
#else
  table[0] = 0.0f;
#endif
  for (uint32_t i = 1; i < hlen; i++) {
#ifndef DEBUG_IOTA
    const float pos = (float)i / hlen;
    const float x = pos * halflen * freq * pi;
    table[i] = sinf(x) / x * window(pos);
#else
    table[i] = i;
#endif
  }
#ifndef DEBUG_IOTA
  // we can't rely on window(1) to be zero, but calc_fir needs/uses non-symmetric half-open interval [-1,1) ...
  // -> assume window to be infinitesimally smaller: [-1+eps,1-eps] -> (-1,1)
  table[hlen] = 0.0f; // used as last coeff in phase=0, and for interpolation
#else
  table[hlen] = hlen;
#endif

  this->halflen = halflen;
  this->num_phases = num_phases;
}

void SincWinTable::calc_fir(float *dst, float phase) const
{
  // assert(phase >= 0.0f && phase < 1.0f);
  phase *= num_phases;
  const uint32_t p = (uint32_t)phase;
  phase -= p;
//printf("%d . %f\n", p, phase);

  const float *s;
  float *d;
#if 1
  s = table.get() + p;
  d = dst + halflen - 1;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = s[0] * (1.0f - phase) + s[1] * phase;
    s += num_phases;
    --d;
  }
#endif

  s = table.get() + (num_phases - 1 - p);
  d = dst + halflen;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = s[1] * (1.0f - phase) + s[0] * phase;
    s += num_phases;
    ++d;
  }
}

static void print_fir(const float *fir, size_t len)
{
  for (size_t i = 0; i < len; i++) {
    printf("%.2f ", fir[i]);
  }
  printf("\n");
}

#define calc_fir  tbl.calc_fir
//#define calc_fir  calc_fir_ref

int main()
{
#define FIR_HLEN  3
  SincWinTable tbl;
  float fir[2 * FIR_HLEN];

  tbl.gen_sinc_table(FIR_HLEN, 10, 1.0f);

  printf("Phase 0:\n");
  calc_fir(fir, 0.0f);
  print_fir(fir, sizeof(fir)/sizeof(*fir));

  printf("Phase 0.01:\n");
  calc_fir(fir, 0.01f);
  print_fir(fir, sizeof(fir)/sizeof(*fir));

  printf("Phase 0.1:\n");
  calc_fir(fir, 0.1f);
  print_fir(fir, sizeof(fir)/sizeof(*fir));

  printf("Phase 0.25:\n");
  calc_fir(fir, 0.25f);
  print_fir(fir, sizeof(fir)/sizeof(*fir));

  printf("Phase 0.9:\n");
  calc_fir(fir, 0.9f);
  print_fir(fir, sizeof(fir)/sizeof(*fir));

  printf("Phase 0.99:\n");
  calc_fir(fir, 0.99f);
  print_fir(fir, sizeof(fir)/sizeof(*fir));

  return 0;
}
