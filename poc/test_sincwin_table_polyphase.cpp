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

  void dump_table() const;
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
  table[(num_phases - 1) * halflen] = 1.0f * window(0.0f); // (window(0) should also be 1.0f ...)
#else
  table[(num_phases - 1) * halflen] = 0.0f;
#endif
#if 0
  for (uint32_t i = 1, phase = 1; i < hlen; i++, phase = (phase + 1) % num_phases) {
#ifndef DEBUG_IOTA
    const float pos = (float)i / hlen;
    const float x = pos * halflen * freq * pi;
    table[i / num_phases + (num_phases - 1 - phase) * halflen] = sinf(x) / x * window(pos);
#else
    table[i / num_phases + (num_phases - 1 - phase) * halflen] = i;
#endif
  }
#elif 0
  for (uint32_t j = 0, phase = 1, i = 1; j < halflen; j++) {
    for (; phase < num_phases; phase++, i++) {
#ifndef DEBUG_IOTA
      // const float pos = ("i" - ((int32_t)halflen - 1) - phase) * num_phases / (2 * halflen);
      const float pos = (float)i / hlen;  // i = j * num_phases + phase
      const float x = pos * halflen * freq * pi;
      table[j + (num_phases - 1 - phase) * halflen] = sinf(x) / x * window(pos);
#else
      table[j + (num_phases - 1 - phase) * halflen] = i;
#endif
    }
    phase = 0;
  }
#else
#ifdef DEBUG_IOTA
#error not implemented yet
#endif
  for (int32_t phase = num_phases - 1, i = 0; phase >= 0; phase--) {
    for (uint32_t j = 0; j < halflen; j++, i++) {
//      const float pos = ((int32_t)j + phase / (float)num_phases) / halflen;
      const float pos = (j * num_phases + phase) / (float)hlen;
      const float x = pos * halflen * freq * pi;
      table[i] = sinf(x) / x * window(pos);
//table[i] = pos;  // debug
    }
  }
  // TODO? remove from loop (unroll) ?
  table[(num_phases - 1) * halflen] = 1.0f * window(0.0f); // (window(0) should also be 1.0f ...)
#endif

#ifndef DEBUG_IOTA
  // we can't rely on window(1) to be zero, but calc_fir needs/uses non-symmetric half-open interval [-1,1) ...
  // -> assume window to be infinitesimally smaller: [-1+eps,1-eps] -> (-1,1)
  table[hlen] = 0.0f; // used as last coeff in phase=0, and for interpolation
  // For interpolation, the "phase after num_phase - 1" is needed. But this is just "phase 0" + 1 full coeff
  // - but this requires an extra coeff table entry for phase 0 (-> hlen).
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

  const float *s0, *s1;
  float *d;
#if 1
  s0 = table.get() + (num_phases - 1 - p) * halflen;
  s1 = (p == num_phases - 1) ? table.get() + (num_phases - 1) * halflen + 1 : s0 - halflen;
  d = dst + halflen - 1;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = *s0 * (1.0f - phase) + *s1 * phase;
    ++s0;
    ++s1;
    --d;
  }
#endif

  s1 = table.get() + p * halflen;
  s0 = (p == 0) ? table.get() + (num_phases - 1) * halflen + 1 : s1 - halflen;
  d = dst + halflen;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = *s0 * (1.0f - phase) + *s1 * phase;
    ++s0;
    ++s1;
    ++d;
  }
}

void SincWinTable::dump_table() const
{
  for (uint32_t j = 0, k = 0; j < num_phases; j++) {
    printf("%d: ", j);
    for (uint32_t i = 0; i < halflen + (j == num_phases - 1 ? 1 : 0); i++, k++) {
      printf("%.2f ", table[k]);
    }
    printf("\n");
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

//tbl.dump_table();
// return 0;

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
