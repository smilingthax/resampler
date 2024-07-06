#include "resampler.h"
#include <stdexcept>
#include <algorithm>
#include <string.h>

// NOTE: don't use separate sincwin_*.c impl, inlining of calc_fir improves performance!

#if __cpp_lib_math_constants >= 201907L
#include <numbers>
constexpr float pi = std::numbers::pi<float>;
#else
#include <math.h>
constexpr float pi = M_PI;
#endif

//#define SCALE_FIR_BY_FREQ

Resampler::Resampler(float ratio)
  : ipos(0), iend(0),
    iphase(0.0f)
{
  if (ratio <= 0.0f) {
    throw std::invalid_argument("ratio must be > 0"); // safe, no allocs yet!
  }
  istep = 1.0f / ratio;

  float f_filt = 0.95f;
  const uint32_t sinc_hperiods = 32; // "number of periods" / 2
  const uint32_t sinc_resolution = 256; // number of "phases" between periods (assuming freq == 1.0)

  // assert(f_filt <= 1.0f);
#ifdef SCALE_FIR_BY_FREQ
  halflen = (uint32_t)ceilf(sinc_hperiods * std::max(1.0f, istep) / f_filt); // istep = 1/ratio
  f_filt = 1.0f;
#else
  halflen = (uint32_t)ceilf(sinc_hperiods * std::max(1.0f, istep));
#endif

  if (istep > 1.0f) { // downsample
    const float rscale = (float)sinc_hperiods / halflen;
    gen_sinc_table(halflen, sinc_resolution * rscale, f_filt * rscale); // (halflen = sinc_hperiods / rscale)
  } else { // upsample
    gen_sinc_table(halflen, sinc_resolution, f_filt);
  }
  fir.reset(new float[2 * halflen]); // uninitialized ("for overwrite")
  // assert(istep < 2 * halflen); // fir.size());

  ibuf.reset(new float[4 * halflen]{}); // zero initialized
}

// pos = [0,1]  (or [-1,1])
static float window(float pos)
{
  // blackman
  return 0.42f + 0.5f * cosf(pi * pos) + 0.08f * cosf(2 * pi * pos);
  // return 0.5f * (1.0f - alpha) + 0.5f * cosf(pi * pos) + 0.5f * \alpha * cosf(2 * pi * pos);
}

void Resampler::gen_sinc_table(uint32_t fir_hlen, uint32_t num_phases, float freq)
{
  const uint32_t hlen = fir_hlen * num_phases;
  // assert(hlen > 0);
  sinc_table.reset(new float[hlen + 1]); // uninitialized ("for overwrite")

  const float atten = freq; // attenuation to prevent clipping (-> all energy above freq could also end up in output)

  float *t = sinc_table.get();
  for (uint32_t phase = num_phases - 1; phase > 0; phase--) {
    for (uint32_t j = 0; j < halflen; j++) {
      const float pos = (j * num_phases + phase) / (float)hlen;
      const float x = pos * halflen * freq * M_PI;
      *t++ = atten * sinf(x) / x * window(pos);
    }
  }
  // phase 0 is special
  *t++ = atten * 1.0f * window(0.0f); // (window_fn(0) should also be 1.0f ...)
  for (uint32_t j = 1; j < halflen; j++) {
    const float pos = (j * num_phases + 0) / (float)hlen;
    const float x = pos * halflen * freq * M_PI;
    *t++ = atten * sinf(x) / x * window(pos);
  }
  // we can't rely on window(1) to be zero, but calc_fir needs/uses non-symmetric half-open interval [-1,1) ...
  // -> assume window to be infinitesimally smaller: [-1+eps,1-eps] -> (-1,1)
  *t++ = 0.0f; // used as last coeff in phase=0
  this->num_phases = num_phases;
}

// phase = [0,1)
void Resampler::calc_fir(float phase)
{
  // assert(phase >= 0.0f && phase < 1.0f);
  const uint32_t halflen = this->halflen;
  const uint32_t num_phases = this->num_phases;

  phase *= num_phases;
  const uint32_t p = (uint32_t)phase;
  phase -= p;

  const float *s0 = sinc_table.get() + (num_phases - 1 - p) * halflen,
              *s1 = (p == num_phases - 1) ? sinc_table.get() + (num_phases - 1) * halflen + 1 : s0 - halflen;
  float *d = fir.get() + halflen - 1;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = *s0 * (1.0f - phase) + *s1 * phase;
    ++s0;
    ++s1;
    --d;
  }

  s1 = sinc_table.get() + p * halflen;
  s0 = (p == 0) ? sinc_table.get() + (num_phases - 1) * halflen + 1 : s1 - halflen;
  d = fir.get() + halflen;
  for (uint32_t i = 0; i < halflen; i++) {
    *d = *s0 * (1.0f - phase) + *s1 * phase;
    ++s0;
    ++s1;
    ++d;
  }
}

static float dot(const float *as, const float *bs, uint32_t len)
{
  // https://www.earlevel.com/main/2012/12/03/a-note-about-de-normalization/
  float ret = 1e-30; // prevent denormals
  for (uint32_t i = 0; i < len; i++) {
    // ret += (*as) * (*bs);
    ret = fmaf((*as), (*bs), ret);
    ++as;
    ++bs;
  }
  return ret - 1e-30;
}

void Resampler::process2(float *&out, uint32_t &olen, const float *&in, uint32_t &ilen)
{
  const uint32_t fir_len = 2 * halflen;
  for (; olen; olen--) {
    // consume enough input for next block (or return)
    if (ipos >= fir_len) {
      // assert(iend - ipos < fir_len); // and thus: (iend - ipos < ipos), no overlap in memcpy
      iend -= ipos; // "fill"
      memcpy(ibuf.get(), ibuf.get() + ipos, iend * sizeof(*ibuf.get()));
      ipos = 0;
    }
    const uint32_t ineed = ipos + fir_len;
    // assert(ineed < 2 * fir_len);
    for (; iend < ineed; iend++) {
      if (!ilen) {
        return;
      }
      ibuf[iend] = *in;
      ++in;
      --ilen;
    }
    // assert(iend == ineed && iend <= 2 * fir_len);

    calc_fir(iphase);

    *out = dot(ibuf.get() + ipos, fir.get(), fir_len);
    ++out;

    // advance ipos + iphase
    iphase += istep;
    const int ip = (int)iphase;
    iphase -= ip;
    ipos += ip;
    // assert(iphase >= 0.0f && iphase < 1.0f);
    // assert(ipos < iend); // because (istep < fir_len)
  }
}

