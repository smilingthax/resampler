#pragma once

#include <memory>

class Resampler {
public:
  Resampler(float ratio); //, uint32_t sinc_hperiods);

  void process2(float *&out, uint32_t &olen, const float *&in, uint32_t &ilen);

  // returns number of used input samples
  uint32_t process(float *out, uint32_t &olen, const float *in, uint32_t ilen) {
    float *o = out;
    const float *i = in;
    process2(o, olen, i, ilen);
    olen = o - out;  // return written, not remaining
    return i - in;
  }

private:
  void gen_sinc_table(uint32_t fir_hlen, uint32_t phases, float freq);

  // phase = [0, 1)
  void calc_fir(float phase);

private:
  float istep;

  uint32_t halflen, num_phases;
  std::unique_ptr<float[]> sinc_table; // size: halflen * num_phases + 1

  std::unique_ptr<float[]> fir; // size: 2 * halflen  // temporary array of current-phase sinc coefficients
  std::unique_ptr<float[]> ibuf; // size: 4 * halflen

  uint32_t ipos, iend;
  float iphase;
};

