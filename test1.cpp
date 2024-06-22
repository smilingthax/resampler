#include "resampler.h"

int main(int argc, char **argv)
{
  Resampler resamp(1.0f);

  float in[128] = {}, out[128];
in[31+4] = 1000.0f;

  uint32_t ilen=63 + 9, olen=10;
  resamp.process(out, olen, in, ilen);
  printf("olen: %d, ilen %d\n", olen, ilen);

  for (size_t i = 0; i < 10; i++) printf("%f ", out[i]); printf("\n");

  return 0;
}
