#include "resampler.h"
//#include <stdio.h>

#define BLOCKSIZE 1024

int main(int argc, char **argv)
{
  FILE *in = stdin;
  FILE *out = stdout;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s [RATIO | OUTFREQ/INFREQ]\n", argv[0]);
    return 1;
  }
  // assert(argv[1]);

  char *end;
  double num = strtod(argv[1], &end);
  // assert(end);
  if (*end == '/') {
    double denom = strtod(end + 1, &end);
    if (!*end) {
      if (denom <= 0.0) {
        fprintf(stderr, "INFREQ must be > 0\n");
        return 2;
      }
      num /= denom;
    }
  }
  if (*end) {
    fprintf(stderr, "Error: expected floating point RATIO or OUTFREQ/INFREQ\n");
    return 1;
  }

  const float ratio = num;
  Resampler resamp(ratio);

  // read f32le, write f32le
  float inbuf[BLOCKSIZE];
  float outbuf[BLOCKSIZE];
  while (!feof(in)) {
    uint32_t ilen = fread(inbuf, sizeof(*inbuf), BLOCKSIZE, in);

    uint32_t ipos = 0;
    while (ilen > ipos) {
      uint32_t olen = BLOCKSIZE;
      ipos += resamp.process(outbuf, olen, inbuf + ipos, ilen - ipos);
      // assert(olen > 0);
      fwrite(outbuf, sizeof(*outbuf), olen, out);
    }
  }

  return 0;
}
