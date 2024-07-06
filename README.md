# Audio Resampler

Based on "Bandlimited interpolation" as described in https://ccrma.stanford.edu/~jos/resample/resample.pdf.

Features:
* Windowed-sinc FIR filter, supports both upsampling and downsampling by an arbitrary ratio.
* zlib-like API: Each call processes data until either input length or output length is exhausted.
  Each samples has to be provided only once (samples needed for overlap are stored internally).

* Written in C++11.

* TBD: Blackman vs. Kaiser window?

Usage:
```
  Resampler resam(1.5f); // upsample

  while (...more data...) {
    iused = resamp.process(out, olen, in, ilen);
    // olen returns "oused"
// or:
    resamp.process2(out, olen, in, ilen);
    // out, olen, in, ilen are all updated.
  }
```

TODO:
* Support fine-grained (relative) ratio adjustments with smooth update (`set_relative_ratio()`...)
* Expose current input -> output delay (useful for conversion between clock-domains via  adaptive resampling,
  cf. https://kokkinizita.linuxaudio.org/papers/adapt-resamp.pdf).
* Multi-channel (planar? interleaved?)...
* Support setting fractional phase-shift with `ratio = 1.0`?

IDEAs:
* Allow double (currently: float) samples?

Copyright (c) 2024 Tobias Hoffmann

License: https://opensource.org/licenses/MIT

