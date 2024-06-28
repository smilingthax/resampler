#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _sincwin_t sincwin_t;

sincwin_t *sincwin_create(uint32_t halflen, uint32_t num_phases, float freq, float (*window_fn)(float pos));
void sincwin_destroy(sincwin_t *sw);

// dst must have (2 * halflen) space
// phase = [0, 1)
void sincwin_calc_fir(sincwin_t *sw, float *dst, float phase);

#ifdef __cplusplus
}
#endif
