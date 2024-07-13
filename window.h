#pragma once

union window_u {
  struct {
    float alpha;
  } blackman;
  struct {
    float beta;
    float _denom;
  } kaiser;
};

// pos = [0,1]  (or [-1,1])
typedef float (*window_func_t)(float pos, void *);

window_func_t window_init_blackman(window_u &win, float alpha=0.16f);

window_func_t window_init_kaiser(window_u &win, float beta=12.0f);

