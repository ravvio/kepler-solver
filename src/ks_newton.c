#include "ks.h"

double_t ks_newton_v0(double_t M, double_t e) {
  double_t E = M;
  double_t delta = 0.0;

  for (int i = 0; i < MAX_ITERATIONS; ++i) {
    delta = -(E - e * sin(E) - M) / (1 - e * cos(E));
    if (fabs(delta) < TOLERANCE) {
      break;
    }
    E += delta;
  }

  return E;
}

