#include "ks.h"

double_t ks_bisection_v1(double_t M, double_t e) {
  double_t E_min = M - e;
  double_t E_max = M + e;

  double_t err_min = M - E_min + e * sin(E_min);
  double_t err_max = M - E_max + e * sin(E_max);

  // Swap if f(E_min) > f(E_max)
  if (err_min > 0) {
    double_t tmp = E_min;
    E_min = E_max;
    E_max = tmp;
  }

  double_t E;
  double_t error = MAXFLOAT;

  for (int i = 0; i < MAX_ITERATIONS; ++i) {
    E = (E_max + E_min) / 2.0;
    error = M - E + e * sin(E);

    if (fabs(error) < TOLERANCE) {
      break;
    }

    if (error > 0) {
      E_max = E;
    } else {
      E_min = E;
    }
  }
  return E;
}
