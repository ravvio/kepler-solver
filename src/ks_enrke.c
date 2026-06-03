// ENRKE based on the paper: "Two fast and accurate routines for solving the
// elliptic Kepler equation for all values of the eccentricity and mean anomaly"
// by Daniele Tommasini and David N. Olivieri.
// https://www.aanda.org/articles/aa/full_html/2022/02/aa41423-21/aa41423-21.html#S18

#include "ks.h"

// Newton-Rhapson solver with a OG131 (Odell & Gooding 1986) iterator for the
// first step.
// https://www.aanda.org/articles/aa/full_html/2022/02/aa41423-21/aa41423-21.html#R30)
double_t ks_enrke_v0(double_t M, double_t e) {
  double_t E = M;

  // First step
  double_t sinE = e * sin(E);
  double_t cosE = e * cos(E);
  double_t f0 = E - sinE - M;
  double_t f1 = 1 - cosE;
  double_t f2 = sinE;
  double_t f3 = cosE;
  double_t f111 = f1 * f1 * f1;
  double_t f012 = f0 * f1 * f2;
  double_t f003 = f0 * f0 * f3;

  double_t delta = -(f0 * (f111 - (1.0 / 2.0) * f012) + (1.0 / 3.0) * f003) /
                   (f1 * (f111 - f012 + (1.0 / 2.0) * f003));
  if (fabs(delta) < TOLERANCE) {
    return E;
  }
  E += delta;

  for (int i = 0; i < MAX_ITERATIONS; ++i) {
    delta = -(E - e * sin(E) - M) / (1 - e * cos(E));
    if (fabs(delta) < TOLERANCE) {
      break;
    }
    E += delta;
  }
  return E;
}

// Newton-Rhapson solver with a OG131 (Odell & Gooding 1986) iterator for the
// first step. With an alternative stopping condition that does not require to
// calculate E_(n+1).
double_t ks_enrke_v1(double_t M, double_t e) {
  double_t E = M;

  // First step
  double_t sinE = e * sin(E);
  double_t cosE = e * cos(E);
  double_t f0 = E - sinE - M;
  double_t f1 = 1 - cosE;
  double_t f2 = sinE;
  double_t f3 = cosE;
  double_t f111 = f1 * f1 * f1;
  double_t f012 = f0 * f1 * f2;
  double_t f003 = f0 * f0 * f3;

  double_t delta = -(f0 * (f111 - (1.0 / 2.0) * f012) + (1.0 / 3.0) * f003) /
                   (f1 * (f111 - f012 + (1.0 / 2.0) * f003));
  if (fabs(delta) < TOLERANCE) {
    return E;
  }
  E += delta;

  double_t stopper = 0.0;
  for (int i = 0; i < MAX_ITERATIONS; ++i) {
    cosE = e * cos(E);

    delta = -(E - e * sin(E) - M) / (1 - cosE);
    stopper = (2.0 * (1.0 - cosE) * E) / (e + DBL_EPSILON);

    E += delta;

    if (delta * delta < stopper) {
      break;
    }
  }
  return E;
}
