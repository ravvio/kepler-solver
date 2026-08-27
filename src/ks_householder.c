#include "./ks.h"
#include <math.h>

// Newton-Rhapson delta
double_t delta_newton(double_t M, double_t e, double_t E) {
  return -(E - e * sin(E) - M) / (1 - e * cos(E));
}

// Halley delta (OG12)
double_t delta_halley(double_t M, double_t e, double_t E) {
  double_t sinE = e * sin(E);
  double_t cosE = e * cos(E);
  double_t f0 = E - sinE - M;
  double_t f1 = 1 - cosE;
  double_t f2 = sinE;

  return -(f0 * f1) / (f1 * f1 - (1.0 / 2.0) * f0 * f2);
}

// Odell-Gooding 131 delta
double_t delta_og131(double_t M, double_t e, double_t E) {
  double_t sinE = e * sin(E);
  double_t cosE = e * cos(E);
  double_t f0 = E - sinE - M;
  double_t f1 = 1 - cosE;
  double_t f2 = sinE;
  double_t f3 = cosE;
  double_t f111 = f1 * f1 * f1;
  double_t f012 = f0 * f1 * f2;
  double_t f003 = f0 * f0 * f3;

  return -(f0 * (f111 - (1.0 / 2.0) * f012) + (1.0 / 3.0) * f003) /
         (f1 * (f111 - f012 + (1.0 / 2.0) * f003));
}

// Prussing-Conway seed
double_t seed_pc(double_t M, double_t e) { return M + 0.5 * e; }

// Householder method with a seed given by SEED first delta given by FN0 and
// successive delta by FN
#define KS_HH(NAME, SEED, FN0, FN)                                             \
  double_t NAME(double_t M, double_t e) {                                      \
    double_t E = SEED(M, e);                                                   \
    double_t delta = FN0(M, e, E);                                             \
    if (fabs(delta) < TOLERANCE) {                                             \
      return E;                                                                \
    }                                                                          \
    E += delta;                                                                \
    for (int i = 0; i < MAX_ITERATIONS; ++i) {                                 \
      delta = FN(M, e, E);                                                     \
      if (fabs(delta) < TOLERANCE) {                                           \
        break;                                                                 \
      }                                                                        \
      E += delta;                                                              \
    }                                                                          \
    return E;                                                                  \
  }

KS_HH(ks_hh_pc_nr_nr, seed_pc, delta_newton, delta_newton);
KS_HH(ks_hh_pc_ha_nr, seed_pc, delta_halley, delta_newton);
KS_HH(ks_hh_pc_og131_nr, seed_pc, delta_og131, delta_newton);

// Enhanced Householder method with a seed given by SEED first step given by FN0
// and an alternative stopping condition that lets us skip a delta calculation.
// Based on the paper: "Two fast and accurate routines for solving the
// elliptic Kepler equation for all values of the eccentricity and mean anomaly"
// by Daniele Tommasini and David N. Olivieri.
// https://www.aanda.org/articles/aa/full_html/2022/02/aa41423-21/aa41423-21.html#S18
#define KS_EHH(NAME, SEED, FN0)                                                \
  double_t NAME(double_t M, double_t e) {                                      \
    double_t E = SEED(M, e);                                                   \
    double_t delta = FN0(M, e, E);                                             \
    if (fabs(delta) < TOLERANCE) {                                             \
      return E;                                                                \
    }                                                                          \
    E += delta;                                                                \
    for (int i = 0; i < MAX_ITERATIONS; ++i) {                                 \
      double_t cosE = e * cos(E);                                              \
      double_t stop = 2 * (1 - cosE) * TOLERANCE / (e + DBL_EPSILON);          \
      delta = -(E - e * sin(E) - M) / (1 - cosE);                              \
      E += delta;                                                              \
      if (delta * delta < stop) {                                              \
        break;                                                                 \
      }                                                                        \
    }                                                                          \
    return E;                                                                  \
  }

KS_EHH(ks_ehh_pc_nr, seed_pc, delta_newton);
KS_EHH(ks_ehh_pc_ha, seed_pc, delta_halley);
KS_EHH(ks_ehh_pc_og131, seed_pc, delta_og131);
