#include "ks.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TOLERANCE 3e-9
#define MAX_ITERATIONS 100

double_t solve_orbit(double_t T, double_t a, double_t e, double_t t,
                     KS_Solver solver) {
  // Mean anomaly
  double_t M = 2 * M_PI * t / T;

  // Eccentric anomaly

  clock_t start = clock();
  double_t E = solver(M, e);
  clock_t end = clock();
  float seconds = (float)(end - start) / CLOCKS_PER_SEC;
  double_t error = fabs(M - E + e * sin(E));
  fprintf(stdout, "%e %e %e %f ms\n", M, E, error, seconds * 1000);

  // True anomaly
  double_t b = e / (1.0 + sqrt(1 - e * e));
  double_t theta = E + 2 * atan((b * sin(E)) / (1 - b * cos(E)));
  // Heliocentric distance
  double_t r = a * (1 - e * cos(E));

  return E;
}

void solve_kepler_full(double_t e, char *prefix, KS_Solver solver) {
  int steps = 100000;
  double_t error_total = 0.0;

  clock_t start = clock();
  for (int i = 0; i < steps; i++) {
    double_t M = 2 * M_PI * i / steps;
    double_t E = solver(M, e);
    error_total += fabs(M - E + e * sin(E));
  }
  clock_t end = clock();
  float_t seconds = (float)(end - start) / CLOCKS_PER_SEC;
  double_t error_avg = error_total / steps;

  fprintf(stdout, "[%s] Time: %fms %frad\n", prefix, seconds, error_avg);
}

int main(int argc, char **argv) {
  char *end_ptr;
  // Period
  double_t T = strtod(argv[1], &end_ptr);
  // Semi-Major Axis
  double_t a = strtod(argv[2], &end_ptr);
  // Eccentricity
  double_t e = strtod(argv[3], &end_ptr);

  // double_t time = strtod(argv[4], &end_ptr);
  // double_t t = 0.0;

  solve_kepler_full(e, "BISEC_v0", ks_bisection_v0);
  solve_kepler_full(e, "NEWTO_v0", ks_newton_v0);
  solve_kepler_full(e, "ENRKE_v0", ks_enrke_v0);
  solve_kepler_full(e, "ENRKE_v1", ks_enrke_v1);
}
