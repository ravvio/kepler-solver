#include "ks.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

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

void report_full_orbits(int iterations, char *prefix, KS_Solver solver) {
  char *filename_report;
  asprintf(&filename_report, "./output/report_%s.csv", prefix);
  fprintf(stderr, "[%-10s] Report file: %s\n", prefix, filename_report);

  char *filename_orbit;
  asprintf(&filename_orbit, "./output/orbit_%s.csv", prefix);
  fprintf(stderr, "[%-10s] Orbit file: %s\n", prefix, filename_orbit);

  FILE *file_report;
  file_report = fopen(filename_report, "w");
  fprintf(file_report,
          "e,time_avg_ms,error_avg_rad,error_min_rad,error_max_rad\n");

  FILE *file_orbit;
  file_orbit = fopen(filename_orbit, "w");
  fprintf(file_orbit, "e,M,time_ms,error_rad\n");

  // Solve the full orbit n (= iterations) times in m (= steps) steps
  // for many values of eccentricity
  int steps = 1000;
  for (double_t e = 0.01; e < 0.99; e += 0.01) {
    double_t error_sum_rad = 0.0;
    double_t error_min_rad = MAXFLOAT;
    double_t error_max_rad = 0.0;
    double_t time_sum_ms = 0.0;

    for (int step = 0; step < steps; step++) {
      double_t E, M, m;
      M = 2 * M_PI * step / steps;
      if (M > M_PI) {
        m = 2 * M_PI - M;
      } else {
        m = M;
      }

      clock_t start = clock();
      for (int i = 0; i < iterations; i++) {
        E = solver(m, e);
      }
      if (M > M_PI) {
        E = 2 * M_PI - E;
      }
      clock_t end = clock();

      double_t err = fabs(M - E + e * sin(E));
      error_sum_rad += err;
      error_min_rad = fmin(error_min_rad, err);
      error_max_rad = fmax(error_max_rad, err);

      float_t ms = 1000 * (float)(end - start) / CLOCKS_PER_SEC;
      time_sum_ms += ms / iterations;

      fprintf(file_orbit, "%e,%e,%e,%e\n", e, M, err, ms);
    }

    float_t time_avg_ms = time_sum_ms / steps;
    double_t error_avg_rad = error_sum_rad / steps;
    fprintf(file_report, "%e,%e,%e,%e,%e\n", e, time_avg_ms, error_avg_rad,
            error_min_rad, error_max_rad);
    fprintf(stdout, "[%-10s] e: %e | Time: %ems | Error: %erad\n", prefix, e,
            time_avg_ms, error_avg_rad);
  }

  fflush(file_report);
  fclose(file_report);
  fflush(file_orbit);
  fclose(file_orbit);
}

int main(int argc, char **argv) {
  // char *end_ptr;

  struct stat st = {0};
  if (stat("./output", &st) == -1) {
    mkdir("./output", 0700);
  }

  report_full_orbits(100, "BISEC_v1", ks_bisection_v1);

  report_full_orbits(100, "HH_PC_NR_NR", ks_hh_pc_nr_nr);
  report_full_orbits(100, "HH_PC_HA_NR", ks_hh_pc_ha_nr);
  report_full_orbits(100, "HH_PC_OG_NR", ks_hh_pc_og131_nr);

  report_full_orbits(100, "EHH_PC_NR", ks_ehh_pc_nr);
  report_full_orbits(100, "EHH_PC_HA", ks_ehh_pc_ha);
  report_full_orbits(100, "EHH_PC_OG", ks_ehh_pc_og131);
}
