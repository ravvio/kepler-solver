#ifndef KS_INCLUDE
#define KS_INCLUDE

#include <float.h>
#include <math.h>

#define TOLERANCE 3e-15
#define MAX_ITERATIONS 1000

typedef double_t (*KS_Solver)(double_t M, double_t e);

double_t ks_bisection_v1(double_t M, double_t e);

double_t ks_enrke_v1(double_t M, double_t e);
double_t ks_enrke_v2(double_t M, double_t e);

typedef double_t (*KS_DeltaFunction)(double_t M, double_t e, double_t E);
double_t delta_newton(double_t M, double_t e, double_t E);
double_t delta_og131(double_t M, double_t e, double_t E);
double_t delta_halley(double_t M, double_t e, double_t E);

double_t ks_hh_pc_nr_nr(double_t M, double_t e);
double_t ks_hh_pc_ha_nr(double_t M, double_t e);
double_t ks_hh_pc_og131_nr(double_t M, double_t e);

double_t ks_ehh_pc_nr(double_t M, double_t e);
double_t ks_ehh_pc_ha(double_t M, double_t e);
double_t ks_ehh_pc_og131(double_t M, double_t e);

#endif
