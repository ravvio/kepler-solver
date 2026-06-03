#ifndef KS_INCLUDE
#define KS_INCLUDE

#include <math.h>
#include <float.h>

#define TOLERANCE 3e-9
#define MAX_ITERATIONS 1000

typedef double_t (*KS_Solver)(double_t M, double_t e);

double_t ks_bisection_v0(double_t M, double_t e);

double_t ks_newton_v0(double_t M, double_t e);

double_t ks_enrke_v0(double_t M, double_t e);
double_t ks_enrke_v1(double_t M, double_t e);

#endif
