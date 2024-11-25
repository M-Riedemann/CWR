#ifndef RNG_H
#define RNG_H

#include "my_numerics.h"


// Typ stochastische DGL
typedef int mlr_sde_func(double, const double[], double[], double[], void*);

// Statistische Methoden: Mittelwert und Standartdabweichung
double mlr_mean(double y[], int dimension);
double mlr_sigma(double y[], int dimension);


// Tuple of 2 gaussian-distributed random numbers with static random-number-generator: Polarmethode
struct mlr_tuple_2 mlr_random_gaussian(gsl_rng* generator);

// Volumen eines Hyperquaders R
double mlr_domain_volume(int D, double R[]);

// 2-Dimensional Integration: Monte-Carlo
double mlr_mc_integrate_2D(gsl_rng* generator, int A(double, double), double a_x, double b_x, double a_y, double b_y, int N, double f(double, double));

// MC-Berechnung der Dichte im Hyperquader R [xi_min, xi_max] mit D Dimensionen für die Dichtefunktion integrand()  
double mlr_mc_integrate(gsl_rng* generator, int D, double integrand(int, double*), double R[], int N);

// Euler-Maruyama Integration stochastischer DGLs
void mlr_eulerMaruyama_step(double t, double delta_t, double y[], mlr_sde_func func, int dimension, void *params);

#endif
