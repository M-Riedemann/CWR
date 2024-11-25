#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>
#include "my_numerics.h"


// Typ stochastische DGL
typedef int mlr_sde_func(double, const double[], double[], double[], void*);


// Statistische Methoden: Mittelwert
double mlr_mean(double y[], int dimension) {
    double mean_sum = 0;
    for (int i = 0; i < dimension; i++) {
        mean_sum += y[i];
    }
    return mean_sum / dimension;
}

// Statistische Methoden: Standartdabweichung
double mlr_sigma(double y[], int dimension) {
    double sigma_sum = 0;
    double y_mean = mlr_mean(y, dimension);
    for (int i = 0; i < dimension; i++) {
        sigma_sum += pow(y[i] - y_mean, 2);
    }
    return sqrt(sigma_sum / dimension);
}



// Tupel aus 2 normalverteile Zuvallsgrößen für gegebenen (GSL_RNG) Zufallsgenerator: Polarmethode
struct mlr_tuple_2 mlr_random_gaussian(gsl_rng* generator) {   
    double u, v, r = 0, m;                                          // Erstellung 2 Zufallszahlen u, v
    while (r > 1 || r == 0){
        u = (gsl_rng_uniform(generator) * 2) - 1;
        v = (gsl_rng_uniform(generator) * 2) - 1;
        r = pow(u, 2) + pow(v, 2);
    }
    m = sqrt(-2*log(r)/r);
    struct mlr_tuple_2 random_2;                                    // Tupel mit 2 normalverteilten Zufallszahlen
    random_2.x1 = u*m;
    random_2.x2 = v*m;
    return random_2;
}  


// Volumen eines Hyperquaders R [xi_min, xi_max] mit D Dimensionen
double mlr_domain_volume(int D, double R[]) {
    double volume_R = 1;
    for (int i_dim = 0; i_dim < D; i_dim++) {
        volume_R *= (R[2*i_dim + 1] - R[2*i_dim]);
    }
    return volume_R;
}


// 2-Dimensionale Integration: MC
double mlr_mc_integrate_2D(gsl_rng* generator, int A(double, double), double a_x, double b_x, double a_y, double b_y, int N, double f(double, double)) {
    double x, y, sum = 0, R_area = (b_x - a_x) * (b_y - a_y);
    for (int i = 0; i < N; i++) {
        x = gsl_rng_uniform(generator) * (b_x - a_x) + a_x;
        y = gsl_rng_uniform(generator) * (b_y - a_y) + a_y;
        sum += f(x, y) * A(x, y);
    }
    return R_area * sum / N;
}


// MC-Berechnung der Dichte im Hyperquader R [xi_min, xi_max] mit D Dimensionen für die Dichtefunktion integrand()  
double mlr_mc_integrate(gsl_rng* generator, int D, double integrand(int, double*), double R[], int N) { 
    double total_volume = mlr_domain_volume(D, R);                      // Gesamtvolumen
    double density_sum = 0;                                         // Laufvariable der Dichtebestimmung
    double *x = (double*) malloc(sizeof(double) * D);               // Erstellung N Zufallsvektoren x[]
    for (int i = 0; i < N; i++) {
        for (int i_dim = 0; i_dim < D; i_dim++) {
            x[i_dim] = gsl_rng_uniform(generator) * (R[2*i_dim + 1] - R[2*i_dim]) + R[2*i_dim];
        }
        density_sum += integrand(D, x);
    }
    free(x);
    return total_volume / N * density_sum;
}


// Euler-Maruyama Integration stochastischer DGLs
void mlr_eulerMaruyama_step(double t, double delta_t, double y[], mlr_sde_func func, int dimension, void *params) {
    double *f = (double*) malloc(sizeof(double) * dimension);       
    double *g = (double*) malloc(sizeof(double) * dimension);
    func(t, y, f, g, params);                                       // Berechnung f und g-Arrays
    static gsl_rng* generator = NULL;
    if (generator == NULL) {
        generator = gsl_rng_alloc(gsl_rng_mt19937);                 // Wahl des Zufallsgenerators: gsl_rng_mt19937
        gsl_rng_set(generator, time(NULL));                      
    }                                                                                
    for (int i = 0; i < dimension; i++) {
        struct mlr_tuple_2 rand_gauss_2 = mlr_random_gaussian(generator);
        double w = rand_gauss_2.x1 * sqrt(delta_t);
        y[i] += f[i] * delta_t + g[i] * w;
    }
    free(f), free(g);
    return;
}