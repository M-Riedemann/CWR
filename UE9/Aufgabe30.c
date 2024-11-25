#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "../my_numerics.h"
#include "../mlr_rand.h"




// zu integrierende Funktion f(x, y)
double f(double x, double y) {
    return pow(x, 2) + 6*x*y + pow(y, 2);
}


// Areal A über das integriert wird
int A(double x, double y) {
    if (pow(x, 2) + pow(y, 2) <= 1) {
        return 1;
    }
    return 0;
}


int main(void) {
    // Abmessung Rechteck R
    double x_a = -3, x_b = 1, y_a = 0, y_b = 4;
    double integral_midpoint, integral_mc;

    // Midpoint-Integration
    FILE* file_midpoint_integration = fopen("A30_Integration_Mittelpunktsregel.csv", "w");
    fprintf(file_midpoint_integration, "Delta x, Integral\n");

    double delta_x = 0.1;
    for (int i = 0; i < 6; i++) {
        delta_x /= 2;
        integral_midpoint = integrate_midpoint_2D(A, x_a, x_b, y_a, y_b, delta_x, f);
        fprintf(file_midpoint_integration, "%g, %g\n", delta_x, integral_midpoint);
        printf("Delta x: %f\n", delta_x);
    }
    fclose(file_midpoint_integration);

    // MC-Integration
    FILE* file_mc_integration = fopen("data/A30_Integration_MC.csv", "w");
    fprintf(file_mc_integration, "N, Integral\n");

    // Wahl des Zufallsgenerators: gsl_rng_mt19937
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL)); 

    int N = 10;
    for (int i = 0; i < 6; i++) {
        N *= 10;
        integral_mc = mlr_mc_integrate_2D(A, x_a, x_b, y_a, y_b, N, f, rng);
        fprintf(file_mc_integration, "%d, %g\n", N, integral_mc);
        printf("Anzahl Stützstellen: %d\n", N);
    }
    fclose(file_mc_integration);

    // Gleichung 4 für f(x, y) = 1 sollte Flächeninhalt des Einheitskreises geben, also Pi / 2
    return 0;
}